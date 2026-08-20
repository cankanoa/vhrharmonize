"""Expose an interactive POSIX shell through a pseudo-terminal.

The helper is a separate process because creating a PTY by forking directly
inside the multithreaded QGIS process is unsafe.
"""

from __future__ import annotations

import errno
import fcntl
import os
import pty
import pwd
import select
import shutil
import signal
import struct
import sys
import termios


def _write_all(fd: int, data: bytes) -> None:
    while data:
        try:
            written = os.write(fd, data)
        except InterruptedError:
            continue
        data = data[written:]


def _shell_path() -> str:
    candidates = (
        shutil.which("bash"),
        os.environ.get("SHELL"),
        pwd.getpwuid(os.getuid()).pw_shell,
        shutil.which("zsh"),
        shutil.which("sh"),
    )
    for candidate in candidates:
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    raise RuntimeError("No interactive shell executable was found")


def _set_terminal_size(master_fd: int, rows: int = 40, columns: int = 160) -> None:
    size = struct.pack("HHHH", rows, columns, 0, 0)
    fcntl.ioctl(master_fd, termios.TIOCSWINSZ, size)


def _emit_echo_state(master_fd: int, previous: bool | None) -> bool:
    enabled = bool(termios.tcgetattr(master_fd)[3] & termios.ECHO)
    if enabled != previous:
        _write_all(sys.stderr.fileno(), b"ECHO:1\n" if enabled else b"ECHO:0\n")
    return enabled


def main() -> int:
    working_directory = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    child_pid, master_fd = pty.fork()
    if child_pid == 0:
        os.chdir(working_directory)
        shell = _shell_path()
        environment = os.environ.copy()
        environment.setdefault("TERM", "xterm-256color")
        environment.setdefault("COLORTERM", "truecolor")
        arguments = [shell, "--noediting", "-i"] if shell.endswith("bash") else [shell, "-i"]
        # The terminal intentionally starts the executable validated by _shell_path().
        os.execvpe(shell, arguments, environment)  # nosec B606

    _set_terminal_size(master_fd)
    os.set_blocking(sys.stdin.fileno(), False)
    last_echo_state: bool | None = None

    def stop_child(_signum, _frame) -> None:
        try:
            foreground_pid = os.tcgetpgrp(master_fd)
            if foreground_pid > 0:
                os.killpg(foreground_pid, signal.SIGHUP)
        except (OSError, ProcessLookupError):
            pass
        try:
            os.killpg(child_pid, signal.SIGHUP)
        except ProcessLookupError:
            pass

    signal.signal(signal.SIGTERM, stop_child)
    signal.signal(signal.SIGHUP, stop_child)

    try:
        while True:
            waited_pid, wait_status = os.waitpid(child_pid, os.WNOHANG)
            if waited_pid == child_pid:
                return os.waitstatus_to_exitcode(wait_status)

            last_echo_state = _emit_echo_state(master_fd, last_echo_state)
            readable, _, _ = select.select([master_fd, sys.stdin.fileno()], [], [], 0.1)
            if master_fd in readable:
                try:
                    output = os.read(master_fd, 65536)
                except OSError as exc:
                    if exc.errno == errno.EIO:
                        continue
                    raise
                if output:
                    _write_all(sys.stdout.fileno(), output)

            if sys.stdin.fileno() in readable:
                try:
                    incoming = os.read(sys.stdin.fileno(), 65536)
                except BlockingIOError:
                    incoming = b""
                if not incoming:
                    stop_child(signal.SIGHUP, None)
                    return 0
                _write_all(master_fd, incoming)
    finally:
        try:
            os.close(master_fd)
        except OSError:
            pass


if __name__ == "__main__":
    raise SystemExit(main())
