"""Persistent interactive system shell used by the Run tab."""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

from qgis.PyQt.QtCore import QObject, QProcess, QTimer, Qt, pyqtSignal
from qgis.PyQt.QtWidgets import QLineEdit


ANSI_PATTERN = re.compile(
    r"\x1b(?:\][^\x07]*(?:\x07|\x1b\\)|\[[0-?]*[ -/]*[@-~]|[@-_])"
)


class TerminalInput(QLineEdit):
    """Single terminal input with local Up/Down command history."""

    submitted = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._history: list[str] = []
        self._history_index = 0

    def submit(self) -> None:
        text = self.text()
        normal_mode = getattr(QLineEdit, "Normal", None)
        if normal_mode is None:
            normal_mode = QLineEdit.EchoMode.Normal
        if text and self.echoMode() == normal_mode:
            if not self._history or self._history[-1] != text:
                self._history.append(text)
            self._history_index = len(self._history)
        self.clear()
        self.submitted.emit(text)

    def keyPressEvent(self, event) -> None:
        if event.key() in {Qt.Key.Key_Return, Qt.Key.Key_Enter}:
            self.submit()
            event.accept()
            return
        if event.key() == Qt.Key.Key_Up and self._history:
            self._history_index = max(0, self._history_index - 1)
            self.setText(self._history[self._history_index])
            self.end(False)
            return
        if event.key() == Qt.Key.Key_Down and self._history:
            self._history_index = min(len(self._history), self._history_index + 1)
            self.setText(
                "" if self._history_index == len(self._history) else self._history[self._history_index]
            )
            self.end(False)
            return
        super().keyPressEvent(event)


class TerminalSession(QObject):
    """Keep one interactive shell alive and expose its PTY input/output."""

    output_received = pyqtSignal(str)
    input_echo_changed = pyqtSignal(bool)
    state_changed = pyqtSignal(str)

    def __init__(self, working_directory: Path, parent: QObject | None = None):
        super().__init__(parent)
        self.working_directory = Path(working_directory)
        self._process: QProcess | None = None
        self._waiting_writes: list[bytes] = []
        self._shutting_down = False
        self._restart_requested = False

    def start(self) -> None:
        if self._process is not None and not self._not_running(self._process):
            return
        self._shutting_down = False
        process = QProcess(self)
        separate_channels = getattr(QProcess, "SeparateChannels", None)
        if separate_channels is None:
            separate_channels = QProcess.ProcessChannelMode.SeparateChannels
        process.setProcessChannelMode(separate_channels)
        process.setWorkingDirectory(str(self.working_directory))
        process.readyReadStandardOutput.connect(self._read_output)
        process.readyReadStandardError.connect(self._read_control)
        process.started.connect(self._started)
        process.finished.connect(self._finished)
        process.errorOccurred.connect(self._process_error)
        self._process = process

        if os.name == "posix":
            host = Path(__file__).with_name("terminal_host.py")
            process.start(sys.executable, ["-u", str(host), str(self.working_directory)])
        else:
            shell = os.environ.get("COMSPEC", "cmd.exe")
            process.start(shell, ["/Q", "/K"])

    def send_line(self, text: str) -> None:
        self.send_bytes(text.encode() + b"\n")

    def send_bytes(self, payload: bytes) -> None:
        self.start()
        if self._process is None or self._not_running(self._process):
            self._waiting_writes.append(payload)
            return
        starting = getattr(QProcess, "Starting", None)
        if starting is None:
            starting = QProcess.ProcessState.Starting
        if self._process.state() == starting:
            self._waiting_writes.append(payload)
            return
        self._process.write(payload)

    def interrupt(self) -> None:
        self.send_bytes(b"\x03")

    def send_eof(self) -> None:
        self.send_bytes(b"\x04")

    def restart(self) -> None:
        process = self._process
        if process is None or self._not_running(process):
            self.start()
            return
        self._restart_requested = True
        self.state_changed.emit("Restarting shell…")
        process.terminate()
        QTimer.singleShot(1500, lambda: self._kill_if_running(process))

    def shutdown(self) -> None:
        self._shutting_down = True
        self._waiting_writes.clear()
        process = self._process
        if process is None or self._not_running(process):
            self._process = None
            return
        process.terminate()
        if not process.waitForFinished(1000):
            process.kill()
            process.waitForFinished(500)
        self._process = None

    def _started(self) -> None:
        self.state_changed.emit("Interactive shell ready")
        self.input_echo_changed.emit(True)
        if self._process is not None:
            for payload in self._waiting_writes:
                self._process.write(payload)
        self._waiting_writes.clear()

    def _read_output(self) -> None:
        if self._process is None:
            return
        output = bytes(self._process.readAllStandardOutput()).decode(errors="replace")
        output = ANSI_PATTERN.sub("", output).replace("\r", "")
        if output:
            self.output_received.emit(output)

    def _read_control(self) -> None:
        if self._process is None:
            return
        control = bytes(self._process.readAllStandardError()).decode(errors="replace")
        for line in control.splitlines():
            if line == "ECHO:1":
                self.input_echo_changed.emit(True)
            elif line == "ECHO:0":
                self.input_echo_changed.emit(False)
            elif line:
                self.output_received.emit(f"{line}\n")

    def _finished(self, exit_code: int, _exit_status) -> None:
        self._read_output()
        self._read_control()
        self._process = None
        self.input_echo_changed.emit(True)
        if self._shutting_down:
            self.state_changed.emit("Shell stopped")
            return
        if not self._restart_requested:
            self.output_received.emit(f"\n[Shell exited with code {exit_code}; restarting]\n")
        self._restart_requested = False
        self.start()

    def _process_error(self, error) -> None:
        failed_to_start = getattr(QProcess, "FailedToStart", None)
        if failed_to_start is None:
            failed_to_start = QProcess.ProcessError.FailedToStart
        if error == failed_to_start:
            self.state_changed.emit("Shell failed to start")
            self.output_received.emit("The interactive shell could not be started.\n")

    @staticmethod
    def _not_running(process: QProcess) -> bool:
        not_running = getattr(QProcess, "NotRunning", None)
        if not_running is None:
            not_running = QProcess.ProcessState.NotRunning
        return process.state() == not_running

    @classmethod
    def _kill_if_running(cls, process: QProcess) -> None:
        if not cls._not_running(process):
            process.kill()
