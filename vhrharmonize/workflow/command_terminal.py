"""Give child commands a private terminal while the dashboard owns the real one."""
from contextlib import contextmanager
from contextvars import ContextVar

_current = ContextVar('command_display', default=None)


@contextmanager
def command_display(display):
    token = _current.set(display if display is not None and display.app is not None else None)
    try:
        yield
    finally:
        _current.reset(token)


def current_display():
    return _current.get()


def run_command(command, display, *, check=True, output_callback=None):
    import codecs
    from collections import deque
    import errno
    import os
    import pty
    import selectors
    import subprocess
    import sys
    import termios
    from threading import Event

    submitted = Event()
    master, slave = pty.openpty()
    process = None
    output = deque(maxlen=50) if output_callback else []
    pending = ''
    decoder = codecs.getincrementaldecoder('utf-8')(errors='replace')
    tty_decoder = codecs.getincrementaldecoder('utf-8')(errors='replace')
    try:
        # A fresh interpreter acquires the controlling terminal before exec.
        # Avoid preexec_fn: the live dashboard already runs another thread.
        process = subprocess.Popen([sys.executable, __file__, *command], stdin=slave,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        os.close(slave)
        slave = None
        def send(data):
            submitted.set()
            os.write(master, data)
        display.begin_command(send)
        with selectors.DefaultSelector() as selector:
            selector.register(master, selectors.EVENT_READ, 'tty')
            selector.register(process.stdout, selectors.EVENT_READ, 'output')
            while selector.get_map():
                for key, _ in selector.select(.1):
                    if submitted.is_set():
                        pending = ''
                        submitted.clear()
                    try:
                        data = os.read(key.fd, 65536)
                    except OSError as exc:
                        if exc.errno != errno.EIO:
                            raise
                        data = b''  # PTY hangup on Linux.
                    if not data:
                        selector.unregister(key.fileobj)
                        continue
                    if key.data == 'tty':
                        echo = bool(termios.tcgetattr(master)[3] & termios.ECHO)
                        display.command_output(tty_decoder.decode(data), secret=not echo)
                    else:
                        text = decoder.decode(data)
                        if output_callback:
                            pending += text
                            import re
                            records = re.split(r'[\r\n]', pending)
                            pending = records.pop()
                            for line in records:
                                if line:
                                    output.append(line + '\n')
                                    output_callback(line)
                                    if not line.startswith('VHR:') and not re.match(r'^\s*[\d,]+\s+(?:\d+|nan)%', line):
                                        display.command_output(line + '\n', secret=not bool(termios.tcgetattr(master)[3] & termios.ECHO))
                            if pending and not pending.startswith('VHR:') and not re.match(r'^\s*[\d,]+\s+\d+%', pending):
                                display.command_output(pending, secret=not bool(termios.tcgetattr(master)[3] & termios.ECHO), replace=True)
                        else:
                            output.append(text)
                            display.command_output(text, secret=not bool(termios.tcgetattr(master)[3] & termios.ECHO))
                if process.poll() is not None and not selector.select(0):
                    break
        if pending:
            output.append(pending)
            output_callback(pending)
        code = process.wait()
        result = ''.join(output)
        if check and code:
            raise subprocess.CalledProcessError(code, command, output=result)
        return subprocess.CompletedProcess(command, code, stdout=result, stderr='')
    finally:
        display.end_command()
        if process is not None:
            if process.poll() is None:
                import signal
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            process.stdout.close()
        if slave is not None:
            os.close(slave)
        os.close(master)


if __name__ == '__main__':
    import fcntl
    import os
    import sys
    import termios
    fcntl.ioctl(0, termios.TIOCSCTTY, 0)
    os.execvp(sys.argv[1], sys.argv[1:])
