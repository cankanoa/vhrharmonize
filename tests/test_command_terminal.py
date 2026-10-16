"""Interactive subprocess prompts must not compete with the dashboard terminal."""
import sys

import pytest

from vhrharmonize.workflow.command_terminal import run_command


class Display:
    def __init__(self, answer):
        self.answer = answer
        self.events = []
        self.writer = None
        self.answered = False

    def begin_command(self, writer):
        self.writer = writer

    def command_output(self, text, *, secret, replace=False):
        self.events.append((text, secret))
        if 'Secret: ' in text or 'Continue? ' in text:
            if not self.answered:
                self.answered = True
                self.writer(self.answer.encode() + b'\n')

    def end_command(self):
        self.writer = None


@pytest.mark.parametrize('tty_prompt', [True, False])
def test_prompt_without_newline_and_progress_are_separate(tty_prompt):
    display = Display('private-answer')
    code = (
        "import getpass; answer = getpass.getpass('Secret: ')" if tty_prompt else
        "answer = input('Continue? ')"
    ) + "; assert answer == 'private-answer'; print('VHR:>f:4:file'); print('4 100%')"
    records = []
    result = run_command([sys.executable, '-c', code], display, output_callback=records.append)
    assert result.returncode == 0
    assert 'VHR:>f:4:file' in records
    assert '4 100%' in records
    assert display.writer is None
    if tty_prompt:
        assert any('Secret: ' in text and secret for text, secret in display.events)
        assert all('private-answer' not in text for text, _ in display.events)
        assert 'private-answer' not in result.stdout


def test_failed_command_restores_input_and_retains_diagnostics():
    import subprocess
    display = Display('')
    with pytest.raises(subprocess.CalledProcessError) as error:
        run_command([sys.executable, '-c', "print('failure'); raise SystemExit(7)"], display)
    assert error.value.returncode == 7
    assert 'failure' in error.value.output
    assert display.writer is None


def test_prompt_is_above_dashboard_and_password_input_is_masked():
    from io import StringIO
    from prompt_toolkit.input import create_pipe_input
    from prompt_toolkit.output import DummyOutput
    from vhrharmonize.workflow.upload_progress import TerminalUploadDisplay
    display = TerminalUploadDisplay(stream=StringIO())
    with create_pipe_input() as input:
        app = display.create_application(input=input, output=DummyOutput())
        sent = []
        display.begin_command(sent.append)
        display.command_output('Secret: ', secret=True)
        assert app.layout.container.children[0].content.children[-1] == display.command_input.window
        display.command_input.buffer.text = 'private-answer'
        display.command_input.buffer.validate_and_handle()
        assert sent == [b'private-answer\n']
        assert display.command_input.buffer.text == ''
        assert not list(display.command_input.buffer.history.get_strings())


def test_direct_controlling_terminal_prompt_routes_through_upload_command_runner():
    from vhrharmonize.slurm import _run_local_command
    from vhrharmonize.workflow.command_terminal import command_display
    display = Display('yes')
    display.app = object()
    code = "import os; fd = os.open('/dev/tty', os.O_RDWR); os.write(fd, b'Continue? '); assert os.read(fd, 100).strip() == b'yes'; os.close(fd)"
    with command_display(display):
        result = _run_local_command([sys.executable, '-c', code])
    assert result.returncode == 0
    assert display.answered
    assert display.writer is None


def test_cancellation_terminates_child_and_releases_prompt():
    class CancelDisplay(Display):
        def command_output(self, *args, **kwargs):
            raise KeyboardInterrupt
    display = CancelDisplay('')
    with pytest.raises(KeyboardInterrupt):
        run_command([sys.executable, '-c', "input('Continue? ')"], display)
    assert display.writer is None


def test_live_dashboard_forwards_keyboard_to_password_prompt(monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from io import StringIO
    from time import monotonic, sleep
    from prompt_toolkit.input import create_pipe_input
    from test_progress_terminal import ScreenOutput
    from vhrharmonize.workflow import progress_terminal as terminal
    from vhrharmonize.workflow.upload_progress import TerminalUploadDisplay, UploadProgress

    class TtyStream(StringIO):
        def isatty(self):
            return True

    with create_pipe_input() as pipe:
        monkeypatch.setattr(terminal, 'create_input', lambda **kw: pipe)
        monkeypatch.setattr(terminal, 'create_output', lambda **kw: ScreenOutput())
        monkeypatch.setattr(terminal.sys, 'stdin', TtyStream())
        monkeypatch.setenv('TERM', 'xterm-256color')
        display = TerminalUploadDisplay(stream=TtyStream())
        display.update(UploadProgress([]).snapshot())
        with display, ThreadPoolExecutor(1) as workers:
            job = workers.submit(run_command, [sys.executable, '-c',
                "import getpass, signal; signal.alarm(5); assert getpass.getpass('Secret: ') == 'private-answer'"], display)
            deadline = monotonic() + 3
            while 'Secret: ' not in display.command_text and monotonic() < deadline:
                sleep(.01)
            assert 'Secret: ' in display.command_text
            assert display.command_secret
            pipe.send_text('private-answer\r')
            assert job.result(timeout=6).returncode == 0
            assert display.command_input.buffer.text == ''
        assert display.error is None
