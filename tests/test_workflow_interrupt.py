"""Ctrl+C must wake a blocked workflow and stop its owned workers."""

import os
import signal
import subprocess
import sys
import textwrap

import pytest


@pytest.mark.skipif(not hasattr(signal, "pthread_kill"), reason="POSIX terminal interrupt")
@pytest.mark.parametrize("mode", ["main_wait", "scene_pool", "image_pool", "stubborn_worker", "shutdown_wait"])
def test_dashboard_ctrl_c_wakes_wait_and_restores_terminal(tmp_path, mode):
    script = tmp_path / "interrupt_run.py"
    script.write_text(textwrap.dedent('''
        from concurrent.futures import ProcessPoolExecutor
        from contextlib import nullcontext
        from io import StringIO
        from multiprocessing import get_context
        import os
        from pathlib import Path
        import signal
        import sys
        from threading import Event, Thread
        from time import monotonic, sleep

        from prompt_toolkit.input import create_pipe_input
        from prompt_toolkit.output import DummyOutput
        from vhrharmonize.workflow import progress_terminal as terminal
        from vhrharmonize.workflow.concurrency import image_executor, interruptible_pool

        def block(marker, stubborn):
            if stubborn:
                signal.signal(signal.SIGTERM, signal.SIG_IGN)
            Path(marker).write_text(str(os.getpid()))
            Event().wait(60)

        class Tty(StringIO):
            def isatty(self):
                return True

        class Output(DummyOutput):
            visible = True
            def hide_cursor(self):
                self.visible = False
            def show_cursor(self):
                self.visible = True

        if __name__ == "__main__":
            mode, marker = sys.argv[1:]
            output = Output()
            os.environ["TERM"] = "xterm"
            terminal.sys.stdin = Tty()
            terminal.create_output = lambda **kw: output
            with create_pipe_input() as pipe:
                terminal.create_input = lambda **kw: pipe
                display = terminal.TerminalProgressDisplay(stream=Tty())
                def cancel():
                    deadline = monotonic() + 10
                    while not Path(marker).exists() and monotonic() < deadline:
                        sleep(.01)
                    pipe.send_text("\\x03")
                sender = Thread(target=cancel, daemon=True)
                try:
                    with display:
                        if mode == "main_wait":
                            pool = nullcontext(None)
                        elif mode == "image_pool":
                            pool = image_executor(2, "process_pool")
                        else:
                            pool = interruptible_pool(ProcessPoolExecutor(
                                max_workers=2, mp_context=get_context("spawn")))
                        with pool as executor:
                            if executor is None:
                                Path(marker).write_text("main")
                            else:
                                future = executor.submit(block, marker, mode == "stubborn_worker")
                            sender.start()
                            if executor is None:
                                Event().wait(60)
                            elif mode != "shutdown_wait":
                                future.result(timeout=60)
                except KeyboardInterrupt:
                    sender.join(2)
                    assert not display.thread.is_alive()
                    assert output.visible
                    assert display.error is None
                    if mode != "main_wait":
                        pid = int(Path(marker).read_text())
                        try:
                            os.kill(pid, 0)
                        except ProcessLookupError:
                            pass
                        else:
                            raise AssertionError("Cancelled worker is still running")
                    print("INTERRUPTED: terminal restored and workers stopped", flush=True)
                    sys.exit(130)
                raise AssertionError("Ctrl+C did not interrupt the workflow")
    '''))
    environment = {**os.environ, "PYTHONPATH": os.pathsep.join(map(str, sys.path))}
    # A separate session keeps test signals away from pytest and makes timeout
    # cleanup include only this test's subprocess and its spawned workers.
    process = subprocess.Popen([sys.executable, str(script), mode, str(tmp_path / "ready")],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                               env=environment, cwd=tmp_path, start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=20)
        assert process.returncode == 130, stdout + stderr
        assert "INTERRUPTED: terminal restored and workers stopped" in stdout
    finally:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=5)
