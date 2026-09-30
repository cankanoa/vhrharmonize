"""Exercise native terminal output and progress-box rendering."""

import asyncio
from copy import deepcopy
from io import BytesIO, StringIO, TextIOWrapper
from threading import current_thread

import pytest
from prompt_toolkit.data_structures import Size
from prompt_toolkit.formatted_text import to_plain_text
from prompt_toolkit.input import create_pipe_input
from prompt_toolkit.output import DummyOutput
from prompt_toolkit.utils import get_cwidth

from vhrharmonize.progress import render_progress, validate_progress_snapshot
from vhrharmonize.workflow import progress_terminal as terminal
from vhrharmonize.workflow.progress import ProgressState


def snapshot(messages=(), sequence=None):
    def row(name, done, active):
        counts = dict(unused=52, reused=192, done=done, run=20, all=264)
        return dict(name=name, **counts, percentages={k: 100 * v / 264 for k, v in counts.items()},
                    fraction_done=done / 20, active=active, pending=False, worker_progress=True,
                    status="running", eta_seconds=None)

    alignment = row("alignment", 3, 1)
    return dict(version=2, run_id="test", job_id=None, updated_at="2026-09-29T00:00:00Z",
                status="running", total=row("total", 3, 1), rows=[alignment],
                active=[dict(task_id="tile", step="alignment", scene="scene_P004",
                             stats=dict(n=120, total=200, prefix="matching tiles", unit="tiles",
                                        elapsed=60, rate=2), eta_seconds=40)],
                messages=list(messages)[-5:],
                message_history=dict(sequence=len(messages) if sequence is None else sequence,
                                     messages=list(messages)))


class ScreenOutput(DummyOutput):
    size = Size(rows=24, columns=120)

    def __init__(self):
        self.entered = self.restored = False
        self.written = []
        self.mouse_enabled = False
        self.cursor_visible = True

    def get_size(self):
        return self.size

    def enter_alternate_screen(self):
        self.entered = True

    def quit_alternate_screen(self):
        self.restored = True

    def write(self, text):
        self.written.append(text)

    def enable_mouse_support(self):
        self.mouse_enabled = True

    def hide_cursor(self):
        self.cursor_visible = False

    def show_cursor(self):
        self.cursor_visible = True


def test_native_logs_print_once_above_progress_and_updates_resize():
    async def exercise(pipe):
        output = ScreenOutput()
        display = terminal.TerminalProgressDisplay(stream=StringIO())
        display.update(snapshot([f"message {i:03d}" for i in range(40)]))
        app = display.create_application(input=pipe, output=output)
        frames = asyncio.Queue()

        def rendered(app):
            screen = app.renderer.last_rendered_screen
            if screen is not None:
                frames.put_nowait(["".join(screen.data_buffer[y][x].char for x in range(output.size.columns))
                                   for y in range(output.size.rows)])

        app.after_render += rendered

        async def frame_when(predicate):
            async def receive():
                while True:
                    frame = await frames.get()
                    if predicate(frame):
                        return frame
            return await asyncio.wait_for(receive(), 5)

        task = asyncio.create_task(app.run_async(handle_sigint=False, set_exception_handler=False))
        try:
            first = await frame_when(lambda f: terminal.TITLE in "\n".join(f)
                                     and "message 039" in "".join(output.written))
            workflow_y = next(i for i, line in enumerate(first) if terminal.TITLE in line)
            active_y = next(i for i, line in enumerate(first) if line.startswith("├")) + 1
            assert "working" in "\n".join(first)
            assert not output.entered and not output.mouse_enabled
            assert not app.full_screen and not app.mouse_support()
            assert "message 039" not in "\n".join(first)  # Logs belong to native scrollback.
            # New snapshots append only the new log line and update the same box.
            newer = snapshot([f"message {i:03d}" for i in range(41)])
            newer["rows"][0]["done"] = newer["total"]["done"] = 4
            newer["rows"][0]["percentages"]["done"] = newer["total"]["percentages"]["done"] = 100 * 4 / 264
            display.update(newer)
            updated = await frame_when(lambda f: "4(2%)" in "\n".join(f)
                                       and "message 040" in "".join(output.written))
            assert terminal.TITLE in updated[workflow_y]
            assert "ID" in updated[active_y]
            text = "".join(output.written)
            assert all(text.count(f"message {i:03d}") == 1 for i in range(41))
            # Resizing does not let panels overflow or hide the overall progress.
            for rows, columns in ((18, 80), (10, 60), (30, 140)):
                output.size = Size(rows=rows, columns=columns)
                app.invalidate()
                resized = await frame_when(lambda f: len(f) == rows and "total" in "\n".join(f))
                assert all(get_cwidth(line) <= columns for line in resized)
                assert "Window too small" not in "\n".join(resized)
            crowded = deepcopy(newer)
            crowded["rows"] = [{**newer["rows"][0], "name": f"step_{i}"} for i in range(30)]
            crowded["active"] = newer["active"] * 12
            display.update(crowded)
            for rows in (24, 12):
                output.size = Size(rows=rows, columns=120)
                app.invalidate()
                resized = await frame_when(lambda f: len(f) == rows and "total" in "\n".join(f))
                assert "Window too small" not in "\n".join(resized)
                if rows == 24:
                    assert "more steps" in "\n".join(resized)
                    assert "other active operations" in "\n".join(resized)
        finally:
            if not app.is_done:
                app.exit()
            await asyncio.wait_for(task, 5)
        assert not output.entered and output.cursor_visible

    with create_pipe_input() as pipe:
        asyncio.run(exercise(pipe))


@pytest.mark.parametrize("failure", [False, True])
def test_ui_thread_restores_terminal_without_moving_processing(monkeypatch, failure):
    class Tty(StringIO):
        def isatty(self):
            return True

    output, stream = ScreenOutput(), Tty()
    monkeypatch.setenv("TERM", "xterm")
    monkeypatch.setattr(terminal.sys, "stdin", Tty())
    monkeypatch.setattr(terminal, "create_output", lambda **kwargs: output)
    monkeypatch.setattr(terminal, "print_formatted_text",
                        lambda value, **kwargs: stream.write(to_plain_text(value)))
    caller = current_thread()
    with create_pipe_input() as pipe:
        monkeypatch.setattr(terminal, "create_input", lambda **kwargs: pipe)
        display = terminal.TerminalProgressDisplay(stream=stream)
        try:
            with display:
                assert current_thread() is caller
                display.update(snapshot(["processing scene"]))
                if failure:
                    raise RuntimeError("backend failed")
        except RuntimeError as exc:
            assert failure and str(exc) == "backend failed"
    assert display.error is None
    assert not display.thread.is_alive()
    assert not output.entered and output.cursor_visible
    assert "".join(output.written).count("processing scene") == 1
    assert stream.getvalue().count(terminal.TITLE) == 1
    assert "processing scene" not in stream.getvalue()  # Do not repeat logs in the final box.


def test_message_history_keeps_bursts_and_repeats_and_is_bounded():
    state = ProgressState()
    for i in range(1050):
        state.handle({"kind": "message", "text": f"line {i}"})
    assert len(state.messages) == 1000 and state.messages.sequence == 1050
    display = terminal.TerminalProgressDisplay(stream=StringIO())
    data = snapshot(list(state.messages), state.messages.sequence)
    display.update(data)
    assert "no longer available" in display.messages[0]
    assert list(display.messages)[1:] == list(state.messages)
    for _ in range(12):
        state.handle({"kind": "message", "text": "same message"})
    newer = snapshot(list(state.messages), state.messages.sequence)
    display.update(newer)
    display.update(newer)
    assert list(display.messages)[-12:] == ["same message"] * 12
    assert len(display.messages) == 1013
    display.update(data)  # An older polled snapshot must not move the message cursor back.
    display.update(newer)
    assert len(display.messages) == 1013
    newer["run_id"] = "another run"
    display.update(newer)
    assert len(display.messages) == 1001


def test_static_colors_plain_output_and_old_snapshots(monkeypatch):
    data = snapshot(["scene [red] written"])
    data.pop("message_history")  # Existing version 2 snapshots remain usable.
    value = render_progress(validate_progress_snapshot(data))
    assert "scene [red] written" in to_plain_text(value)  # Never interpret log markup.
    lines = to_plain_text(value).splitlines()
    assert lines[0].startswith("scene [red] written")
    assert "Recent messages" not in to_plain_text(value)
    assert sum(line.startswith("╭") for line in lines) == 1
    assert sum(line.startswith("├") for line in lines) == 1
    assert sum(line.startswith("╰") for line in lines) == 1
    assert to_plain_text(value).count(terminal.TITLE) == 1
    assert "Active operation" not in to_plain_text(value)
    assert [line.strip("│ ").split() for line in lines if line.strip("│ ").startswith("Step ")] == [
        ["Step", "Unused", "Loaded", "Done", "Run", "All", "Active", "ETA", "Progress"],
        ["Step", "ID", "Status", "Elapsed", "ETA", "Progress"],
    ]
    for name, color in (("Unused", "#4b87b9"), ("Loaded", "#956ac3"),
                        ("Done", "#298845"), ("Run", "#858585")):
        assert any(color in style and text == name for style, text in value)
        assert any(color in style and "▬" in text for style, text in value)
    assert "TBD" in to_plain_text(value) and "~0m 40s" in to_plain_text(value)
    stream = TextIOWrapper(BytesIO(), encoding="ascii")
    monkeypatch.setattr(terminal.TerminalProgressDisplay, "create_application",
                        lambda *args, **kwargs: pytest.fail("Noninteractive output started an app"))
    data["messages"] = ["scene 海岸 written"]
    with terminal.TerminalProgressDisplay(stream=stream) as display:
        display.update(data)
    stream.seek(0)
    output = stream.read()
    assert terminal.TITLE in output and "Loaded" in output and "scene" in output
    assert "\x1b" not in output and output.isascii()


@pytest.mark.parametrize("history", [{"sequence": -1, "messages": []},
                                      {"sequence": 0, "messages": ["extra"]},
                                      {"sequence": 1.5, "messages": []}])
def test_reject_invalid_history_cursor(history):
    data = deepcopy(snapshot())
    data["message_history"] = history
    with pytest.raises(ValueError, match="Malformed"):
        validate_progress_snapshot(data)
