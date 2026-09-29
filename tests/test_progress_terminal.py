"""Exercise the terminal UI with real key/mouse input and screen rendering."""

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

    def get_size(self):
        return self.size

    def enter_alternate_screen(self):
        self.entered = True

    def quit_alternate_screen(self):
        self.restored = True


def test_scrolling_pins_progress_and_keeps_receiving_updates():
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
            first = await frame_when(lambda f: any("message 039" in line for line in f))
            workflow_y = next(i for i, line in enumerate(first) if "Workflow progress" in line)
            active_y = next(i for i, line in enumerate(first) if "Active operation" in line)
            assert "working" in "\n".join(first) and output.entered
            assert all(line.startswith("message ") for line in first[:workflow_y])
            # PageUp affects only the message viewport.
            pipe.send_text("\x1b[5~")
            paused = await frame_when(lambda f: f[0] != first[0])
            assert paused[workflow_y:] == first[workflow_y:]
            assert "message 039" not in "\n".join(paused)
            # New snapshots keep the paused viewport in place while counts update.
            newer = snapshot([f"message {i:03d}" for i in range(41)])
            newer["rows"][0]["done"] = newer["total"]["done"] = 4
            newer["rows"][0]["percentages"]["done"] = newer["total"]["percentages"]["done"] = 100 * 4 / 264
            display.update(newer)
            updated = await frame_when(lambda f: "4(2%)" in "\n".join(f))
            assert updated[:workflow_y] == paused[:workflow_y]
            assert "Active operation" in updated[active_y]
            # End resumes following, and SGR mouse-wheel input scrolls messages.
            pipe.send_text("\x1b[F")
            following = await frame_when(lambda f: "message 040" in "\n".join(f))
            pipe.send_text("\x1b[<64;4;3M")
            wheel = await frame_when(lambda f: f[0] != following[0])
            assert wheel[workflow_y:] == updated[workflow_y:]
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
        assert output.restored

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
    assert output.entered and output.restored
    assert stream.getvalue().count("Workflow progress") == 1


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
    assert len(display.messages) == 1001 and display.follow


def test_static_colors_plain_output_and_old_snapshots(monkeypatch):
    data = snapshot(["scene [red] written"])
    data.pop("message_history")  # Existing version 2 snapshots remain usable.
    value = render_progress(validate_progress_snapshot(data))
    assert "scene [red] written" in to_plain_text(value)  # Never interpret log markup.
    lines = to_plain_text(value).splitlines()
    assert lines[0].startswith("scene [red] written")
    assert "Recent messages" not in to_plain_text(value)
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
    assert "Workflow progress" in output and "Loaded" in output and "scene" in output
    assert "\x1b" not in output and output.isascii()


@pytest.mark.parametrize("history", [{"sequence": -1, "messages": []},
                                      {"sequence": 0, "messages": ["extra"]},
                                      {"sequence": 1.5, "messages": []}])
def test_reject_invalid_history_cursor(history):
    data = deepcopy(snapshot())
    data["message_history"] = history
    with pytest.raises(ValueError, match="Malformed"):
        validate_progress_snapshot(data)
