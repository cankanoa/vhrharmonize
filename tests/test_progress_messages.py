"""Compatibility with standard Python streams and logging during progress capture."""

from concurrent.futures import ThreadPoolExecutor
from io import StringIO
import logging
import sys
from threading import Barrier

import pytest

from vhrharmonize.io.progress import capture_messages, progress_context


def test_standard_stream_writes_preserve_lines_and_flush(monkeypatch):
    stdout, stderr, messages = StringIO(), StringIO(), []
    monkeypatch.setattr(sys, "stdout", stdout)
    monkeypatch.setattr(sys, "stderr", stderr)
    with progress_context(messages=messages.append), capture_messages():
        print("normal print")
        assert sys.stdout.write("partial ") == 8
        assert sys.stdout.writelines(("line", "\n", "\n")) is None
        sys.stderr.write("first\r\nsecond\r")
        sys.stderr.write("")
        sys.stderr.write("\nthird\n")
        sys.stdout.write("x" * 5000 + "\n")
        sys.stderr.write("tail")
        sys.stderr.flush()
        assert messages[-1] == "tail"
        sys.stdout.write("  ")
        assert sys.stdout.writable() and not sys.stdout.isatty()
        with pytest.raises(TypeError):
            sys.stdout.write(b"bytes")
    assert messages == ["normal print", "partial line", "", "first", "second", "third", "x" * 5000, "tail", "  "]
    assert (sys.stdout, sys.stderr) == (stdout, stderr)
    with capture_messages():
        print("outside a workflow")
    assert stdout.getvalue() == "outside a workflow\n" and stderr.getvalue() == ""


def test_logging_handlers_keep_format_levels_and_file_output_and_restore_on_error(tmp_path, monkeypatch):
    stdout, stderr, messages = StringIO(), StringIO(), []
    monkeypatch.setattr(sys, "stdout", stdout)
    monkeypatch.setattr(sys, "stderr", stderr)
    logger = logging.getLogger("vhr-message-compatibility")
    late_logger = logging.getLogger("vhr-message-created-during-capture")
    dynamic_logger = logging.getLogger("vhr-message-dynamic-stderr")
    console = logging.StreamHandler(stderr)  # Configured before the workflow starts.
    console.setLevel(logging.WARNING)
    console.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    file_handler = logging.FileHandler(tmp_path / "messages.log")
    file_stream = file_handler.stream
    monkeypatch.setattr(logger, "handlers", [console, file_handler])
    monkeypatch.setattr(logger, "propagate", False)
    monkeypatch.setattr(late_logger, "handlers", [])
    monkeypatch.setattr(late_logger, "propagate", False)
    monkeypatch.setattr(dynamic_logger, "handlers", [logging.lastResort])
    monkeypatch.setattr(dynamic_logger, "propagate", False)
    old_level = logger.level
    logger.setLevel(logging.INFO)
    late = None
    try:
        with pytest.raises(RuntimeError, match="processing failed"):
            with progress_context(messages=messages.append), capture_messages():
                logger.info("file only")
                logger.warning("console message")
                with capture_messages():
                    logger.warning("nested capture")
                late = logging.StreamHandler(sys.stdout)
                late_logger.addHandler(late)
                late_logger.warning("created during capture")
                dynamic_logger.warning("dynamic stderr")
                try:
                    raise ValueError("example error")
                except ValueError:
                    logger.exception("details")
                raise RuntimeError("processing failed")
        assert console.stream is stderr and late.stream is stdout
        assert file_handler.stream is file_stream
        assert (sys.stdout, sys.stderr) == (stdout, stderr)
        assert messages[:3] == ["WARNING: console message", "WARNING: nested capture", "created during capture"]
        assert "dynamic stderr" in messages
        assert "ERROR: details" in messages and "ValueError: example error" in messages
        assert "file only" not in messages
        assert stdout.getvalue() == stderr.getvalue() == ""
        assert (tmp_path / "messages.log").read_text().count("console message") == 1
        assert "file only" in (tmp_path / "messages.log").read_text()
        logger.warning("after workflow")
        assert stderr.getvalue() == "WARNING: after workflow\n"
    finally:
        logger.setLevel(old_level)
        console.close()
        file_handler.close()
        if late is not None:
            late.close()


def test_concurrent_workers_keep_partial_lines_and_logging_in_their_own_context(monkeypatch):
    stdout, stderr = StringIO(), StringIO()
    monkeypatch.setattr(sys, "stdout", stdout)
    monkeypatch.setattr(sys, "stderr", stderr)
    logger = logging.getLogger("vhr-concurrent-messages")
    handler = logging.StreamHandler(stderr)
    monkeypatch.setattr(logger, "handlers", [handler])
    monkeypatch.setattr(logger, "propagate", False)
    barrier = Barrier(2)

    def worker(name):
        messages = []
        with progress_context(messages=messages.append), capture_messages():
            sys.stdout.write(name)
            barrier.wait(timeout=5)
            sys.stdout.writelines((" complete", "\n"))
            logger.warning("%s logged", name)
        return messages

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            first, second = pool.submit(worker, "first"), pool.submit(worker, "second")
            assert first.result(timeout=5) == ["first complete", "first logged"]
            assert second.result(timeout=5) == ["second complete", "second logged"]
        assert handler.stream is stderr and (sys.stdout, sys.stderr) == (stdout, stderr)
        assert stdout.getvalue() == stderr.getvalue() == ""
    finally:
        handler.close()
