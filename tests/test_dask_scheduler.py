"""The engine and footprint workers share the plugin scheduler format."""
from unittest.mock import Mock, patch
from types import ModuleType

import pytest

from vhrharmonize.workflow.concurrency import _make_dask_client


@pytest.mark.parametrize("kind", ["file", "address"])
def test_scheduler_connection(kind, tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    client = Mock()
    distributed = ModuleType("dask.distributed")
    distributed.Client = client
    target = "~/scheduler.json" if kind == "file" else " tcp://scheduler:8786 "
    with patch.dict("sys.modules", {"dask.distributed": distributed}):
        assert _make_dask_client([kind, target]) is client.return_value
    if kind == "file":
        client.assert_called_once_with(scheduler_file=str(tmp_path / "scheduler.json"))
    else:
        client.assert_called_once_with("tcp://scheduler:8786")


@pytest.mark.parametrize("value", [None, "address", [], ["file"], ["bad", "target"], ["address", " "]])
def test_invalid_connection(value):
    with pytest.raises(ValueError, match="dask_scheduler"):
        _make_dask_client(value)
