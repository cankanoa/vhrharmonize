"""Provider-agnostic scene pipeline contracts."""

from dataclasses import dataclass
from typing import Any, Dict


@dataclass
class ScenePipelineResult:
    """Standard return object for provider scene pipelines."""

    scene_id: str
    outputs: Dict[str, Any]
