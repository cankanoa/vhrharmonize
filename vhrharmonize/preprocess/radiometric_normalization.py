"""Radiometric normalization entrypoints (RRN/ARN)."""

from typing import Any


def spectralmatch_rrn(reference_raster: str, target_raster: str, output_raster: str, **kwargs: Any) -> str:
    """Apply relative radiometric normalization.

    This function is a scaffold that validates dependency presence and defines
    the target callable contract.
    """
    _ = (reference_raster, target_raster, output_raster, kwargs)
    try:
        import spectralmatch  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "spectralmatch is not installed. Install extras with `pip install -e \".[spectralmatch]\"`."
        ) from exc

    raise NotImplementedError(
        "spectralmatch_rrn scaffold added. Implement project-specific RRN strategy."
    )


def apply_rrn(reference_raster: str, target_raster: str, output_raster: str, **kwargs: Any) -> str:
    """Apply relative radiometric normalization (RRN)."""
    return spectralmatch_rrn(reference_raster, target_raster, output_raster, **kwargs)


def apply_arn(input_raster: str, output_raster: str, method: str = "py6s", **kwargs: Any) -> str:
    """Apply absolute radiometric normalization (ARN) scaffold.

    Intended future backends include Py6S-driven atmospheric normalization.
    """
    _ = (input_raster, output_raster, method, kwargs)
    raise NotImplementedError("ARN scaffold added; implement backend-specific normalization strategy.")


__all__ = ["spectralmatch_rrn", "apply_rrn", "apply_arn"]
