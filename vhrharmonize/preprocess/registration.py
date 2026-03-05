"""Image registration utilities."""

from typing import Any, Optional, Sequence


def _require_itk():
    try:
        import itk
    except ImportError as exc:
        raise RuntimeError(
            "itk-elastix is not installed. Install extras with `pip install -e \".[elastix]\"`."
        ) from exc
    return itk


def estimate_elastix_transform(
    fixed_image_path: str,
    moving_image_path: str,
    *,
    parameter_map: str = "rigid",
    parameter_file_paths: Optional[Sequence[str]] = None,
    fixed_mask_path: Optional[str] = None,
    moving_mask_path: Optional[str] = None,
    log_to_console: bool = False,
) -> Any:
    """Estimate transform parameters that map a moving image to a fixed image.

    Args:
        fixed_image_path: Path to the fixed/reference image.
        moving_image_path: Path to the moving image to be aligned.
        parameter_map: Default elastix parameter map name if no parameter files are provided.
            Typical values include ``"rigid"``, ``"affine"``, and ``"bspline"``.
        parameter_file_paths: Optional elastix parameter file path(s). If provided,
            these are used instead of ``parameter_map``.
        fixed_mask_path: Optional fixed-image mask path.
        moving_mask_path: Optional moving-image mask path.
        log_to_console: If ``True``, emit elastix logs to stdout.

    Returns:
        ITK transform parameter object produced by elastix.

    Raises:
        RuntimeError: If ``itk-elastix`` is not installed.
    """
    itk = _require_itk()

    fixed = itk.imread(fixed_image_path, itk.F)
    moving = itk.imread(moving_image_path, itk.F)

    parameter_object = itk.ParameterObject.New()
    if parameter_file_paths:
        for path in parameter_file_paths:
            parameter_object.AddParameterMap(parameter_object.ReadParameterFile(path))
    else:
        parameter_object.AddParameterMap(parameter_object.GetDefaultParameterMap(parameter_map))

    kwargs = {
        "fixed_image": fixed,
        "moving_image": moving,
        "parameter_object": parameter_object,
        "log_to_console": log_to_console,
    }

    if fixed_mask_path:
        kwargs["fixed_mask"] = itk.imread(fixed_mask_path, itk.UC)
    if moving_mask_path:
        kwargs["moving_mask"] = itk.imread(moving_mask_path, itk.UC)

    _, transform_parameter_object = itk.elastix_registration_method(**kwargs)
    return transform_parameter_object


def apply_elastix_transform(
    moving_image_path: str,
    output_image_path: str,
    transform_parameter_object: Any,
    *,
    log_to_console: bool = False,
) -> str:
    """Apply a precomputed elastix transform to an image and write the result.

    Args:
        moving_image_path: Path to the moving image to warp.
        output_image_path: Output path for the transformed image.
        transform_parameter_object: Transform parameter object from
            :func:`estimate_elastix_transform`.
        log_to_console: If ``True``, emit transformix logs to stdout.

    Returns:
        The written ``output_image_path``.

    Raises:
        RuntimeError: If ``itk-elastix`` is not installed.
    """
    itk = _require_itk()
    moving = itk.imread(moving_image_path, itk.F)
    transformed = itk.transformix_filter(
        moving,
        transform_parameter_object=transform_parameter_object,
        log_to_console=log_to_console,
    )
    itk.imwrite(transformed, output_image_path)
    return output_image_path


def run_elastix_registration(
    fixed_image_path: str,
    moving_image_path: str,
    output_image_path: str,
    *,
    parameter_map: str = "rigid",
    parameter_file_paths: Optional[Sequence[str]] = None,
    fixed_mask_path: Optional[str] = None,
    moving_mask_path: Optional[str] = None,
    log_to_console: bool = False,
) -> str:
    """Register a moving image to a fixed image using itk-elastix.

    Args:
        fixed_image_path: Path to the fixed/reference image.
        moving_image_path: Path to the moving image that will be warped.
        output_image_path: Path where the registered image will be written.
        parameter_map: Default elastix parameter map name when no parameter files are provided.
            Common values: ``"rigid"``, ``"affine"``, ``"bspline"``.
        parameter_file_paths: Optional parameter file path(s). When provided, these are used
            instead of ``parameter_map``.
        fixed_mask_path: Optional fixed-image mask path.
        moving_mask_path: Optional moving-image mask path.
        log_to_console: Whether elastix should print logs to stdout.

    Returns:
        The written ``output_image_path``.
    """
    transform_parameter_object = estimate_elastix_transform(
        fixed_image_path=fixed_image_path,
        moving_image_path=moving_image_path,
        parameter_map=parameter_map,
        parameter_file_paths=parameter_file_paths,
        fixed_mask_path=fixed_mask_path,
        moving_mask_path=moving_mask_path,
        log_to_console=log_to_console,
    )
    apply_elastix_transform(
        moving_image_path=moving_image_path,
        output_image_path=output_image_path,
        transform_parameter_object=transform_parameter_object,
        log_to_console=log_to_console,
    )
    return output_image_path


__all__ = [
    "estimate_elastix_transform",
    "apply_elastix_transform",
    "run_elastix_registration",
]
