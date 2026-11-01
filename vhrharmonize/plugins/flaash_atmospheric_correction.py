"""ENVI FLAASH atmospheric correction."""

from __future__ import annotations

import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from shapely.geometry.base import BaseGeometry

from vhrharmonize.io.geospatial import get_image_percentile_value
from vhrharmonize.io.logging import _log, _logged_operation
from vhrharmonize.io.progress import current_callback, progress, reports_progress
from vhrharmonize.io.progress_transport import local_worker_progress
from vhrharmonize.io.workflow_utils import remove_output_files
from .base import FunctionPlugin


@dataclass(frozen=True)
class FLAASHRunResult:
    """Summary of a shared FLAASH run."""

    output_raster: str
    params: Dict[str, Any]
    params_output_path: str


FLAASH_ALLOWED_PARAMS = {
    "SENSOR_TYPE",
    "INPUT_SCALE",
    "OUTPUT_SCALE",
    "CALIBRATION_FILE",
    "CALIBRATION_FORMAT",
    "CALIBRATION_UNITS",
    "LAT_LONG",
    "SENSOR_ALTITUDE",
    "DATE_TIME",
    "USE_ADJACENCY",
    "DEFAULT_VISIBILITY",
    "USE_POLISHING",
    "POLISHING_RESOLUTION",
    "SENSOR_AUTOCALIBRATION",
    "SENSOR_CAL_PRECISION",
    "SENSOR_CAL_FEATURE_LIST",
    "GROUND_ELEVATION",
    "SOLAR_AZIMUTH",
    "SOLAR_ZENITH",
    "LOS_AZIMUTH",
    "LOS_ZENITH",
    "IFOV",
    "MODTRAN_ATM",
    "MODTRAN_AER",
    "MODTRAN_RES",
    "MODTRAN_MSCAT",
    "CO2_MIXING",
    "WATER_ABS_CHOICE",
    "WATER_MULT",
    "WATER_VAPOR_PRESET",
    "USE_AEROSOL",
    "AEROSOL_SCALE_HT",
    "AER_BAND_RATIO",
    "AER_BAND_WAVL",
    "AER_REFERENCE_VALUE",
    "AER_REFERENCE_PIXEL",
    "AER_BANDLOW_WAVL",
    "AER_BANDLOW_MAXREFL",
    "AER_BANDHIGH_WAVL",
    "AER_BANDHIGH_MAXREFL",
    "INPUT_RASTER",
    "OUTPUT_RASTER_URI",
    "CLOUD_RASTER_URI",
    "WATER_RASTER_URI",
}


def _init_envi_engine(envi_engine_path: str) -> Any:
    """Initialize the ENVI task engine used by FLAASH."""

    import envipyengine.config
    from envipyengine import Engine

    envipyengine.config.set("engine", envi_engine_path)
    envi_engine = Engine("ENVI")
    envi_engine.tasks()
    return envi_engine


def _validate_flaash_params(flaash_params: Dict[str, Any]) -> Dict[str, Any]:
    """Validate FLAASH params against the supported parameter list."""
    unknown_params = sorted(set(flaash_params) - FLAASH_ALLOWED_PARAMS)
    if unknown_params:
        raise ValueError("Unsupported FLAASH parameter(s): " + ", ".join(unknown_params))
    return flaash_params


def _build_flaash_kwargs_from_variables(
    input_raster: str,
    dem_file_path: str,
    footprint_geometry: BaseGeometry,
    variables: Any,
    output_raster: str,
    *,
    dem_ground_percentile: float,
    modtran_atm: str,
    modtran_aer: str,
    use_aerosol: str,
    default_visibility: Optional[float],
    custom_params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build validated shared FLAASH parameters.
    Args:
        input_raster: Input raster path.
        dem_file_path: DEM raster path.
        footprint_geometry: Shapely footprint in EPSG:4326 for DEM sampling.
        variables: Mapped variables JSON object.
        output_raster: Output raster path.
        dem_ground_percentile: DEM percentile used for ground elevation estimation.
        modtran_atm: MODTRAN atmosphere profile name.
        modtran_aer: MODTRAN aerosol profile name.
        use_aerosol: FLAASH aerosol handling mode.
        default_visibility: Optional default visibility override.
        custom_params: Optional custom FLAASH parameter overrides.
    Returns:
        Validated FLAASH parameter dictionary.
    """

    ground_elevation_m = get_image_percentile_value(
        dem_file_path,
        percentile=dem_ground_percentile,
        mask=footprint_geometry,
    )
    flaash_params = {
        "INPUT_RASTER": {"url": input_raster, "factory": "URLRaster"},
        "MODTRAN_ATM": modtran_atm,
        "MODTRAN_AER": modtran_aer,
        "MODTRAN_RES": 5.0,
        "MODTRAN_MSCAT": "DISORT",
        "USE_AEROSOL": use_aerosol,
        "DEFAULT_VISIBILITY": default_visibility,
        "AER_BAND_RATIO": 0.5,
        "AER_BANDLOW_WAVL": 425,
        "AER_BANDHIGH_WAVL": 660,
        "AER_BANDHIGH_MAXREFL": 0.2,
        "GROUND_ELEVATION": ground_elevation_m / 1000.0,
        "SOLAR_AZIMUTH": variables["solar_azimuth"],
        "SOLAR_ZENITH": variables["solar_zenith"],
        "LOS_AZIMUTH": variables["line_of_sight_azimuth"],
        "LOS_ZENITH": variables["line_of_sight_zenith"],
        "OUTPUT_RASTER_URI": output_raster,
    }
    if custom_params:
        flaash_params.update(custom_params)
    flaash_params = {key: value for key, value in flaash_params.items() if value is not None}
    return _validate_flaash_params(flaash_params)


def _wsl_path_to_windows_for_envi(path: str) -> str:
    """Convert `/mnt/<drive>/...` WSL paths into Windows drive paths for ENVI."""
    match = re.match(r"^/mnt/([a-zA-Z])/(.*)$", path)
    if not match:
        raise ValueError(f"Path is not Windows-drive-backed via /mnt/<drive>/: {path}")
    drive = match.group(1).upper()
    rest = match.group(2).replace("/", "\\")
    return f"{drive}:\\{rest}"


def _convert_flaash_params_paths_for_windows(flaash_params: Dict[str, Any]) -> Dict[str, Any]:
    """Convert FLAASH path params to Windows form for ENVI-on-Windows execution."""
    converted = dict(flaash_params)
    input_raster = converted.get("INPUT_RASTER")
    if isinstance(input_raster, dict) and "url" in input_raster:
        updated_input = dict(input_raster)
        updated_input["url"] = _wsl_path_to_windows_for_envi(updated_input["url"])
        converted["INPUT_RASTER"] = updated_input

    for key in ("OUTPUT_RASTER_URI", "CLOUD_RASTER_URI", "WATER_RASTER_URI"):
        if converted.get(key):
            converted[key] = _wsl_path_to_windows_for_envi(converted[key])
    return converted


def _execute_flaash_task(
    flaash_params: Dict[str, Any],
    output_params_path: str,
    envi_engine: Any,
    output_image_path_to_delete: Optional[str] = None,
    *,
    log_to_console: bool = False,
) -> None:
    """Execute ENVI FLAASH with the provided parameter dictionary."""
    _log("Processing", enabled=log_to_console, step="flaash")
    if output_image_path_to_delete:
        remove_output_files([output_image_path_to_delete])
    task = envi_engine.task("FLAASH")
    task.execute(flaash_params)
    with open(output_params_path, "w", encoding="utf-8") as file:
        file.write(str(flaash_params))
    _log("Wrote output", enabled=log_to_console, step="flaash")


@reports_progress
def _run_flaash_wrapper(args: tuple[Dict[str, Any], str, Any]) -> str:
    """Executor wrapper for FLAASH grid runs."""
    test_params, test_output_params_path, envi_engine = args
    _execute_flaash_task(test_params, test_output_params_path, envi_engine)
    return test_params["OUTPUT_RASTER_URI"]


@reports_progress(worker_progress=True)
def parallel_flaash(
    test_flaash_params_array: List[tuple[Dict[str, Any], str]],
    envi_engine: Any,
    max_workers: int = 4,
) -> List[str]:
    """Execute `run_flaash` in parallel on parameterized FLAASH runs."""
    all_output_paths = []
    tasks = [
        (test_params, test_output_params_path, envi_engine)
        for test_params, test_output_params_path in test_flaash_params_array
    ]

    with (
        local_worker_progress(current_callback(), processes=True) as reporter,
        ProcessPoolExecutor(max_workers=max_workers) as executor,
    ):
        callback_kwargs = {"progress_callback": reporter} if reporter is not None else {}
        futures = [executor.submit(_run_flaash_wrapper, task, **callback_kwargs) for task in tasks]
        for future in progress(as_completed(futures), total=len(futures), desc="FLAASH", disable=False):
            output_uri = future.result()
            all_output_paths.append(output_uri)

    return all_output_paths


class FLAASHCorrector:
    """Adapter around existing ENVI FLAASH execution helper."""

    def __init__(self, envi_engine: Any) -> None:
        """Initialize the FLAASH adapter.
        Args:
            envi_engine: Initialized ENVI engine instance.
        Returns:
            None.
        """
        self.envi_engine = envi_engine

    def run(self, input_raster: str, output_raster: str, **kwargs: Any) -> str:
        """Run FLAASH correction and return the output raster path."""
        params: Dict[str, Any] = dict(kwargs)
        params.setdefault("INPUT_RASTER", {"url": input_raster, "factory": "URLRaster"})
        params.setdefault("OUTPUT_RASTER_URI", output_raster)
        params_path = kwargs.get("params_path")
        if not params_path:
            raise ValueError("params_path is required for FLAASHCorrector.run")
        _execute_flaash_task(params, params_path, self.envi_engine, output_raster)
        return output_raster


@_logged_operation("flaash_atmospheric_correction", inputs=("input_raster",), outputs=("output_raster",))
def run_flaash(
    input_raster: str,
    output_raster: str,
    *,
    variables: Any = None,
    dem_file_path: Optional[str] = None,
    footprint_geometry: dict | BaseGeometry | None = None,
    envi_engine_path: Optional[str] = None,
    envi_engine: Any = None,
    output_params_path: Optional[str] = None,
    convert_paths_for_windows: bool = False,
    delete_output_before_run: Optional[str] = None,
    params: Optional[Dict[str, Any]] = None,
    dem_ground_percentile: float = 50.0,
    modtran_atm: str = "Mid-Latitude Summer",
    modtran_aer: str = "Rural",
    use_aerosol: str = "AUTO",
    default_visibility: Optional[float] = None,
    custom_params: Optional[Dict[str, Any]] = None,
    log_to_console: bool = False,
) -> FLAASHRunResult:
    """Run FLAASH using mapped variables JSON or explicit params.
    Args:
        input_raster: Input raster path.
        output_raster: Output raster path.
        variables: Optional mapped variables JSON object.
        dem_file_path: Optional DEM raster path.
        footprint_geometry: Optional Shapely footprint in EPSG:4326.
        envi_engine_path: Optional ENVI engine executable path.
        envi_engine: Optional initialized ENVI engine.
        output_params_path: Optional executed-params output path.
        convert_paths_for_windows: Whether to convert paths for Windows ENVI execution.
        delete_output_before_run: Optional output path to delete before execution.
        params: Optional explicit FLAASH parameter dictionary.
        dem_ground_percentile: DEM percentile used for ground elevation estimation.
        modtran_atm: MODTRAN atmosphere profile name.
        modtran_aer: MODTRAN aerosol profile name.
        use_aerosol: FLAASH aerosol handling mode.
        default_visibility: Optional default visibility override.
        custom_params: Optional custom FLAASH parameter overrides.
        log_to_console: Whether to emit console logs.
    Returns:
        FLAASH run summary.
    """
    if not output_raster:
        raise ValueError("output_raster must be an explicit output path")
    output_raster = str(output_raster)
    if footprint_geometry is None and variables is not None:
        footprint_geometry = variables.get("geometry")
    if isinstance(footprint_geometry, dict):
        from vhrharmonize.io.metadata import materialize_geometry

        footprint_geometry = materialize_geometry(footprint_geometry)
    if params is None:
        if variables is None or dem_file_path is None or footprint_geometry is None:
            raise ValueError(
                "variables, dem_file_path, and footprint_geometry are required when params is not provided."
            )
        params = _build_flaash_kwargs_from_variables(
            input_raster=input_raster,
            dem_file_path=dem_file_path,
            footprint_geometry=footprint_geometry,
            variables=variables,
            output_raster=output_raster,
            dem_ground_percentile=dem_ground_percentile,
            modtran_atm=modtran_atm,
            modtran_aer=modtran_aer,
            use_aerosol=use_aerosol,
            default_visibility=default_visibility,
            custom_params=custom_params,
        )
    else:
        params = _validate_flaash_params(dict(params))
        params.setdefault("INPUT_RASTER", {"url": input_raster, "factory": "URLRaster"})
        params.setdefault("OUTPUT_RASTER_URI", output_raster)

    params_output_path = output_params_path or f"{output_raster}.flaash_params.txt"
    params_to_run = (
        _convert_flaash_params_paths_for_windows(params) if convert_paths_for_windows else params
    )

    resolved_engine = envi_engine
    if resolved_engine is None:
        if not envi_engine_path:
            raise ValueError("envi_engine or envi_engine_path is required for FLAASH.")
        resolved_engine = _init_envi_engine(envi_engine_path)

    _execute_flaash_task(
        params_to_run,
        params_output_path,
        resolved_engine,
        delete_output_before_run or output_raster,
        log_to_console=log_to_console,
    )
    return FLAASHRunResult(
        output_raster=output_raster,
        params=params,
        params_output_path=params_output_path,
    )


class FLAASHAtmosphericCorrection(FunctionPlugin):
    target = "vhrharmonize.plugins.flaash_atmospheric_correction:run_flaash"
    input_dependency_paths = frozenset({"dem_file_path", "input_raster"})
    input_existence_check_paths = input_dependency_paths
    input_protection_paths = input_dependency_paths
    input_hpc_staging_paths = input_dependency_paths
    output_path_resolution_paths = frozenset({"output_params_path", "output_raster"})
    output_dependency_paths = output_path_resolution_paths
    output_parent_creation_paths = output_path_resolution_paths
    output_collision_check_paths = output_path_resolution_paths
    output_reuse_paths = output_path_resolution_paths
    output_validation_paths = output_path_resolution_paths
    output_invalid_removal_paths = output_path_resolution_paths
    output_overview_calculation_paths = frozenset({"output_raster"})
    output_temporary_cleanup_paths = output_path_resolution_paths
    output_hpc_staging_paths = output_path_resolution_paths
    output_hpc_download_paths = output_path_resolution_paths


__all__ = [
    "FLAASHAtmosphericCorrection", "FLAASHRunResult", "FLAASHCorrector",
    "FLAASH_ALLOWED_PARAMS", "run_flaash", "parallel_flaash",
]
