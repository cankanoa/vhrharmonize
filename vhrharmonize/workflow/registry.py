"""Processing plugins registered through Python package entry points."""

from importlib import import_module, metadata


GROUP = "vhrharmonize.plugins"
BUILTINS = {
    "joint_coregistration": "joint_coregistration:JointCoregistration",
    "global_regression": "global_regression:GlobalRegression",
    "local_block_adjustment": "local_block_adjustment:LocalBlockAdjustment",
    "align_rasters": "align_rasters:AlignRasters",
    "create_footprints": "create_footprints:CreateFootprints",
    "postprocess_footprints": "postprocess_footprints:PostprocessFootprints",
    "voronoi_center_seamline": "voronoi_center_seamline:VoronoiCenterSeamline",
    "weighted_seamline": "weighted_seamline:WeightedSeamline",
    "markov_triangles": "markov_triangles:MarkovTriangles",
    "mask_rasters": "mask_rasters:MaskRasters",
    "merge_rasters": "merge_rasters:MergeRasters",
    "merge_vectors": "merge_vectors:MergeVectors",
    "band_math": "band_math:BandMath",
    "create_cloud_mask_with_omnicloudmask": "create_cloud_mask_with_omnicloudmask:CreateCloudMaskWithOmnicloudmask",
    "process_raster_values_to_vector_polygons": "process_raster_values_to_vector_polygons:ProcessRasterValuesToVectorPolygons",
    "compute_overviews": "compute_overviews:ComputeOverviews",
    "search_paths": "search_paths:SearchPaths",
    "create_paths": "create_paths:CreatePaths",
    "match_paths": "match_paths:MatchPaths",
    "import_files": "import_files:ImportFiles",
    "file_source": "file_source:FileSource",
    "fetch_atmosphere": "fetch_atmosphere:FetchAtmosphere",
    "fetch_dem": "fetch_dem:FetchDEM",
    "atmospheric_correction": "atmospheric_correction:AtmosphericCorrection",
    "orthorectification": "orthorectification:Orthorectification",
    "pansharpen": "pansharpen:Pansharpen",
    "cloud_mask": "cloud_mask:CloudMask",
    "alignment": "alignment:Alignment",
    "seamline_metadata": "seamline_metadata:SeamlineMetadata",
}


def _entries():
    entries = metadata.entry_points()
    selected = entries.select(group=GROUP) if hasattr(entries, "select") else entries.get(GROUP, [])
    return [
        entry
        for entry in selected
        if entry.name in BUILTINS
        or not getattr(entry, "value", "").startswith("vhrharmonize.plugins.")
    ]


def plugin_names():
    """List registered plugins without importing their processing dependencies."""
    return sorted(set(BUILTINS) | {entry.name for entry in _entries()})


def load_plugin(name):
    if name is None:
        from .settings import ContextStep

        return ContextStep()
    matches = [entry for entry in _entries() if entry.name == name]
    if len(matches) > 1:
        raise ValueError(f"Multiple plugins registered as {name!r}")
    if matches:
        return matches[0].load()()
    # Source checkouts work before reinstalling their entry-point metadata.
    if name in BUILTINS:
        module, cls = BUILTINS[name].split(":")
        return getattr(import_module(f"vhrharmonize.plugins.{module}"), cls)()
    raise ValueError(f"Unknown workflow plugin {name!r}; register it in {GROUP}")
