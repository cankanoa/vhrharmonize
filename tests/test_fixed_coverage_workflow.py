"""Joint registration exposes native reference-coverage policies through YAML."""

from unittest.mock import Mock

import pytest

from vhrharmonize.parameters import function_parameters
from vhrharmonize.plugins import joint_coregistration as adapter
from vhrharmonize.workflow.registry import load_plugin


@pytest.mark.parametrize("mode", ["error", "no_influence", "partial_influence"])
def test_fixed_reference_coverage_policy_is_exposed_and_forwarded(monkeypatch, mode):
    parameter = function_parameters(adapter.joint_coregistration)["fixed_on_no_or_partial_overlap"]
    assert parameter.default == "partial_influence"
    native = Mock(return_value=["moving_out.tif"])
    monkeypatch.setattr(adapter, "_upstream", lambda: native)
    load_plugin("joint_coregistration").run(
        params={"input_images": ["moving.tif"],
                "output_images": ["moving_out.tif"],
                "fixed_images": ["reference.tif"], "fixed_on_no_or_partial_overlap": mode},
        shared={},
    )
    assert native.call_args.kwargs["fixed_on_no_or_partial_overlap"] == mode


def test_fixed_reference_paths_use_core_input_features():
    plugin = load_plugin("joint_coregistration")
    for name in ("input_dependency_paths", "input_existence_check_paths",
                 "input_protection_paths", "input_hpc_staging_paths"):
        assert "fixed_images" in getattr(plugin, name)


def test_matching_band_lists_are_forwarded_without_tuple_conversion(monkeypatch):
    native = Mock(return_value=["out.tif"])
    monkeypatch.setattr(adapter, "_upstream", lambda: native)
    load_plugin("joint_coregistration").run(
        params={"input_images":["moving.tif"], "fixed_images":["reference.tif"],
                "output_images":["out.tif"], "input_images_bands":[2,4], "fixed_images_bands":[3]},
        shared={})
    assert native.call_args.kwargs["input_images_bands"] == [2,4]
    assert native.call_args.kwargs["fixed_images_bands"] == [3]


def test_fixed_reference_paths_are_uploaded_and_remapped_for_hpc(tmp_path, make_test_raster):
    from vhrharmonize.workflow.staging import stage_workflow
    from vhrharmonize.workflow.engine import Workflow
    moving = make_test_raster(tmp_path/'moving.tif', crs='EPSG:3857')
    reference = make_test_raster(tmp_path/'reference.tif', crs='EPSG:3857')
    config = {
        'shared': {'plugin':'shared', 'core:run':True, 'core:save_statistics_path':None,
                   'core:load_statistics_path':None, 'const:root':str(tmp_path)},
        'joint': {'plugin':'joint_coregistration', 'core:run':True,
                  'core:require_outputs':'param:output_images', 'param:input_images':[str(moving)],
                  'param:fixed_images':[str(reference)], 'param:output_images':[str(tmp_path/'aligned.tif')]},
    }
    staged, uploads, _ = stage_workflow(
        config, config_dir=tmp_path, remote_work_dir=str(tmp_path/'remote'),
        path_mappings={'const:root':str(tmp_path/'remote')})
    assert str(reference) in uploads and str(moving) in uploads
    node = Workflow(staged, config_dir=tmp_path).plan().nodes[0]
    assert node.params['fixed_images'] == [uploads[str(reference)]]
    assert node.params['input_images'] == [uploads[str(moving)]]


@pytest.mark.parametrize('moving_method', ['orb', 'elastix'])
@pytest.mark.parametrize('fixed_method', ['orb', 'elastix'])
def test_independent_matcher_options_are_forwarded(monkeypatch, moving_method, fixed_method):
    native = Mock(return_value=['out.tif'])
    monkeypatch.setattr(adapter, '_upstream', lambda: native)
    load_plugin('joint_coregistration').run(
        params={'input_images':['a.tif','b.tif'], 'fixed_images':['reference.tif'],
                'output_images':['a_out.tif','b_out.tif'],
                'moving_tie_method':moving_method, 'fixed_tie_method':fixed_method}, shared={})
    assert native.call_args.kwargs['moving_tie_method'] == moving_method
    assert native.call_args.kwargs['fixed_tie_method'] == fixed_method
