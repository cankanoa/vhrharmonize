"""Native overlap placement and pair quotas survive VHR argument resolution."""
import importlib
from pathlib import Path
from unittest.mock import Mock

import pytest
import yaml

from vhrharmonize.parameters import function_parameters
from vhrharmonize.workflow.registry import load_plugin


@pytest.mark.parametrize('name,default', [('joint_coregistration', 'overlap_spread'), ('global_regression', 'overlap_spread')])
def test_native_sampling_defaults_and_quota_forwarding(name, default, monkeypatch):
    module = importlib.import_module('vhrharmonize.plugins.' + name)
    parameters = function_parameters(getattr(module, name))
    assert parameters['tie_sampling_method'].default == default
    assert parameters['tie_max_samples'].default == '500crs'
    assert parameters['tie_min_samples'].default == '1000crs'
    native = Mock(return_value=['out.tif'])
    monkeypatch.setattr(module, '_upstream', lambda: native)
    params = dict(input_images=['a.tif', 'b.tif'], output_images=['out.tif'],
                  tie_sampling_method='overlap_spread', tie_max_samples='500crs',
                  tie_min_samples='1cnt', tie_max_attempts=100)
    assert load_plugin(name).run(params=params, shared={}) == ['out.tif']
    native.assert_called_once_with(**params)


def test_worldview_pair_quotas_are_separate_from_pif_pixels():
    root = Path(__file__).parents[1]
    paths = [root/'configs/example.worldview.yml', *sorted((root/'configs/drylands_groups').glob('*/*.worldview.yml'))]
    for path in paths:
        config = yaml.safe_load(path.read_text())
        for step in ('joint_coregistration', 'global_regression'):
            for name in ('tie_sampling_method', 'tie_max_samples', 'tie_min_samples', 'tie_max_attempts'):
                assert 'param:' + name not in config[step], (path, step, name)
            assert '# param:tie_sampling_method: overlap_spread' in path.read_text()
        assert config['global_regression']['param:pif_region_radius'] == 30, path


def test_accepted_tie_progress_shows_count_and_updates_eta_while_waiting(monkeypatch):
    from prompt_toolkit.formatted_text import to_plain_text
    from vhrharmonize.workflow import substeps
    from vhrharmonize.workflow.progress_terminal import TerminalProgressDisplay
    from spectralmatch.substeps import Substep

    now = [100.0]
    monkeypatch.setattr(substeps, 'monotonic', lambda: now[0])
    state = substeps.Substeps()
    phase = Substep(lambda **stats: state.handle(stats), 'Matching overlaps', total=1, unit='pairs', kind='jobs')
    phase.emit('start')
    reporter = phase.reporter(lambda **stats: state.handle(stats), 0, 'image A ↔ image B')
    reporter.lifecycle('job_start')
    reporter(n=2, total=8, prefix='Accepted tie pairs', unit='pairs', elapsed=10, rate=.2, operation='ties')
    now[0] += 10
    active = state.operations('task', 'joint_coregistration')
    assert active[0]['stats']['n'] == 2 and active[0]['stats']['total'] == 8
    assert active[0]['eta_seconds'] == 60  # Includes the wait since the last successful match.
    display = TerminalProgressDisplay(width=140)
    text = '\n'.join(to_plain_text(line) for line in display.operations(active, width=140))
    assert '2/8 pairs' in text
    reporter.lifecycle('job_finished')
    assert state.operations('task', 'joint_coregistration') == []
