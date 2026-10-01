from pathlib import Path

import pytest

from test_explicit_context import shared, importer
from workflow_helpers import install_function, transfer
from vhrharmonize.workflow.engine import Workflow
from vhrharmonize.workflow.staging import stage_workflow
from vhrharmonize.workflow.yaml_document import hpc_preparation_config


def prepare(config, tmp_path, cutoff):
    workflow = Workflow(hpc_preparation_config(config, cutoff), config_dir=tmp_path, preparing=True)
    workflow.run()
    return workflow


def stage(config, workflow, tmp_path):
    return stage_workflow(config, config_dir=tmp_path, workflow=workflow,
                          remote_work_dir=str(tmp_path / 'remote'),
                          path_mappings={'const:root': str(tmp_path / 'remote')})


def test_import_called_once_without_context_cache(monkeypatch, tmp_path):
    import importlib
    module = importlib.import_module('vhrharmonize.plugins.import_files')
    original = module.import_files
    calls = []
    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)
    # Keep the function signature used by argument validation.
    import functools
    monkeypatch.setattr(module, 'import_files', functools.wraps(original)(counted))
    source = tmp_path / 'raw.txt'
    source.write_text('raw')
    config = {'settings': shared(**{'const:root': str(tmp_path)}), 'import': importer(source),
              'copy': {'plugin': 'file_source', 'core:run': True, 'core:require_outputs': True,
                       'param:input_path': 'var:raw', 'param:output_path': "expr:const.root & '/result.txt'"}}
    workflow = prepare(config, tmp_path, 'import')
    assert calls == [1]
    assert not (tmp_path / 'result.txt').exists()
    staged, uploads, _ = stage(config, workflow, tmp_path)
    assert calls == [1]
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert calls == [1]
    assert (tmp_path / 'remote/result.txt').read_text() == 'raw'


def test_constants_only_null_cutoff_never_calls_locally(monkeypatch, tmp_path):
    calls = []
    def create(output_path, value):
        calls.append(value)
        Path(output_path).write_text(str(value))
    install_function(monkeypatch, 'create', create, output_paths=('output_path',))
    config = {'settings': shared(**{'const:root': str(tmp_path), 'const:value': 7}),
              'create': {'plugin': 'create', 'core:run': True, 'core:require_outputs': True,
                         'param:output_path': "expr:const.root & '/result.txt'", 'param:value': 'const:value'}}
    workflow = prepare(config, tmp_path, None)
    assert calls == []
    assert 'var' not in workflow.initial_context
    assert workflow.nodes[0].status == 'processing'
    staged, uploads, _ = stage(config, workflow, tmp_path)
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert calls == [7]


def test_local_return_resolves_remote_paths_without_repeating_producer(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, 'choose', lambda: calls.append('choose') or {'name': 'chosen'})
    def create(output_path):
        calls.append('create')
        Path(output_path).write_text('result')
    install_function(monkeypatch, 'create', create, output_paths=('output_path',))
    config = {'settings': shared(**{'const:root': str(tmp_path)}),
              'choose': {'plugin': 'choose', 'core:run': True, 'const:name': 'returned:name'},
              'create': {'plugin': 'create', 'core:run': True, 'core:require_outputs': True,
                         'param:output_path': "expr:const.root & '/' & const.name & '.txt'"}}
    workflow = prepare(config, tmp_path, 'choose')
    assert calls == ['choose']
    assert workflow.nodes[0].params['output_path'] == str(tmp_path / 'chosen.txt')
    staged, uploads, _ = stage(config, workflow, tmp_path)
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert calls == ['choose', 'create']
    assert (tmp_path / 'remote/chosen.txt').exists()


def test_unknown_remote_output_stays_needed(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, 'choose', lambda: calls.append('choose') or {'name': 'chosen'})
    install_function(monkeypatch, 'create', lambda output_path: calls.append('create'), output_paths=('output_path',))
    config = {'settings': shared(**{'const:root': str(tmp_path)}),
              'choose': {'plugin': 'choose', 'core:run': True, 'const:name': 'returned:name'},
              'create': {'plugin': 'create', 'core:run': True, 'core:require_outputs': True,
                         'param:output_path': "expr:const.root & '/' & const.name & '.txt'"}}
    workflow = prepare(config, tmp_path, None)
    assert calls == []
    assert all(node.needed and node.status == 'processing' for node in workflow.nodes)
    stage(config, workflow, tmp_path)


def test_var_requires_initialized_records(tmp_path):
    with pytest.raises(ValueError, match='var is unavailable'):
        Workflow({'settings': shared(), 'bad': {'core:run': True, 'var:value': 1}}, config_dir=tmp_path)


def test_canonical_var_plugin_contract_and_scope(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, 'records', lambda: {'items': [{'value': 1}, {'value': 2}]},
                     var_records_return='items')
    install_function(monkeypatch, 'consume', lambda value: calls.append(value))
    Workflow({'settings': shared(), 'records': {'plugin': 'records', 'core:run': True},
              'consume': {'plugin': 'consume', 'core:run': True, 'core:scope': 'var',
                          'core:require_outputs': True, 'param:value': 'var:value'}}, config_dir=tmp_path).run()
    assert calls == [1, 2]


def test_local_file_result_is_uploaded_and_prefix_is_not_repeated(tmp_path):
    source = tmp_path / 'raw.txt'
    source.write_text('raw')
    config = {'settings': shared(**{'const:root': str(tmp_path)}), 'import': importer(source),
              'local': {'plugin': 'file_source', 'core:run': True,
                        'param:input_path': 'var:raw', 'var:prepared': "expr:const.root & '/local.txt'",
                        'param:output_path': 'var:prepared'},
              'remote': {'plugin': 'file_source', 'core:run': True, 'core:require_outputs': True,
                         'param:input_path': 'var:prepared', 'param:output_path': "expr:const.root & '/result.txt'"}}
    workflow = prepare(config, tmp_path, 'local')
    assert (tmp_path / 'local.txt').read_text() == 'raw'
    assert not (tmp_path / 'result.txt').exists()
    staged, uploads, _ = stage(config, workflow, tmp_path)
    assert str(tmp_path / 'local.txt') in uploads
    assert str(source) not in uploads
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert (tmp_path / 'remote/result.txt').read_text() == 'raw'


def test_unknown_remote_outputs_execute_and_download_by_directory(monkeypatch, tmp_path):
    install_function(monkeypatch, 'choose', lambda: {'name': 'chosen'})
    def create(output_path):
        Path(output_path).write_text('dynamic')
    install_function(monkeypatch, 'create', create, output_paths=('output_path',))
    config = {'settings': shared(**{'const:root': str(tmp_path)}),
              'choose': {'plugin': 'choose', 'core:run': True, 'const:name': 'returned:name'},
              'create': {'plugin': 'create', 'core:run': True, 'core:require_outputs': True,
                         'param:output_path': "expr:const.root & '/' & const.name & '.txt'"}}
    workflow = prepare(config, tmp_path, None)
    staged, uploads, downloads = stage(config, workflow, tmp_path)
    assert downloads[str(tmp_path)] == str(tmp_path / 'remote')
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert (tmp_path / 'remote/chosen.txt').read_text() == 'dynamic'


def test_unresolved_external_input_requires_explicit_preparation(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, 'choose', lambda: calls.append(1) or str(tmp_path / 'external.txt'))
    config = {'settings': shared(**{'const:root': str(tmp_path)}),
              'choose': {'plugin': 'choose', 'core:run': True, 'const:external': 'returned:$'},
              'copy': {'plugin': 'file_source', 'core:run': True, 'core:require_outputs': True,
                       'param:input_path': 'const:external', 'param:output_path': "expr:const.root & '/result.txt'"}}
    workflow = prepare(config, tmp_path, None)
    with pytest.raises(ValueError, match='unresolved HPC input'):
        stage(config, workflow, tmp_path)
    assert calls == []


def test_null_cutoff_does_not_call_record_producer(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, 'records', lambda: calls.append(1) or [{'value': 1}], var_records_return='$')
    config = {'settings': shared(**{'const:root': str(tmp_path)}),
              'records': {'plugin': 'records', 'core:run': True}}
    workflow = prepare(config, tmp_path, None)
    assert calls == []
    with pytest.raises(ValueError, match='run_to_step_before_prepare'):
        stage(config, workflow, tmp_path)
    assert calls == []


def test_explicit_var_scope_without_records_is_rejected(monkeypatch, tmp_path):
    install_function(monkeypatch, 'consume', lambda: None)
    with pytest.raises(ValueError, match='requires initialized var records'):
        Workflow({'settings': shared(), 'consume': {'plugin': 'consume', 'core:run': True,
                                                    'core:scope': 'var'}}, config_dir=tmp_path)


def test_local_prefix_collects_return_values_even_without_file_targets(monkeypatch, tmp_path):
    calls = []
    install_function(monkeypatch, 'first', lambda: calls.append('first') or 3)
    install_function(monkeypatch, 'second', lambda: calls.append('second') or 4)
    config = {'settings': shared(),
              'first': {'plugin': 'first', 'core:run': True, 'const:first': 'returned:$'},
              'second': {'plugin': 'second', 'core:run': True, 'const:second': 'returned:$'}}
    workflow = prepare(config, tmp_path, 'second')
    assert calls == ['first', 'second']
    assert workflow.preparation_state['const'] == {'first': 3, 'second': 4}


def test_null_cutoff_can_restore_explicit_var_context_without_import(monkeypatch, tmp_path):
    source = tmp_path / 'raw.txt'
    source.write_text('raw')
    snapshot = tmp_path / 'saved.json'
    config = {'settings': shared(**{'const:root': str(tmp_path)}),
              'import': {**importer(source), 'core:save_context': {str(snapshot): 'all'},
                         'core:load_context': {str(snapshot): 'all'}},
              'copy': {'plugin': 'file_source', 'core:run': True, 'core:require_outputs': True,
                       'param:input_path': 'var:raw', 'param:output_path': "expr:const.root & '/result.txt'"}}
    Workflow(config, config_dir=tmp_path)  # Explicitly saved state from a prior run.
    import importlib
    module = importlib.import_module('vhrharmonize.plugins.import_files')
    import functools
    original = module.import_files
    @functools.wraps(original)
    def unexpected(*args, **kwargs):
        pytest.fail('import_files called despite restored state and a null cutoff')
    monkeypatch.setattr(module, 'import_files', unexpected)
    workflow = prepare(config, tmp_path, None)
    staged, uploads, _ = stage(config, workflow, tmp_path)
    transfer(uploads)
    Workflow(staged, config_dir=tmp_path).run()
    assert (tmp_path / 'remote/result.txt').read_text() == 'raw'


def test_prepared_import_preserves_satisfied_output_bindings(monkeypatch, tmp_path):
    from test_discovery_plugin_contract import discovery_recipe
    config, _, mappings, calls, _ = discovery_recipe(tmp_path, monkeypatch)
    config['discover_inputs']['core:satisfies'] = {'mask': 'output_path'}
    finish = config.pop('copy')
    config['mask'] = {'plugin': 'file_source', 'core:run': True,
                      'param:input_path': 'var:missing_upstream',
                      'var:masked': "expr:const.root & '/' & var.identity.key & '.txt'",
                      'param:output_path': 'var:masked'}
    finish['param:input_path'] = 'var:masked'
    config['copy'] = finish
    workflow = prepare(config, tmp_path, 'discover_inputs')
    staged, uploads, _ = stage_workflow(config, config_dir=tmp_path, workflow=workflow,
                                       remote_work_dir=str(tmp_path / 'remote'), path_mappings=mappings)
    transfer(uploads)
    remote = Workflow(staged, config_dir=tmp_path)
    remote.run()
    assert len(calls) == 1
    assert remote.counts()['mask']['processing'] == 0
    assert [(tmp_path / f'remote/products/{key}/done.txt').read_text() for key in ('a', 'b')] == ['a', 'b']


@pytest.mark.parametrize('declaration', ['scene_records_return', 'scene_records_mode',
                                        'scene_id_return', 'scene_path_return'])
@pytest.mark.parametrize('on_class', [False, True])
def test_removed_scene_declarations_are_rejected(declaration, on_class):
    from vhrharmonize.plugins.base import FunctionPlugin
    if on_class:
        plugin = type('OldPlugin', (FunctionPlugin,), {declaration: 'items'})()
    else:
        plugin = FunctionPlugin()
        setattr(plugin, declaration, 'items')
    with pytest.raises(ValueError, match='Removed plugin declarations: ' + declaration):
        plugin.file_features()
    assert plugin.var_records_return is None


def test_removed_scene_scope_is_rejected():
    from vhrharmonize.plugins.base import FunctionPlugin
    from vhrharmonize.workflow.config import validate_config
    with pytest.raises(ValueError, match='core:scope must be var or aggregate'):
        validate_config({'step': {'core:run': True, 'core:scope': 'scene'}})
    plugin = FunctionPlugin()
    plugin.scope = 'scene'
    with pytest.raises(ValueError, match='Plugin scope must be var or aggregate'):
        plugin.file_features()


def test_discovery_cutoff_keeps_existing_graph_and_validation(monkeypatch, tmp_path):
    source = tmp_path / 'raw.txt'
    source.write_text('raw')
    config = {'settings': shared(**{'const:root': str(tmp_path), 'core:report_progress': True}),
              'import': importer(source),
              'copy': {'plugin': 'file_source', 'core:run': True, 'core:require_outputs': True,
                       'param:input_path': 'var:raw', 'param:output_path': str(tmp_path / 'result.txt')}}
    workflow = Workflow(hpc_preparation_config(config, 'import'), config_dir=tmp_path, preparing=True)
    workflow.plan()
    nodes = list(workflow.nodes)
    monkeypatch.setattr(workflow, '_build', lambda *args: pytest.fail('Duplicate build'))
    monkeypatch.setattr(workflow, '_valid', lambda *args: pytest.fail('Repeated file validation'))
    workflow.run()
    assert all(before is after for before, after in zip(nodes, workflow.nodes))
    row = next(row for row in workflow.get_progress()['rows'] if row['name'] == 'copy')
    assert row['disabled'] is True
    assert row['active'] == 0


def test_remote_suffix_is_not_built_before_local_return(monkeypatch, tmp_path):
    install_function(monkeypatch, 'choose', lambda: {'name': 'chosen'})
    install_function(monkeypatch, 'create', lambda output_path: None, output_paths=('output_path',))
    config = {'settings': shared(**{'const:root': str(tmp_path)}),
              'choose': {'plugin': 'choose', 'core:run': True, 'const:name': 'returned:name'},
              'create': {'plugin': 'create', 'core:run': True, 'core:require_outputs': True,
                         'param:output_path': "expr:const.root & '/' & const.name & '.txt'"}}
    workflow = Workflow(hpc_preparation_config(config, 'choose'), config_dir=tmp_path, preparing=True)
    assert [node.step['name'] for node in workflow.nodes] == ['choose']
    workflow.run()
    assert [node.step['name'] for node in workflow.nodes] == ['create']
    assert workflow.nodes[0].params['output_path'] == str(tmp_path / 'chosen.txt')
