from copy import deepcopy

import pytest
import yaml

from vhrharmonize.workflow.yaml_document import rewrite_yaml, preparation_config


def test_changed_block_scalar_keeps_header_comment_and_following_lines():
    source = "root: >- # root location\n  /local/project\n# next setting\nlabel: original\n"
    original = yaml.safe_load(source)
    result = rewrite_yaml(source, original, {**original, "root": "/remote/project"})
    assert "# root location\n# next setting\nlabel: original\n" in result
    assert yaml.safe_load(result)["root"] == "/remote/project"


def test_crlf_and_unchanged_quoted_values_survive_added_settings():
    source = '# recipe\r\nstep:\r\n  plugin: "file_source" # keep\r\n'
    original = yaml.safe_load(source)
    updated = {"step": {**original["step"], "core:run": False}}
    result = rewrite_yaml(source, original, updated)
    assert result.startswith(source)
    assert '\n' not in result.replace('\r\n', '')


def test_selective_edits_preserve_comments_quotes_and_unmodified_text():
    source = '''# Heading
settings:
  plugin: shared
  core:run: true
  const:root: '/local/project' # relocate this

# Keep the expression and comments exactly.
process:
  plugin: file_source
  core:run: true # enabled
  var:image: "expr:const.root & '/scene.tif'"
  param:output_path: var:image
'''
    original = yaml.safe_load(source)
    updated = deepcopy(original)
    updated["settings"]["const:root"] = "path:/remote/project"
    updated["process"]["core:run"] = False
    actual = rewrite_yaml(source, original, updated)
    assert actual == source.replace("'/local/project'", "'path:/remote/project'").replace("true # enabled", "false # enabled")


def test_context_filename_key_and_new_controls_preserve_comments():
    source = '''save:
  core:run: true
  core:save_context:
    "path:./local.json": # snapshot
      - defined # fields
# Last comment
'''
    original = yaml.safe_load(source)
    updated = deepcopy(original)
    updated["save"]["core:save_context"] = {"path:/remote/context.json": ["defined"]}
    updated["save"]["core:require_outputs"] = True
    actual = rewrite_yaml(source, original, updated)
    assert '# snapshot' in actual and '# fields' in actual and '# Last comment' in actual
    assert yaml.safe_load(actual) == updated


def test_alias_override_does_not_modify_the_anchor():
    source = 'first: &root "/local" # anchor\nsecond: *root # alias\n'
    original = yaml.safe_load(source)
    actual = rewrite_yaml(source, original, {"first": "/local", "second": "/remote"})
    assert actual.startswith('first: &root "/local" # anchor\n')
    assert '# alias' in actual


def test_preparation_disables_only_later_steps_and_keeps_parameters():
    original = {
        "off": {"plugin": "import_files", "core:run": False},
        "import": {"plugin": "import_files", "core:run": True, "param:search_glob": "*.tif"},
        "later": {"plugin": "file_source", "core:run": True, "param:input_path": "var:raw"},
    }
    expected = deepcopy(original)
    expected["later"]["core:run"] = False
    assert preparation_config(original, "import") == expected
    assert original["later"]["core:run"] is True
    with pytest.raises(ValueError, match="Unknown preparation step"):
        preparation_config(original, "missing")
