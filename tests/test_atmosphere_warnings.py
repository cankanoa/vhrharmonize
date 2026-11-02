"""Atmosphere lookup warnings, unchanged retrieval order, and persisted provenance."""

from dataclasses import asdict
from datetime import date, datetime, timedelta, timezone
import importlib
import json
from pathlib import Path

import pytest
import requests
import yaml

from vhrharmonize.workflow.values import resolve

fetch = importlib.import_module('vhrharmonize.plugins.fetch_atmosphere')
DEFAULTS = dict(aot550=.2, water_vapor=2.5, ozone_cm_atm=.3)
BBOX = dict(min_lon=0, min_lat=0, max_lon=1, max_lat=1)
DAY = date(2020, 1, 2)
TIME = datetime(2020, 1, 2, 12, tzinfo=timezone.utc)


def codes(result, code):
    return [item for item in result.warnings if item['code'] == code]


def nasa(monkeypatch, callback, **kwargs):
    monkeypatch.setattr(fetch, '_fetch_power_daily_point', callback)
    return fetch.fetch_power_atmosphere_for_bbox(day_utc=DAY, grid_size=1, **BBOX, **kwargs)


def candidate(collection, *, parameter='water_vapor', raw=None, timestamp=None, band_index=0):
    bands = (['Atmospheric_Water_Vapor_Mean'] if parameter == 'water_vapor' else
             ['Aerosol_Optical_Depth_Land_Ocean_Mean',
              'Aerosol_Optical_Depth_Land_Ocean_Mean_Mean',
              'Aerosol_Optical_Depth_Land_Ocean'])
    timestamp = timestamp or TIME
    return dict(collection=collection, band_found=True,
        bands=bands[band_index:], selected_band=bands[band_index],
        raw_value=(2000 if parameter == 'water_vapor' else 100) if raw is None else raw,
        image_time_utc=timestamp.isoformat(),
        abs_time_diff_hours=abs((timestamp - TIME).total_seconds()) / 3600)


def modis_mocks(monkeypatch, water=None, aerosol=None):
    monkeypatch.setattr(fetch, '_fetch_collection_value',
        water or (lambda ee, collection, *a, **kw: candidate(collection)))
    monkeypatch.setattr(fetch, '_fetch_first_available_band_value',
        aerosol or (lambda ee, collection, *a, **kw: candidate(collection, parameter='aot550')))


def modis(**kwargs):
    return fetch.fetch_modis_water_vapor_for_bbox(scene_datetime_utc=TIME, ee=object(), **BBOX, **kwargs)


def test_nasa_complete_data_has_no_warnings_or_substitutions(monkeypatch, capsys):
    result = nasa(monkeypatch, lambda *a, **kw: dict(DEFAULTS), fallback_values=DEFAULTS,
                  log_to_console=True)
    assert not result.warnings and result.defaults_used == {}
    assert result.lookup_attempts[0]['status'] == 'ok'
    assert result.date_used == result.date_requested == DAY.isoformat()
    assert 'WARNING' not in capsys.readouterr().out


def test_nasa_missing_fields_warn_for_every_point_and_only_default_missing_fields(monkeypatch, capsys):
    calls = []
    def lookup(day, lon, lat, **kwargs):
        calls.append((day, lon, lat))
        return dict(aot550=None if len(calls) == 1 else .1, water_vapor=1., ozone_cm_atm=None)
    monkeypatch.setattr(fetch, '_fetch_power_daily_point', lookup)
    result = fetch.fetch_power_atmosphere_for_bbox(day_utc=DAY, grid_size=2, **BBOX,
        fallback_values=DEFAULTS, log_to_console=True, search_days=0)
    assert len(calls) == 4 and {call[0] for call in calls} == {DAY}
    missing = codes(result, 'lookup_missing_values')
    assert len(missing) == 4 and missing[0]['parameters'] == ['aot550', 'ozone_cm_atm']
    assert result.aot550 == pytest.approx(.1) and result.water_vapor == 1
    assert result.defaults_used == {'ozone_cm_atm': .3}
    assert result.ozone_cm_atm == .3
    assert not codes(result, 'alternate_date_lookup')
    output = capsys.readouterr().out
    assert all('WARNING | ' + item['message'] in output for item in result.warnings)


@pytest.mark.parametrize('offset,sequence', [(1, [0, 1]), (-1, [0, 1, -1])])
def test_nasa_alternate_date_attempt_and_use_are_both_recorded(monkeypatch, offset, sequence, capsys):
    calls = []
    def lookup(day, *args, **kwargs):
        calls.append((day - DAY).days)
        return dict(DEFAULTS) if day == DAY + timedelta(days=offset) else dict.fromkeys(DEFAULTS)
    result = nasa(monkeypatch, lookup, fallback_values=DEFAULTS, log_to_console=True)
    assert calls == sequence
    assert len(codes(result, 'alternate_date_lookup')) == len(sequence) - 1
    assert len(codes(result, 'lookup_missing_values')) == len(sequence) - 1
    used = codes(result, 'alternate_date_used')
    assert len(used) == 1 and used[0]['date_used'] == (DAY + timedelta(days=offset)).isoformat()
    assert result.date_used == used[0]['date_used']
    assert result.date_requested == DAY.isoformat() and not result.defaults_used
    output = capsys.readouterr().out
    assert all(item['message'] in output for item in result.warnings)


def test_nasa_no_match_records_every_lookup_and_all_three_defaults(monkeypatch):
    result = nasa(monkeypatch, lambda *a, **kw: dict.fromkeys(DEFAULTS), fallback_values=DEFAULTS)
    assert len(result.lookup_attempts) == 3
    assert len(codes(result, 'lookup_missing_values')) == 3
    assert len(codes(result, 'alternate_date_lookup')) == 2
    assert len(codes(result, 'no_atmosphere_match')) == 1
    assert len(codes(result, 'default_used')) == 3
    assert result.defaults_used == DEFAULTS
    assert {k: getattr(result, k) for k in DEFAULTS} == DEFAULTS
    assert all(all(v is None for v in attempt['values'].values()) for attempt in result.lookup_attempts)


def test_nasa_no_configured_defaults_remains_null_with_explicit_warnings(monkeypatch, capsys):
    result = nasa(monkeypatch, lambda *a, **kw: dict.fromkeys(DEFAULTS), search_days=0)
    assert all(getattr(result, k) is None for k in DEFAULTS)
    assert len(codes(result, 'parameter_unavailable')) == 3
    assert result.defaults_used == {} and not codes(result, 'default_used')
    assert capsys.readouterr().out == ''


@pytest.mark.parametrize('values', [{'unknown': 1}, {'aot550': -1}, {'water_vapor': float('nan')},
                                  {'ozone_cm_atm': float('inf')}, {'aot550': True}])
def test_invalid_fallbacks_fail_before_any_lookup(monkeypatch, values):
    def unexpected(*args, **kwargs):
        pytest.fail('Must validate before making any requests')
    with pytest.raises(ValueError):
        nasa(monkeypatch, unexpected, fallback_values=values)


def test_modis_records_missing_terra_and_alternate_collection_band_and_day(monkeypatch, capsys):
    tomorrow = TIME + timedelta(days=1)
    def water(ee, collection, *args, **kwargs):
        if 'MOD08' in collection:
            return None
        return candidate(collection, timestamp=tomorrow)
    def aerosol(ee, collection, *args, **kwargs):
        if 'MOD08' in collection:
            return dict(collection=collection, band_found=False, bands=[], raw_value=None)
        return candidate(collection, parameter='aot550', band_index=1, timestamp=tomorrow)
    modis_mocks(monkeypatch, water, aerosol)
    result = modis(fallback_values=DEFAULTS, log_to_console=True)
    assert result.status == 'ok'
    assert result.water_vapor == 2 and result.aot550 == pytest.approx(.1)
    assert result.defaults_used == {'ozone_cm_atm': .3}
    assert len(codes(result, 'lookup_no_images')) == 1
    assert len(codes(result, 'band_not_found')) == 4
    assert len(codes(result, 'alternate_collection_used')) == 2
    assert len(codes(result, 'alternate_band_used')) == 1
    assert {w['parameter'] for w in codes(result, 'alternate_date_used')} == {'water_vapor', 'aot550'}
    assert all(w['date_used'] == tomorrow.date().isoformat() for w in codes(result, 'alternate_date_used'))
    assert len(result.lookup_attempts) == 4
    output = capsys.readouterr().out
    assert all('WARNING | ' + w['message'] in output for w in result.warnings)


def test_modis_same_utc_date_different_hour_is_not_an_alternate_day(monkeypatch):
    modis_mocks(monkeypatch, water=lambda ee, c, *a, **kw: candidate(c, timestamp=TIME - timedelta(hours=2)))
    result = modis(fallback_values=DEFAULTS)
    assert not codes(result, 'alternate_date_used')
    assert result.scene_datetime_utc == TIME.isoformat()
    assert result.defaults_used == {'ozone_cm_atm': .3}


def test_modis_aod_date_is_reported_independently_from_water_vapor(monkeypatch):
    modis_mocks(monkeypatch, aerosol=lambda ee, c, *a, **kw: candidate(
        c, parameter='aot550', timestamp=TIME - timedelta(days=1)))
    result = modis(fallback_values=DEFAULTS)
    records = codes(result, 'alternate_date_used')
    assert len(records) == 1 and records[0]['parameter'] == 'aot550'
    assert records[0]['date_used'] == (TIME - timedelta(days=1)).date().isoformat()


def test_modis_no_match_preserves_raw_missing_values_and_records_defaults(monkeypatch):
    def must_not_query(*args, **kwargs):
        pytest.fail('AOD lookup order must remain unchanged when water vapor is missing')
    modis_mocks(monkeypatch, water=lambda *a, **kw: None, aerosol=must_not_query)
    result = modis(fallback_values=DEFAULTS)
    assert result.status == 'no_modis_match'
    assert result.modis_scaled_g_cm2 is None and result.aod_scaled is None
    assert result.defaults_used == DEFAULTS
    assert len(codes(result, 'lookup_no_images')) == 2
    assert len(codes(result, 'lookup_skipped')) == 1
    assert len(codes(result, 'default_used')) == 3
    assert len(result.lookup_attempts) == 2


def test_modis_missing_regional_values_are_reported(monkeypatch):
    def water(ee, collection, *args, **kwargs):
        value = candidate(collection)
        value['raw_value'] = None
        return value
    modis_mocks(monkeypatch, water=water)
    result = modis()
    assert result.status == 'no_modis_match'
    assert len(codes(result, 'lookup_missing_values')) == 2
    assert {a['status'] for a in result.lookup_attempts} == {'no_valid_pixels'}
    assert result.water_vapor is None and result.defaults_used == {}


def test_modis_missing_timestamp_is_reported(monkeypatch):
    def water(ee, collection, *args, **kwargs):
        value = candidate(collection)
        value['abs_time_diff_hours'] = None
        return value
    modis_mocks(monkeypatch, water=water)
    result = modis()
    assert result.status == 'no_modis_match'
    assert {a['status'] for a in result.lookup_attempts} == {'missing_timestamp'}


def test_modis_missing_aod_is_defaulted_without_replacing_good_water_vapor(monkeypatch):
    modis_mocks(monkeypatch, aerosol=lambda *a, **kw: None)
    result = modis(fallback_values=DEFAULTS)
    assert result.water_vapor == result.modis_scaled_g_cm2 == 2
    assert result.aod_scaled is None and result.aot550 == .2
    assert result.defaults_used == {'aot550': .2, 'ozone_cm_atm': .3}
    assert {w['parameter'] for w in codes(result, 'lookup_no_images')} == {'aot550'}


@pytest.mark.parametrize('provider', ['nasa_power', 'modis_gee'])
def test_adapter_saved_json_contains_identical_warnings_and_effective_values(tmp_path, monkeypatch, capsys, provider):
    params = dict(source=provider, output_path=str(tmp_path/'atmosphere.json'), **BBOX,
                  fallback_values=DEFAULTS)
    if provider == 'nasa_power':
        params.update(day_utc=DAY.isoformat(), grid_size=1, search_days=0)
        monkeypatch.setattr(fetch, '_fetch_power_daily_point', lambda *a, **kw: dict.fromkeys(DEFAULTS))
    else:
        params.update(scene_datetime_utc=TIME.isoformat(), ee=object())
        modis_mocks(monkeypatch, water=lambda *a, **kw: None)
    result = fetch.FetchAtmosphere().run(params=params, shared={'log_to_console': True})
    saved = json.loads(Path(params['output_path']).read_text())
    assert saved == result
    assert saved['defaults_used'] == DEFAULTS
    assert saved['lookup_attempts'] and len([w for w in saved['warnings'] if w['code'] == 'default_used']) == 3
    output = capsys.readouterr().out
    assert all('WARNING | ' + w['message'] in output for w in saved['warnings'])
    # The example passes the persisted values to Py6S, with no hidden YAML substitution.
    config = yaml.safe_load((Path(__file__).parents[1]/'configs/example.worldview.yml').read_text())
    assert config['fetch_atmosphere']['param:fallback_values'] == DEFAULTS
    expected = dict(aot550=.2, water_vapor=2.5, ozone=.3)
    for param, value in expected.items():
        assert resolve(config['py6s_atmospheric_correction']['param:'+param],
            {'const': {}, 'var': {'atmosphere': saved}}) == value


@pytest.mark.parametrize('provider', ['nasa_power', 'modis_gee'])
def test_request_errors_are_logged_saved_and_reraised_without_defaulting(tmp_path, monkeypatch, capsys, provider):
    params = dict(source=provider, output_path=str(tmp_path/'atmosphere.json'), **BBOX,
                  fallback_values=DEFAULTS)
    calls = []
    failure = requests.Timeout('test timeout')
    def fail(*args, **kwargs):
        calls.append(1)
        raise failure
    if provider == 'nasa_power':
        params.update(day_utc=DAY.isoformat(), grid_size=1)
        monkeypatch.setattr(fetch, '_fetch_power_daily_point', fail)
    else:
        params.update(scene_datetime_utc=TIME.isoformat(), ee=object())
        monkeypatch.setattr(fetch, '_fetch_collection_value', fail)
    with pytest.raises(requests.Timeout) as caught:
        fetch.FetchAtmosphere().run(params=params, shared={'log_to_console': True})
    assert caught.value is failure and calls == [1]
    saved = json.loads(Path(params['output_path']).read_text())
    assert saved['status'] == 'error' and saved['source'] == provider
    assert saved['lookup_attempts'][0]['status'] == 'error'
    assert [w['code'] for w in saved['warnings']] == ['lookup_error']
    assert saved['warnings'][0]['message'] in capsys.readouterr().out


def test_earth_engine_initialization_failure_is_also_saved(tmp_path, monkeypatch):
    def fail(**kwargs):
        raise RuntimeError('no Earth Engine credentials')
    monkeypatch.setattr(fetch, '_init_ee_client', fail)
    path = tmp_path/'atmosphere.json'
    with pytest.raises(RuntimeError, match='credentials'):
        fetch.FetchAtmosphere().run(params=dict(source='modis_gee', output_path=str(path),
            scene_datetime_utc=TIME.isoformat(), **BBOX), shared={})
    saved = json.loads(path.read_text())
    assert saved['lookup_attempts'][0]['operation'] == 'initialize_earth_engine'
    assert saved['warnings'][0]['code'] == 'lookup_error'


def test_nasa_searches_missing_parameters_without_replacing_acquisition_day_values(monkeypatch):
    calls = []
    def lookup(day, *args, **kwargs):
        calls.append((day - DAY).days)
        if day == DAY:
            return dict(aot550=None, water_vapor=1., ozone_cm_atm=.28)
        return dict(aot550=.1, water_vapor=99., ozone_cm_atm=.99)
    result = nasa(monkeypatch, lookup, fallback_values=DEFAULTS)
    assert calls == [0, 1]
    assert (result.aot550, result.water_vapor, result.ozone_cm_atm) == (.1, 1., .28)
    assert not result.defaults_used
    assert result.parameter_dates == dict(aot550=(DAY + timedelta(days=1)).isoformat(),
                                         water_vapor=DAY.isoformat(), ozone_cm_atm=DAY.isoformat())
    assert codes(result, 'alternate_date_used')[0]['parameters'] == ['aot550']
    assert result.date_used == DAY.isoformat()


def test_nasa_defaults_only_after_searching_all_dates_for_missing_parameters(monkeypatch):
    result = nasa(monkeypatch, lambda *a, **kw: dict(aot550=None, water_vapor=1., ozone_cm_atm=.28),
                  fallback_values=DEFAULTS)
    assert len(result.lookup_attempts) == 3
    assert result.defaults_used == {'aot550': .2}
    assert result.parameter_dates['aot550'] is None
    assert result.parameter_dates['water_vapor'] == DAY.isoformat()
