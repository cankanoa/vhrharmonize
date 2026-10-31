"""Spectral integration and real-raster Py6S output validity."""
import importlib
from pathlib import Path
import sys
from types import SimpleNamespace as NS

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box
import yaml

from vhrharmonize.io.geospatial import get_image_percentile_value
from vhrharmonize.workflow.engine import Workflow
from workflow_helpers import import_settings

atmos = importlib.import_module('vhrharmonize.plugins.atmospheric_correction')


@pytest.fixture
def model(monkeypatch):
    state = NS(coefficients=(1., -.1, 0.), wavelengths=[])
    class SixS:
        def __init__(self, **kwargs):
            self.altitudes = NS(set_sensor_satellite_level=lambda:None,
                                set_target_custom_altitude=lambda x:None)
        def run(self):
            state.wavelengths.append(self.wavelength)
            self.outputs = NS(**dict(zip(('coef_xa','coef_xb','coef_xc'),state.coefficients)))
    fake = NS(SixS=SixS, Geometry=NS(User=lambda:NS()),
              AeroProfile=NS(Maritime=1,Continental=2,Urban=3,Desert=4,BiomassBurning=5,
                             Stratospheric=6,PredefinedType=lambda x:x),
              AtmosProfile=NS(UserWaterAndOzone=lambda *x:x),
              AtmosCorr=NS(AtmosCorrLambertianFromReflectance=lambda x:x),
              Wavelength=lambda *args:args)
    monkeypatch.setitem(sys.modules,'Py6S',fake)
    return state


def correct(source, output, **kwargs):
    params = dict(solar_zenith=30,solar_azimuth=120,view_zenith=5,view_azimuth=80,
                  day=1,month=6,ground_elevation_km=0,atmosphere_profile='user',
                  wavelengths_center_um=[.55],sixs_executable='mocked')
    params.update(kwargs)
    return atmos.Py6SCorrector().run(str(source),str(output),**params)


def write_raster(path, values, *, nodata=None, mask=None):
    values=np.asarray(values,dtype='float32')
    if values.ndim==2:values=values[None]
    with rasterio.open(path,'w',driver='GTiff',width=values.shape[2],height=values.shape[1],
                       count=len(values),dtype='float32',crs='EPSG:4326',
                       transform=from_origin(0,values.shape[1],1,1),nodata=nodata) as ds:
        ds.write(values)
        if mask is not None:ds.write_mask(np.asarray(mask,dtype='uint8'))
    return path


@pytest.mark.parametrize('params,args',[
    ({'wavelengths_center_um':.55},[(.55,)]),
    ({'wavelengths_center_um':[.55,.66]},[(.55,),(.66,)]),
    ({'wavelengths_bounds_um':[.4,.6]},[(.4,.6)]),
    ({'wavelengths_bounds_um':[[.4,.6],[.7,.8]]},[(.4,.6),(.7,.8)]),
])
def test_point_and_bounds_map_to_native_wavelength(params,args):
    assert atmos._py6s_spectral_bands(params)==args


def response_csv(tmp_path, text='wavelength_nm,BLUE\n400,0\n405,1\n410,0\n'):
    path=tmp_path/'response.csv';path.write_text(text)
    return str(path)


def test_csv_interpolates_actual_irregular_wavelengths_without_modifying_source(tmp_path):
    path=response_csv(tmp_path,'wavelength_nm,BLUE,GREEN\n400,0,2\n402,1,0\n410,0,0\n')
    original=Path(path).read_bytes()
    spectra=atmos._py6s_spectral_bands({'wavelengths_response_csv':path})
    # At 402.5 nm the BLUE response is .9375; normalize its sampled peak.
    np.testing.assert_allclose(spectra[0][2],np.array([0,.9375,.625,.3125,0])/.9375)
    np.testing.assert_allclose(spectra[1][2],[1,0,0,0,0])
    assert spectra[0][:2]==pytest.approx((.4,.41))
    assert Path(path).read_bytes()==original


def test_csv_pads_to_native_grid_and_normalizes_at_runtime(tmp_path):
    path=response_csv(tmp_path,'wavelength_nm,BLUE\n401,2\n404,2\n')
    spectra=atmos._py6s_spectral_bands({'wavelengths_response_csv':path})
    assert spectra[0][:2]==pytest.approx((.4,.405))
    assert spectra[0][2]==[0,1,0]


@pytest.mark.parametrize('params',[
    {'wavelengths_center_um':True}, {'wavelengths_center_um':'0.55'},
    {'wavelengths_center_um':.1}, {'wavelengths_center_um':float('nan')},
    {'wavelengths_center_um':[]}, {'wavelengths_center_um':[[.5,.6]]},
    {'wavelengths_bounds_um':[.6,.4]}, {'wavelengths_bounds_um':[.4,.4]},
    {'wavelengths_bounds_um':[.4]}, {'wavelengths_bounds_um':[]},
])
def test_invalid_centers_or_bounds_are_rejected(params):
    with pytest.raises(ValueError):atmos._py6s_spectral_bands(params)


@pytest.mark.parametrize('text',[
    '', 'wavelength,BLUE\n400,0\n405,1\n',
    'wavelength_nm,BLUE,BLUE\n400,0,0\n405,1,1\n',
    'wavelength_nm,\n400,0\n405,1\n',
    'wavelength_nm,BLUE\n400,1\n',
    'wavelength_nm,BLUE\n400,1\n400,1\n',
    'wavelength_nm,BLUE\n410,1\n400,1\n',
    'wavelength_nm,BLUE\n100,1\n405,1\n',
    'wavelength_nm,BLUE\n400,1\n4001,1\n',
    'wavelength_nm,BLUE\n400,0\n405,0\n',
    'wavelength_nm,BLUE\n400,-1\n405,1\n',
    'wavelength_nm,BLUE\n400,nan\n405,1\n',
    'wavelength_nm,BLUE\n400,1\ninf,1\n',
    'wavelength_nm,BLUE\n400,\n405,1\n',
    'wavelength_nm,BLUE\n400,0,1\n405,1\n',
])
def test_invalid_csvs_are_rejected(tmp_path,text):
    with pytest.raises(ValueError):
        atmos._py6s_spectral_bands({'wavelengths_response_csv':response_csv(tmp_path,text)})


@pytest.mark.parametrize('selected',[(),(0,1),(0,2),(1,2),(0,1,2)])
def test_exactly_one_representation_is_required_in_wrapper_and_adapter(tmp_path,model,selected):
    values=([.55],[[.5,.6]],str(tmp_path/'not-read.csv'))
    params=dict(zip(atmos._PY6S_SPECTRAL_PARAMETERS,[None]*3))
    params.update({atmos._PY6S_SPECTRAL_PARAMETERS[i]:values[i] for i in selected})
    with pytest.raises(ValueError,match='exactly one'):
        correct('not-read.tif',tmp_path/'out.tif',**params)
    with pytest.raises(ValueError,match='exactly one'):
        atmos.run_py6s('not-read.tif',str(tmp_path/'out.tif'),**params)


def test_exclusivity_is_checked_after_mapped_variable_fallback(tmp_path):
    with pytest.raises(ValueError,match='exactly one'):
        atmos.run_py6s('not-read.tif',str(tmp_path/'out.tif'),wavelengths_center_um=[.55],
                      variables={'wavelengths_response_csv':'not-read.csv'})


def test_csv_can_be_supplied_through_mapped_variables(tmp_path,model):
    source=write_raster(tmp_path/'in.tif',[[.1,.2],[.3,.4]])
    path=response_csv(tmp_path)
    result=atmos.run_py6s(str(source),str(tmp_path/'out.tif'),
        variables=dict(solar_zenith=30,solar_azimuth=120,view_zenith=5,view_azimuth=80,
                       day=1,month=6,wavelengths_response_csv=path),
        atmosphere_profile='user',sixs_executable='mocked')
    assert result.effective_params['wavelengths_response_csv']==path
    assert len(model.wavelengths)==1


@pytest.mark.parametrize('params',[
    {'wavelengths_center_um':[.5,.6]},
    {'wavelengths_bounds_um':[[.4,.5],[.6,.7]]},
    {'wavelengths_response_csv':'csv'},
])
def test_spectral_band_count_must_match_raster_exactly(tmp_path,model,params):
    if 'wavelengths_response_csv' in params:
        params={'wavelengths_response_csv':response_csv(tmp_path,
            'wavelength_nm,BLUE,GREEN\n400,0,1\n405,1,0\n')}
    source=write_raster(tmp_path/'in.tif',[[.1,.2],[.3,.4]])
    with pytest.raises(ValueError,match='band count exactly'):
        correct(source,tmp_path/'out.tif',**{'wavelengths_center_um':None,**params})
    assert not model.wavelengths


def test_mask_only_input_and_nan_values_remain_invalid(tmp_path,model):
    source=write_raster(tmp_path/'in.tif',[[.2,np.nan],[.4,.5]],mask=[[0,255],[255,255]])
    output=tmp_path/'out.tif';correct(source,output)
    with rasterio.open(output) as ds:
        assert np.isnan(ds.nodata)
        assert ds.read_masks(1).tolist()==[[0,0],[255,255]]
        np.testing.assert_allclose(ds.read(1)[1],[.5,.6])


def test_output_nodata_does_not_mask_a_valid_input_dn(tmp_path,model):
    source=write_raster(tmp_path/'in.tif',[[0,.2],[.3,-9999]],nodata=-9999)
    output=tmp_path/'out.tif';correct(source,output,custom_nodata_value=0)
    with rasterio.open(output) as ds:
        assert ds.read(1)[0,0]==pytest.approx(.1)
        assert ds.read_masks(1)[0,0]==255 and ds.read_masks(1)[1,1]==0


def test_inherited_zero_nodata_does_not_hide_valid_clipped_zero(tmp_path,model):
    model.coefficients=(1.,.5,0.)
    source=write_raster(tmp_path/'in.tif',[[0,.2],[.7,.8]],nodata=0)
    output=tmp_path/'out.tif';correct(source,output)
    with rasterio.open(output) as ds:
        assert np.isnan(ds.nodata)
        assert ds.read(1)[0,1]==0
        assert ds.read_masks(1)[0,1]==255 and ds.read_masks(1)[0,0]==0


def test_invalid_division_is_masked_before_clipping(tmp_path,model):
    model.coefficients=(1.,0.,1.)
    source=write_raster(tmp_path/'in.tif',[[-1,-2],[.2,.3]])
    output=tmp_path/'out.tif';correct(source,output)
    with rasterio.open(output) as ds:
        assert ds.read_masks(1).tolist()==[[0,0],[255,255]]


@pytest.mark.parametrize('scale',[None,1,0,float('nan')])
def test_integer_output_rejects_missing_or_invalid_scale(tmp_path,model,scale):
    source=write_raster(tmp_path/'in.tif',[[.1,.2],[.3,.4]])
    with pytest.raises(ValueError,match='scale|output_scale_factor'):
        correct(source,tmp_path/'out.tif',output_dtype='int16',output_scale_factor=scale)


def test_scaled_integer_rounds_preserves_nodata_and_records_physical_units(tmp_path,model):
    source=write_raster(tmp_path/'in.tif',[[.123456,.2],[.3,.4]],mask=[[255,0],[255,255]])
    output=tmp_path/'out.tif';correct(source,output,output_dtype='int16',output_scale_factor=10000)
    with rasterio.open(output) as ds:
        assert ds.read(1)[0,0]==2235
        assert ds.read_masks(1)[0,1]==0
        assert ds.scales==(0.0001,) and ds.offsets==(0.,)
        assert ds.units==('reflectance',)
        assert ds.tags()['ATMOSPHERIC_CORRECTION']=='Py6S'


@pytest.mark.parametrize('kwargs,message',[
    ({'output_dtype':'uint8','output_scale_factor':10000,'custom_nodata_value':255},'exceeds'),
    ({'output_dtype':'uint16','output_scale_factor':10000,'custom_nodata_value':-9999},'representable'),
    ({'ground_elevation_km':float('nan')},'finite'),
    ({'dn_to_radiance_factors':[float('nan')]},'finite'),
])
def test_invalid_output_or_model_settings_fail_clearly(tmp_path,model,kwargs,message):
    source=write_raster(tmp_path/'in.tif',[[.1,.2],[.3,.4]])
    with pytest.raises(ValueError,match=message):correct(source,tmp_path/'out.tif',**kwargs)


def test_nonfinite_coefficients_fail_before_output_is_created(tmp_path,model):
    model.coefficients=(float('nan'),0,0)
    source=write_raster(tmp_path/'in.tif',[[.1,.2],[.3,.4]])
    output=tmp_path/'out.tif'
    with pytest.raises(ValueError,match='nonfinite coefficients'):correct(source,output)
    assert not output.exists()


def test_dem_sampling_crops_to_the_footprint_without_changing_default_callers(tmp_path,monkeypatch):
    source=write_raster(tmp_path/'dem.tif',np.arange(100).reshape(10,10))
    original=rasterio.mask.mask;shapes=[]
    def record(*args,**kwargs):
        result=original(*args,**kwargs);shapes.append(result[0].shape)
        return result
    monkeypatch.setattr(rasterio.mask,'mask',record)
    assert get_image_percentile_value(str(source),mask=box(0,8,2,10),crop=True)==5.5
    assert get_image_percentile_value(str(source),mask=box(0,8,2,10))==5.5
    assert shapes==[(2,2),(10,10)]


@pytest.mark.parametrize('sensor,count,upper',[('WV02',751,1100),('WV03',2101,2450)])
def test_example_selects_original_sensor_csv_in_raster_band_order(tmp_path,make_test_raster,sensor,count,upper):
    import csv
    config_dir=Path(__file__).parents[1]/'configs'
    recipe=yaml.safe_load((config_dir/'example.worldview.yml').read_text())
    discover=import_settings(make_test_raster(tmp_path/'in.tif',count=8),tmp_path)
    discover['const:spectral_response_csv']=recipe['discover_inputs']['const:spectral_response_csv']
    discover['var:sensor_id']=sensor
    discover['var:wavelengths_response_csv']=recipe['discover_inputs']['var:wavelengths_response_csv']
    discover['var:current_image_paths']='returned:file_path'
    workflow=Workflow({'shared':{'plugin':'shared','core:run':True,'core:save_statistics_path':None},
        'discover':discover,'correct':{'plugin':'atmospheric_correction','core:run':True,
            'core:require_outputs':'param:output_raster','param:input_raster':'var:current_image_paths',
            'param:output_raster':str(tmp_path/'corrected.tif'),
            'param:wavelengths_response_csv':'var:wavelengths_response_csv'}},config_dir=config_dir)
    node=next(n for n in workflow.plan().nodes if n.step['name']=='correct')
    csv_path=Path(node.params['wavelengths_response_csv'])
    assert csv_path==config_dir/f'example.worldview{sensor[-1]}.response.csv'
    with csv_path.open() as handle:
        rows=list(csv.reader(handle))
    assert rows[0]==['wavelength_nm',*recipe['discover_inputs']['const:band_order']]
    assert len(rows)==count+1
    assert [int(row[0]) for row in rows[1:]]==list(range(350,upper+1))
    assert len(atmos._py6s_spectral_bands(node.params))==8


@pytest.mark.parametrize('flag,raw',[
    ('wavelengths-center-um','0.55'),('wavelengths-center-um','[0.55, 0.66]'),
    ('wavelengths-bounds-um','[[0.4, 0.6], [0.7, 0.8]]'),
    ('wavelengths-response-csv','response.csv'),
])
def test_cli_decodes_new_spectral_parameters(flag,raw):
    import json
    from vhrharmonize.cli.functions import build_parser
    args=build_parser(atmos.run_py6s).parse_args(['--'+flag,raw])
    expected=raw if flag.endswith('csv') else json.loads(raw)
    assert getattr(args,flag.replace('-','_'))==expected


def test_hpc_stages_and_remaps_response_csv_as_a_required_input(tmp_path,model):
    from vhrharmonize.workflow.staging import stage_workflow
    from workflow_helpers import transfer
    source=write_raster(tmp_path/'in.tif',[[.1,.2],[.3,.4]])
    write_raster(tmp_path/'other.tif',[[.2,.3],[.4,.5]])
    path=response_csv(tmp_path)
    discover=import_settings(tmp_path/'*.tif',tmp_path)
    discover['var:wavelengths_response_csv']='path:./response.csv'
    recipe={'shared':{'plugin':'shared','core:run':True,'core:save_statistics_path':None,
                      'core:load_statistics_path':None,'const:root':str(tmp_path)},
        'discover':discover,'correct':{'plugin':'atmospheric_correction','core:run':True,
        'core:require_outputs':'param:output_raster','param:input_raster':'var:mul',
        'param:output_raster':"expr:const.root & '/' & var.basename & '_out.tif'",
        'param:wavelengths_response_csv':'var:wavelengths_response_csv',
        'param:solar_zenith':30,'param:solar_azimuth':120,'param:view_zenith':5,
        'param:view_azimuth':80,'param:day':1,'param:month':6,
        'param:atmosphere_profile':'user','param:sixs_executable':'mocked'}}
    remote=tmp_path/'remote'
    staged,uploads,_=stage_workflow(recipe,config_dir=tmp_path,remote_work_dir=str(remote),
        path_mappings={'const:root':str(remote),'var:wavelengths_response_csv':str(remote/'references')})
    assert uploads[path]==str(remote/'references/response.csv')
    assert sum(key==path for key in uploads)==1
    transfer(uploads)
    workflow=Workflow(staged,config_dir=tmp_path).plan()
    nodes=[n for n in workflow.nodes if n.step['name']=='correct']
    assert len(nodes)==2
    assert all(n.params['wavelengths_response_csv']==uploads[path] for n in nodes)
    workflow.run()
    assert (remote/'in_out.tif').exists() and (remote/'other_out.tif').exists()
