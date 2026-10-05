"""Analytical, file-format, and regression tests. Run python -m pytest -q."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import erf

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import coq_core as c


@pytest.mark.parametrize("model,params",[("Linear",(100,1500,1)),("Log(K)",(1000,120,1)),
                                         ("Square root",(50,700,1)),("Log(1+cK)",(200,900,.3))])
def test_response_round_trip_and_fit(model,params):
    k=np.array([.5,1,2,4,8,12])
    original=c.Response(model,*params)
    fit=c.fit_response(k,original.forward(k),model)
    np.testing.assert_allclose(fit.inverse(original.forward(k)),k,rtol=1e-5)
    assert fit.r2> .999999


def test_invalid_response_and_roi():
    with pytest.raises(ValueError):c.fit_response([1,2,3,4],[1,1,1,1])
    with pytest.raises(ValueError):c.Response("Linear",0,0).inverse([1])
    with pytest.raises(ValueError):c.crop(np.ones((10,10)),(8,8,4,4))
    with pytest.raises(ValueError):c.Response("Square root",10,1).inverse([9])


@pytest.mark.parametrize("method",["edge","slit"])
@pytest.mark.parametrize("horizontal",[False,True])
def test_gaussian_mtf(method,horizontal):
    y,x=np.indices((512,512));sigma=1.1;pitch=.2
    distance=(x-256-.055*(y-256))/np.sqrt(1+.055**2)
    if method=="edge":image=1000+9000*.5*(1+erf(distance/(np.sqrt(2)*sigma)))
    else:image=1000+9000*np.exp(-.5*(distance/sigma)**2)
    if horizontal:image=image.T
    result=c.mtf(image,pitch,(128,128,256,256),method=method)
    expected50=np.sqrt(2*np.log(2))/(2*np.pi*sigma*pitch)
    expected10=np.sqrt(2*np.log(10))/(2*np.pi*sigma*pitch)
    assert result["mtf"][0]==pytest.approx(1)
    assert result["mtf50_lp_mm"]==pytest.approx(expected50,rel=.025)
    assert result["mtf10_lp_mm"]==pytest.approx(expected10,rel=.025)
    assert result["orientation"]==("horizontal" if horizontal else "vertical")
    assert result["frequency"].max()<=.5/pitch


def test_non_slanted_edge_is_rejected():
    y,x=np.indices((128,128))
    with pytest.raises(ValueError,match="angle"):c.mtf((x>64)*1000+100,.2)


@pytest.mark.parametrize("window",["none","hann"])
def test_white_noise_nps_normalization(window):
    rng=np.random.default_rng(1234);mean=10000;sd=100;pitch=.2
    arrays=[mean+rng.normal(0,sd,(512,512)) for _ in range(6)]
    result=c.nps(arrays,pitch,fragment=128,window=window)
    expected=sd**2*pitch**2/mean**2
    assert result["nnps2d"].mean()==pytest.approx(expected,rel=.025)
    assert np.median(result["horizontal"])==pytest.approx(expected,rel=.12)
    # Integral of the two-sided NPS recovers variance for an unwindowed tile.
    if window=="none":
        df=1/(128*pitch)
        assert result["nps2d"].sum()*df**2==pytest.approx(sd**2,rel=.025)
    assert result["accepted_fragments"]==6*49


def test_polynomial_detrending_removes_gradient():
    y,x=np.indices((128,128));image=10000+2*x+3*y+.01*x*x
    result=c.nps([image],.2,fragment=128,detrend="quadratic",reject_coq=False)
    assert result["nnps2d"].max()<1e-20


def test_nps_rejects_invalid_inputs():
    with pytest.raises(ValueError):c.nps([np.zeros((64,64))],.2,fragment=64)
    with pytest.raises(ValueError):c.nps([np.ones((64,64))],(.2,.3),fragment=64)
    with pytest.raises(ValueError):c.nps([np.ones((64,64))],.2,fragment=128)


def test_dqe_known_and_overlap():
    f=np.linspace(0,2,51);m=np.exp(-f);nf=np.linspace(.1,2.5,49);nnps=np.full(len(nf),1e-5)
    result=c.dqe(f,m,nf,nnps,2.5,30174)
    expected=np.interp(result["frequency"],f,m)**2/(1e-5*2.5*30174)
    np.testing.assert_allclose(result["dqe"],expected)
    assert result["points_outside_overlap"]>0
    with pytest.raises(ValueError):c.dqe(f,m,nf,nnps,0,30174)
    with pytest.raises(ValueError):c.dqe(f,m,nf,-nnps,1,30174)


def test_bad_pixel_coordinates_and_partial_tiles():
    a=np.ones((57,63))*1000;a[10,20]=1600;a[30,40]=500
    result=c.bad_pixels(a,1,tile_mm=10,margin_mm=0)
    assert result["count"]==2
    assert set(map(tuple,result["coordinates_xy"]))=={(20,10),(40,30)}
    assert result["tested_pixels"]==a.size


def test_uniformity_global_and_local_absolute_max():
    rng=np.random.default_rng(5);a=1000+rng.normal(0,5,(180,180))
    a[60:120,60:120]*=.8
    r=c.uniformity(a,1,tile_mm=30,overlap=0)
    means=r["signal_map"]
    assert r["GSNU_percent"]==pytest.approx(200*(means.max()-means.min())/(means.max()+means.min()))
    vals=[]
    for y in range(1,5):
        for x in range(1,5):
            neighbor=(means[y-1:y+2,x-1:x+2].sum()-means[y,x])/8
            vals.append(abs((means[y,x]-neighbor)/neighbor)*100)
    assert r["LSNU_percent"]==pytest.approx(max(vals))
    assert r["LSNU_percent"]>=0


def test_dark_five_rois():
    a=np.arange(200*200,dtype=float).reshape(200,200)
    r=c.dark_analysis(a,1,tile_mm=20,margin_mm=10)
    assert len(r["rois"])==5
    central=a[90:110,90:110]
    assert r["rois"][0]["mean"]==pytest.approx(central.mean())
    assert r["rois"][0]["sd"]==pytest.approx(central.std(ddof=1))


def test_cnr_background_pooling():
    rng=np.random.default_rng(6);a=1000+rng.normal(0,10,(100,100));a[10:30,10:30]-=100
    sr=(10,10,20,20);bgs=[(40,10,20,20),(70,10,20,20),(40,50,20,20),(70,50,20,20)]
    r=c.cnr(a,sr,bgs);bg=np.concatenate([c.crop(a,b).ravel() for b in bgs])
    assert r["cnr"]==pytest.approx(abs(bg.mean()-c.crop(a,sr).mean())/bg.std(ddof=1))
    with pytest.raises(ValueError):c.cnr(a,sr,[sr])


def test_lag_equations_and_order():
    inside=(5,5,20,20);outside=(60,5,20,20)
    y,x=np.indices((100,100));test=np.where(x<50,1000.,10000.)
    dark=np.where(x<50,502.,500.)
    assert c.lag([test,dark],inside,outside)["lag_percent"]==pytest.approx(.02)
    before=np.full((100,100),10000.);after=np.where(x<50,9970.,10000.)
    assert c.lag([test,before,after],inside,outside,"multiplicative")["lag_percent"]==pytest.approx(.3)


def test_agd_legacy_coefficients_manual_and_inverse_square():
    corrected=c.agd(45,.36,"Mo/Mo",5)
    assert corrected["g"]==pytest.approx(.11*.36**2+.38*.36+.031)
    assert corrected["c"]==pytest.approx(-.024*.36**2-.048*.36+1.1)
    assert corrected["agd_mGy"]==pytest.approx(.982747026688)
    legacy=c.agd(45,.36,"Mo/Mo",5,legacy_index=True)
    assert legacy["g"]==pytest.approx(.055*.36**2+.47*.36+.034)
    assert legacy["agd_mGy"]!=pytest.approx(corrected["agd_mGy"])
    manual=c.agd(25,.4,"Custom",mas=100,output_uGy_mAs=156,reference_distance_mm=500,entrance_distance_mm=600,factors=(.25,1.1,1))
    assert manual["agd_mGy"]==pytest.approx(100*156/1000*(500/600)**2*.25*1.1)
    with pytest.raises(ValueError):c.agd(45,.8,"Mo/Mo",5)


def test_tube_table_and_unavailable_zeros():
    table=c.AGDTables()
    assert table.lookup("Mo/Mo",28)==pytest.approx((.36,156))
    assert table.lookup("Mo/Mo",28.5)==pytest.approx((.365,165))
    with pytest.raises(ValueError):table.lookup("W/Al",28)
    with pytest.raises(ValueError):table.lookup("Mo/Mo",35)


def test_raw_endian_offset_gap(tmp_path):
    a=np.arange(24,dtype=">u2").reshape(4,6);b=a+10
    path=tmp_path/"images.raw";path.write_bytes(b"HEADER"+a.tobytes()+b"xx"+b.astype(">u2").tobytes())
    records=c.load_raw(path,4,6,"uint16",">",offset=6,frames=2,gap=2)
    np.testing.assert_array_equal(records[0].pixels,a)
    np.testing.assert_array_equal(records[1].pixels,b)
    with pytest.raises(ValueError):c.load_raw(path,4,6,"uint16",">",offset=6,frames=3,gap=2)


def test_dicom_spacing_signal_and_dose_units(tmp_path):
    from pydicom.dataset import Dataset,FileDataset,FileMetaDataset
    from pydicom.uid import ExplicitVRLittleEndian,generate_uid,DigitalXRayImageStorageForProcessing
    meta=FileMetaDataset();meta.TransferSyntaxUID=ExplicitVRLittleEndian
    meta.MediaStorageSOPClassUID=DigitalXRayImageStorageForProcessing;meta.MediaStorageSOPInstanceUID=generate_uid()
    ds=FileDataset(str(tmp_path/"qc.dcm"),{},file_meta=meta,preamble=b"\0"*128)
    ds.SOPClassUID=meta.MediaStorageSOPClassUID;ds.SOPInstanceUID=meta.MediaStorageSOPInstanceUID
    ds.Rows=4;ds.Columns=6;ds.SamplesPerPixel=1;ds.PhotometricInterpretation="MONOCHROME1"
    ds.BitsAllocated=16;ds.BitsStored=16;ds.HighBit=15;ds.PixelRepresentation=0
    ds.ImagerPixelSpacing=[.2,.3];ds.PixelSpacing=[.1,.15];ds.RescaleSlope=2;ds.RescaleIntercept=10
    ds.OrganDose=.025;ds.EntranceDoseInmGy=5
    a=np.arange(24,dtype=np.uint16).reshape(4,6);ds.PixelData=a.tobytes();ds.save_as(tmp_path/"qc.dcm",enforce_file_format=True)
    rec=c.load_images(tmp_path/"qc.dcm")[0]
    assert rec.spacing==(.2,.3)
    np.testing.assert_array_equal(rec.pixels,a)
    assert rec.metadata["organ_dose_mGy"]==pytest.approx(2.5)
    rescaled=c.load_images(tmp_path/"qc.dcm",True)[0]
    np.testing.assert_array_equal(rescaled.pixels,a*2+10)


def test_tiff_npy_and_export_round_trip(tmp_path):
    import tifffile
    a=np.arange(100,dtype=np.float32).reshape(10,10)
    tifffile.imwrite(tmp_path/"test.tif",a);np.save(tmp_path/"test.npy",a)
    np.testing.assert_array_equal(c.load_images(tmp_path/"test.tif")[0].pixels,a)
    np.testing.assert_array_equal(c.load_images(tmp_path/"test.npy")[0].pixels,a)
    f=np.linspace(0,2,100);r={"frequency":f,"mtf":np.exp(-f)}
    c.export_result(tmp_path/"mtf.csv","MTF",r,{"test":True})
    rf,rv=c.read_curve(tmp_path/"mtf.csv","mtf")
    np.testing.assert_array_equal(rf,f);np.testing.assert_array_equal(rv,r["mtf"])
    c.export_result(tmp_path/"mtf.json","MTF",r)
    assert json.loads((tmp_path/"mtf.json").read_text())["result"]["mtf"][0]==1


def test_nnps_export_import_explicit_column(tmp_path):
    f=np.arange(1,5,dtype=float);r={"frequency":f,"horizontal":f*1e-5,"vertical":f*2e-5,"average":f*1.5e-5,"radial":f*1.4e-5}
    c.export_result(tmp_path/"nnps.csv","NNPS",r)
    rf,rv=c.read_curve(tmp_path/"nnps.csv","average")
    np.testing.assert_array_equal(rv,r["average"])
