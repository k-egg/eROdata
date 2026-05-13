from make_psfmap import *
from make_exposure_1 import *
from make_exposure_2 import *
from SherpaSpectralModel import *

#from SkyModel_unit_fix import *

import gammapy
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from regions import SkyRegion,RectangleSkyRegion, Regions, CircleSkyRegion, EllipseSkyRegion

# %matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl
from IPython.display import display
from gammapy.data import EventList, Observation, GTI, FixedPointingInfo
from gammapy.datasets import Datasets, MapDataset, MapDatasetOnOff, SpectrumDatasetOnOff, SpectrumDataset
from gammapy.irf import EDispKernelMap, PSFMap,EDispKernel, PSF3D
from gammapy.maps import Map, MapAxis, WcsGeom,RegionNDMap,RegionGeom,MapAxes,WcsNDMap
from gammapy.estimators import FluxPointsEstimator, ExcessMapEstimator
from gammapy.makers import MapDatasetMaker, utils
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    Models,
    PointSpatialModel,
    PowerLawNormSpectralModel,
    PowerLawSpectralModel,
    SkyModel,
    TemplateSpatialModel,
    TemplateSpectralModel,
    create_fermi_isotropic_diffuse_model,
    GaussianSpatialModel,
    GaussianSpectralModel,
    EBLAbsorptionNormSpectralModel,
    ConstantSpatialModel
)
from gammapy.utils.scripts import make_name, make_path
import numpy as np
from os import path
import os, os.path, time, subprocess
import glob
import math


import astropy.units as u
import numpy as np
from gammapy.modeling.models import SpectralModel
from gammapy.modeling import Parameter, Parameters


def Evtfile_converter(evtfile,TM, add_pointing=None, suffix=None):
    if suffix==None:
        suffix=""
    # create output filename
    outfile=evtfile.replace(".fits","")+suffix+".fits"
    # read GTI & events
    hdulist = fits.open(evtfile)
    hdu=hdulist["GTI"+str(TM)]
    hdu1=hdulist["EVENTS"]

    # change EVENTS extension keywords
    hdu1.header["MJDREFI"]=math.modf(hdu.header["MJDREF"])[1]
    hdu1.header["MJDREFF"]=math.modf(hdu.header["MJDREF"])[0]
    hdu1.header["GEOLON"]=0
    hdu1.header["GEOLAT"]=0
    hdu1.header["ALTITUDE"]=0
    hdu1.columns["PI"].name = "ENERGY"
    hdu1.columns["ENERGY"].unit = "eV"
    
    # change GTI extension keywords
    hdu.name="GTI"
    hdu.header["MJDREFI"]=math.modf(hdu.header["MJDREF"])[1]
    hdu.header["MJDREFF"]=math.modf(hdu.header["MJDREF"])[0]
    hdu.header["RADECSYS"]="ICRS"
    hdu.header["RADECSYSa"]="ICRS"

    # add pointing direction if necessary
    if add_pointing:
        hdu1.header["RA_PNT"]=add_pointing[0]
        hdu1.header["DEC_PNT"]=add_pointing[1]

    # get average deadtime correction & write in header
    hdu2=hdulist["DEADCOR"+str(TM)]
    dc=[hdu2.data[i][1] for i in range(len(hdu2.data))]
    deadc=np.average(dc)
    hdu1.header["DEADC"]=deadc

    # write output
    hdulist.writeto(outfile,overwrite=True)
    print("Eventfile successfully converted!")

def RMF_converter(rmf_file, suffix=None):
    # open file and get matrix and energy bounds data
    rmf_hdulist = fits.open(rmf_file)
    
    rmf_matrix=rmf_hdulist[1]
    rmf_matrix_data=rmf_matrix.data
    rmf_ebounds=rmf_hdulist["EBOUNDS"]
    ebounds=rmf_ebounds.data

    # lower F_CHAN column values by 1
    for i in range(len(rmf_matrix_data["F_CHAN"])):
        k=rmf_matrix_data["F_CHAN"][i]
        rmf_matrix_data["F_CHAN"][i]=[k[0]-1]

    # write output
    outfile_rmf=rmf_file.replace(".fits","")+suffix+".fits"
    fits.writeto(outfile_rmf,data=rmf_matrix_data,overwrite=True,header=rmf1.header)
    fits.append(outfile_rmf,data=ebounds,header=rmf_ebounds.header)
    print("RMF successfully converted!")


class eROdata:
    '''
    Container for creating eROSITA dataset in Gammapy.
    '''
    def __init__(self,txt_file,reg_file,out_path,tms=[1,2,3,4,6],pointed=False,catalog_path=""):
        self.txt_file=txt_file
        #allow for both region strings and region files
        try:
            ds9_reg=open(reg_file,"r").read()
        except:
            ds9_reg=reg_file
        self.ds9_reg=ds9_reg
        self.reg_coords=[float(i) for i in ds9_reg.split("(")[1].split(")")[0].split(",")]

        self.region=Regions.parse(ds9_reg, format="ds9")[0]

        
        if out_path[-1]!="/":
            out_path=out_path+"/"
        self.out_path=out_path
        self.tms=tms
        self.pointed=pointed
        self.catalog_path=catalog_path

    def get_eRASS1_point_src(self):
        try:
            ds9_reg=Regions.read(self.out_path+"point_src_list.reg", format='ds9')
            return ds9_reg
        except:
            erass = Table(fits.open(self.catalog_path+"eRASS1_Main.v1.1.fits")[1].data)
            catalog=SkyCoord(ra=erass["RA"].data*u.deg,dec=erass["DEC"].data*u.deg)
            c = SkyCoord(ra=[self.reg_coords[0]]*u.degree, dec=[self.reg_coords[1]]*u.degree)
            idx1, idx2, sep, dist=catalog.search_around_sky(c,self.reg_coords[2]*u.deg)
            regions=Regions([CircleSkyRegion(center=i,radius=30*u.arcsec) for i in catalog[idx2]])
            regions.write(self.out_path+"point_src_list.reg", overwrite=True)
            return regions
    
    def eSASS_data_products(self):
        # generate eventfiles and expmaps for all TMs
        txt_file=self.txt_file
        reg=self.ds9_reg
        out_path=self.out_path
        self.evtfiles=[]
        self.evtfiles_conv=[]
        self.expmaps=[]

        self.get_eRASS1_point_src()

        if not self.pointed:
            pnt=self.reg_coords
        else:
            pnt=None
        
        x="auto"
        for n,i in enumerate(self.tms):
            self.evtfiles.append(out_path+"evt_TM"+str(i)+".fits")
            cmd=["evtool", "eventfiles=@"+txt_file, "outfile= "+out_path+"evt_TM"+str(i)+".fits", "image=yes", "rebin=80","pattern=15", "size="+str(x), "gti=FLAREGTI", "center_position= auto","flag=0xc00f7f30","repair_gtis=yes","region="+reg,"telid="+str(i)]
            subprocess.run(cmd)
            if n==0:
                hdu = fits.open(out_path+"evt_TM"+str(i)+".fits")[0]
                img=hdu.data
                x=max(img.shape[0],img.shape[1])
                x=round(x/10)*10
                cmd=["evtool", "eventfiles=@"+txt_file, "outfile= "+out_path+"evt_TM"+str(i)+".fits", "image=yes", "rebin=80","pattern=15", "size="+str(x), "gti=FLAREGTI", "center_position= auto","flag=0xc00f7f30","repair_gtis=yes","region="+reg,"telid="+str(i)]
                subprocess.run(cmd)
            Evtfile_converter(self.evtfiles[n],i,pnt,suffix="_conv")
            self.evtfiles_conv.append(self.evtfiles[n].replace(".fits","")+"_conv.fits")

            
            self.expmaps.append(out_path+"expmap_TM"+str(i)+".fits")
            cmd_exp=["expmap", "inputdatasets="+out_path+"evt_TM"+str(i)+".fits","templateimage="+out_path+"evt_TM"+str(i)+".fits","emin=0.2","emax=10.0","singlemaps="+out_path+"expmap_TM"+str(i)+".fits", "gtitype=FLAREGTI", "withdetmaps=yes","withsinglemaps=yes","withmergedmaps=no","withvignetting=no","withweights=no"]
            subprocess.run(cmd_exp)
    
    def find_data_products(self):
        # find already created data products
        self.evtfiles=[]
        self.expmaps=[]
        self.q_arf=[]
        if not self.pointed:
            self.psfmaps=[]
        self.evtfiles_conv=[]
        
        for i in self.tms:
            if path.isfile(self.out_path+"evt_TM"+str(i)+".fits"):
                self.evtfiles.append(self.out_path+"evt_TM"+str(i)+".fits")
            else:
                raise FileNotFoundError("File "+self.out_path+"evt_TM"+str(i)+".fits not found.")
                
            if path.isfile(self.out_path+"expmap_TM"+str(i)+".fits"):
                self.expmaps.append(self.out_path+"expmap_TM"+str(i)+".fits")
            else:
                raise FileNotFoundError("File "+self.out_path+"expmap_TM"+str(i)+".fits not found.")
            
            q_arf_map_path=glob.glob(self.out_path+"quantized_arf_exp_b*_"+str(i)+".fits.gz")
            if path.isfile(q_arf_map_path[0]):
                self.q_arf.append(q_arf_map_path[0])

            if not self.pointed:
                if path.isfile(self.out_path+"psfmap_"+str(i)+".fits.gz"):
                    self.psfmaps.append(self.out_path+"psfmap_"+str(i)+".fits.gz")
            
            if path.isfile(self.out_path+"evt_TM"+str(i)+"_conv.fits"):     
                self.evtfiles_conv.append(self.out_path+"evt_TM"+str(i)+"_conv.fits")
                
    
    def make_exposure(self,downsample=20):
        # make 
        self.q_arf=[]
        for n,i in enumerate(self.tms):
            make_arfs(input_evtfile=self.evtfiles[n], input_expmap=self.expmaps[n], output_path=self.out_path+"quantized_arf_exp_b"+str(downsample)+"_"+str(i)+".fits.gz", downsample=downsample, tm=i)
            make_arf_map(input_evtfile=self.evtfiles[n], input_expmap=self.expmaps[n], output_path=self.out_path+"quantized_arf_exp_b"+str(downsample)+"_"+str(i)+".fits.gz", downsample=downsample, tm=i)
            self.q_arf.append(self.out_path+"quantized_arf_exp_b"+str(downsample)+"_"+str(i)+".fits.gz")
            
        
    def make_psf(self,downsample=20,downsample_time=10):
        self.psfmaps=[]
        if self.pointed:
            return "Pointed data, PSFMap creation not necessary."
        else:
            for n,i in enumerate(self.tms):
                make_psfmap(evtfile=self.evtfiles[n], tm=i, outfile=self.out_path+"psfmap_"+str(i)+".fits.gz", downsample=downsample, downsample_time=downsample_time)
                self.psfmaps.append(self.out_path+"psfmap_"+str(i)+".fits.gz")

            
    def load_data(self,onoff= True, wid=0.6,binsize=8,e_reco_axis="original",e_true_axis="original",geom_file=None,geom_file_cutout_cen=None,binsize_irf=3,srctool_like=False):
        self.obs=[]
        self.datasets=[]
        self.binsize=binsize
        
        #energy_axis = MapAxis.from_bounds(0.2, 10, nbin=100, name="energy", unit="keV", interp="log")#1024

        #e_true_axis = MapAxis.from_bounds(
           # 0.2,
           # 10,
           # nbin=200,
           # unit="keV",
           # name="energy_true",
           # interp="log")

        # load data into a Dataset/DatasetOnOff
        for i,val in enumerate(self.tms):
            
            #Events & GTI
            events = EventList.read(self.evtfiles_conv[i])
            gti=GTI.read(self.evtfiles_conv[i])

            obs = Observation(events=events,gti=gti)
            obs.obs_id=val #needed for dataset creation
            obs._pointing=FixedPointingInfo(fixed_icrs=self.region.center.icrs)#v 2.0.1 needed for dataset creation
            
            #RMF
            rmf = EDispKernel.read(os.path.dirname(os.path.abspath(__file__))+"/data/RMF_ARF/120_RMF_00001_conv.fits") # change to RMF provided with script!
            if e_reco_axis == "original":
                energy_axis = MapAxis.from_edges(edges=rmf.axes["energy"].edges,unit="keV",name="energy",interp="log")
            if e_true_axis == "original":
                e_true_axis = MapAxis.from_edges(edges=rmf.axes["energy_true"].edges,unit="keV",name="energy_true",interp="log")
                
            if geom_file:
                hdulist=fits.open(geom_file)
                hdulist[0].header["RADECSYS"]="ICRS"
                geom = WcsGeom.from_header(hdulist[0].header).to_cube([energy_axis])

                if wid:
                    if type(geom_file_cutout_cen) == str:
                        geom_file_cutout_cen=Regions.parse(ds9_reg, format="ds9")[0]
                    geom = geom.cutout(position=geom_file_cutout_cen,width=wid*u.deg)
            else:
                geom = WcsGeom.create(skydir=(self.reg_coords[0],self.reg_coords[1]),axes=[energy_axis],width=wid * u.deg,binsz=binsize * u.arcsec,frame="icrs")
            
            dataset_empty = MapDatasetOnOff.create(geom=geom,energy_axis_true=e_true_axis,name="dataset"+str(val),binsz_irf=(binsize_irf*u.arcmin).to("deg").value)
            
            #PSF
            if not self.pointed:
                obs.psf = PSFMap.read(self.psfmaps[i],hdu="PSF",format="gadf")
            else:
                psf_kernel = PSF3D.read(os.path.dirname(os.path.abspath(__file__))+"/data/PSF/PSFrad_TM"+str(val)+".fits",hdu='PSF', format='gadf-dl3')
                obs.psf=utils.make_psf_map(psf=psf_kernel,pointing=obs.pointing.fixed_icrs,geom=dataset_empty.psf.psf_map.geom,)
            
            #read in an unrelated Aeff for TELESCOP keyword
            obs.aeff = RegionNDMap.read(os.path.dirname(os.path.abspath(__file__))+"/data/RMF_ARF/120_ARF_00001.fits", format="ogip-arf") #needed for TELESCOP keywords
            
            self.obs.append(obs)
            
            #run dataset maker
            if onoff:
                dataset_empty = MapDatasetOnOff.create(geom=geom,energy_axis_true=e_true_axis,rad_axis=obs.psf.psf_map.geom.axes["rad"],name="dataset"+str(val),binsz_irf=(binsize_irf*u.arcmin).to("deg").value)
            else:
                dataset_empty = MapDataset.create(geom=geom,energy_axis_true=e_true_axis,rad_axis=obs.psf.psf_map.geom.axes["rad"],name="dataset"+str(val),binsz_irf=(binsize_irf*u.arcmin).to("deg").value)

            maker = MapDatasetMaker(selection=['counts','psf'] ,background_interp_missing_data=False)
            dataset = maker.run(dataset_empty, obs)
            
            #load exposure
            
            #test_map=gammapy.maps.WcsNDMap.read(self.q_arf[i], hdu='PRIMARY', hdu_bands='PRIMARY_BANDS',map_type="wcs")
            #hdulist=fits.open(self.q_arf[i])
            #hdulist[0].header["RADECSYS"]="ICRS"
            dataset.exposure=gammapy.maps.WcsNDMap.read(self.q_arf[i], hdu='PRIMARY', hdu_bands='PRIMARY_BANDS',map_type="wcs").interp_to_geom(dataset_empty.exposure.geom,preserve_counts=False) #gammapy.maps.WcsNDMap.from_hdulist(hdulist, hdu='PRIMARY', hdu_bands='PRIMARY_BANDS').interp_to_geom(dataset_empty.exposure.geom,preserve_counts=False)


            #get original counts geom from eventfile
            pix_size_orig=round(WcsGeom.from_header(fits.open(self.evtfiles_conv[i])[0].header)._cdelt[0][0]*3600)
            pix_size=round(dataset.counts.geom._cdelt[0][0]*3600)
            
            if srctool_like:
                # apply half-pixel shift to counts, full pixel shift to exposure
                dataset.exposure=dataset.exposure.upsample(int(pix_size/(2*pix_size_orig)),preserve_counts=False)
                dataset.exposure.data=np.roll(dataset.exposure.data,(1,1),axis=(-2,-1))
                dataset.exposure=dataset.exposure.downsample(int(pix_size/(2*pix_size_orig)),preserve_counts=False)
                
                counts_shift = Map.from_geom(dataset.counts.geom.upsample(int(pix_size/(pix_size_orig))))
                counts_shift.fill_events(obs.events)
                counts_shift.data=np.roll(counts_shift.data,(1,1),axis=(-2,-1))
                dataset.counts=counts_shift.downsample(int(pix_size/(pix_size_orig)),preserve_counts=True)

            #else:
                # apply half-pixel shift to exposure
                #dataset.exposure=dataset.exposure.upsample(int(pix_size/(pix_size_orig)),preserve_counts=False)
                #dataset.exposure.data=np.roll(dataset.exposure.data,(1,1),axis=(-2,-1))
                #dataset.exposure=dataset.exposure.downsample(int(pix_size/(pix_size_orig)),preserve_counts=False)
                

            dataset.edisp=EDispKernelMap.from_edisp_kernel(rmf,geom=dataset_empty.edisp.edisp_map.geom)
            dataset.psf.exposure_map=dataset.exposure.interp_to_geom(dataset.psf.psf_map.geom.squash("rad"),preserve_counts=False)
            dataset.edisp.exposure_map=dataset.exposure.interp_to_geom(dataset_empty.edisp.edisp_map.geom.squash("energy"),preserve_counts=False)
            rmf=None

            if onoff:
                dataset.counts_off=dataset.counts.copy()
                dataset.counts_off.data=np.zeros_like(dataset.counts.data)
                dataset.acceptance=dataset.counts.copy()
                dataset.acceptance.data=np.ones_like(dataset.acceptance.data)
                dataset.acceptance_off=dataset.counts.copy()
                dataset.acceptance_off.data=np.zeros_like(dataset.acceptance_off.data)
            else:
                dataset.background=dataset.counts.copy()
                dataset.background.data=np.zeros_like(dataset.counts.data)

            self.datasets.append(dataset)
            continue
            

    def stack(self):
        stacked=self.datasets[0].copy()
        back_spec=self.backgr_spectra[0].copy()
        for i in range(1,len(self.tms)):
            stacked.stack(self.datasets[i])
            back_spec.stack(self.backgr_spectra[i])
        self.stacked=stacked
        self.backgr_stacked=back_spec

    def background_spectrum(self,backgr,region=None,to_mask_out=None,on_off=False,exp=True,ext_mask=None):
        if self.datasets==None:
            raise KeyError("Define a dataset via load_data()!")
            
        if backgr.datasets==None:
            raise KeyError("Define a background dataset via load_data()!")

        if (backgr.datasets[0].counts.geom.axes["energy"]!=self.datasets[0].counts.geom.axes["energy"]) or (backgr.datasets[0].exposure.geom.axes["energy_true"]!=self.datasets[0].exposure.geom.axes["energy_true"]):
            print("Incompatible energy axes!")

        if self.binsize!=backgr.binsize:
            print("Dataset and background binsize does not match!")

        self.backgr_spectra=[]
        self.backgr_ratios=[]
        self.backgr_masks=[]
        self.backgr_exp_ratios=[]
        for j,val1 in enumerate(self.tms):
            if ext_mask:
                # use pre-existing mask if specified
                bin_num=backgr.datasets[j].counts.data.shape[0]
                mask=backgr.datasets[j].mask_safe.copy()

                hdulist=fits.open(ext_mask)
                hdulist[0].header["RADECSYS"]="ICRS"
                mask_image = WcsNDMap.from_hdulist(hdulist)
                #mask_image=WcsNDMap.read(ext_mask)
                mask_image=mask_image.interp_to_geom(mask.geom.to_image(),method="nearest")
                mask.data=np.stack([mask_image.data]*bin_num,axis=0)
            else:
                # get point src regions in region
                regions=backgr.get_eRASS1_point_src()
                for i, val in enumerate(regions):
                    r=backgr.datasets[j].psf.containment_radius(0.9,energy_true=2*u.keV,position=val.center)
                    val.radius=r[0]
                if to_mask_out:
                    if type(to_mask_out[0])==str:
                        to_mask_out=list(Regions.parse(to_mask_out, format="ds9"))   
                    regions.extend(to_mask_out)
                geom=backgr.datasets[j].counts.geom
                mask_exclude= ~geom.region_mask(regions)
        
                #parse backgr region and create include mask:
        
                if type(region)== str:
                    region= Regions.parse(region, format="ds9")
        
                mask_include= geom.region_mask([region])
                #add masks with logical and:
        
                data_and=np.logical_and(mask_exclude.data,mask_include.data)
                mask=mask_exclude.copy()
                mask.data=data_and
            
            
            backgr.datasets[j].mask_safe=mask

            area_factor=np.sum(mask.data[0])*(backgr.datasets[j].counts.geom.cdelt/self.datasets[j].counts.geom.cdelt).decompose().value
            #print(data_and[0])
    
            #Extract spectrum:
            backgr_spectrum=to_spectrum_dataset_xray(self=backgr.datasets[j],on_region=region)
            hdulist=fits.open(backgr.expmaps[j])
            hdulist[0].header["RADECSYS"]="ICRS"
            backgr_expmap=WcsNDMap.from_hdulist(hdulist).interp_to_geom(backgr.datasets[j].counts.geom.to_image())
            #backgr_expmap=WcsNDMap.read(backgr.expmaps[j]).interp_to_geom(backgr.datasets[j].counts.geom.to_image())
            exp_ratio=np.mean(backgr_expmap.data, where=backgr.datasets[j].mask_image.data) # the where used to be given as weights!!!
            
            self.backgr_spectra.append(backgr_spectrum)
            self.backgr_ratios.append(area_factor)
            self.backgr_masks.append(mask.copy())
            self.backgr_exp_ratios.append(exp_ratio)

            backgr.datasets[j].mask_safe.data=np.ones_like(backgr.datasets[j].mask_safe.data,dtype=bool)

        if on_off:
            for i, val in enumerate(self.tms):
                self.datasets[i].counts_off=self.datasets[i].counts.copy()
                self.datasets[i].counts_off.data=np.tile(self.backgr_spectra[i].counts.data.flatten()[:,np.newaxis,np.newaxis],(1,self.datasets[i].counts.data.shape[1],self.datasets[i].counts.data.shape[2]))
                
                if exp:
                    expmap=WcsNDMap.read(self.expmaps[i]).interp_to_geom(self.datasets[i].counts.geom.to_image())
                    self.datasets[i].acceptance.data[:]=expmap.data
                    expmap=None
                    self.datasets[i].acceptance_off.data.fill(self.backgr_exp_ratios[i]*self.backgr_ratios[i])

                else:
                    acc_on=gammapy.estimators.utils.estimate_exposure_reco_energy(self.datasets[i])
                    self.datasets[i].acceptance=gammapy.maps.WcsNDMap(acc_on.geom, data=acc_on.data, unit="")
                    
                    acc_off=acc_on.copy()
                    acc_off_1d=gammapy.estimators.utils.estimate_exposure_reco_energy(self.backgr_spectra[i])
                    acc_off.data=np.tile(acc_off_1d.data.flatten()[:,np.newaxis,np.newaxis],(1,acc_on.data.shape[1],acc_on.data.shape[2]))*self.backgr_ratios[i]*((1*acc_off_1d.unit/acc_on.unit).decompose().value)
                    self.datasets[i].acceptance_off=gammapy.maps.WcsNDMap(acc_off.geom, data=acc_off.data, unit="")


    def write_to_files(self, filename=None,tms=None,stacked=False,overwrite=True):
        if tms:
            for i, val in enumerate(self.tms):
                if filename:
                    filename1=filename.replace(".fits","_TM"+str(val)+".fits")
                else:
                    filename1=self.out_path+"dataset_TM"+str(val)+".fits.gz"
                if val in tms:
                    self.datasets[i].write(filename1, overwrite=True, checksum=False)
                    self.backgr_spectra[i].write(filename.replace(".fits.gz","_backgr.fits.gz"), format="gadf", overwrite=True)

        if stacked:
            if not filename:
                filename=self.out_path+"dataset_stacked.fits.gz"
            self.stacked.write(filename, overwrite=True, checksum=False)
            self.backgr_stacked.write(filename.replace(".fits.gz","_backgr.fits.gz"), format="gadf", overwrite=True)


def make_energy_mask(dataset,energy_bounds,reset=False):
    if reset:
        mask=np.ones_like(dataset.mask.data, dtype=bool)
    else:
        mask = dataset.mask.data
    edges = dataset.counts.geom.axes["energy"].edges
    el = energy_bounds[0]
    eu = energy_bounds[1]
    
    for i in range(mask[:,0,0].shape[0]):
        if not ((edges[i] >= el and edges[i] <= eu) and (edges[i+1] >=el and edges[i+1]<=eu)):
            mask[i]=np.zeros_like(mask[i], dtype=bool)

    dataset.mask_safe.data = mask

def to_spectrum_dataset_xray(self, on_region, name=None):#, containment_correction=False, name=None):
    #No longer works since 2.0:
    #dataset = MapDataset.to_spectrum_dataset(self,
    #    on_region=on_region,
    #    containment_correction=containment_correction,
    #    name=name,
    #)

    name = make_name(name)
    kwargs = {"gti": self.gti, "name": name, "meta_table": self.meta_table}

    if self.mask_safe:
        kwargs["mask_safe"] = self.mask_safe.to_region_nd_map(on_region, func=np.any)

    if self.mask_fit:
        kwargs["mask_fit"] = self.mask_fit.to_region_nd_map(on_region, func=np.any)

    if self.counts:
        kwargs["counts"] = self.counts.to_region_nd_map(
            on_region, np.sum, weights=self.mask_safe
        )

    if self.stat_type == "cash" and self.background:
        kwargs["background"] = self.npred_background().to_region_nd_map(
            on_region, func=np.sum, weights=self.mask_safe
        )

    if self.exposure:
        kwargs["exposure"] = None

    if self.psf:
        kwargs["psf"] = self.psf.to_region_nd_map(on_region)

    if self.edisp is not None:
        kwargs["edisp"] = self.edisp.to_region_nd_map(on_region)

    dataset = SpectrumDataset(**kwargs)
    # correct all the mean:

    kwargs = {"name": name}
    area_factor = self.mask_safe.get_spectrum(on_region, np.sum, weights=self.mask_safe)

    mask_e_true = self.exposure.copy()
    
    #for i in range(self.mask_safe.data.shape[0]):
    #    if np.any(self.mask_safe.data[i]):
    #        mask_slice = self.mask_safe.data[i]
    mask_e_true.data[:] = self.mask_image.data#mask_slice

    area_factor_e_true = mask_e_true.get_spectrum(on_region, func=np.sum, weights=mask_e_true)#*(self.exposure.sum_over_axes()).data)#(self.counts.sum_over_axes().data))

    #self.exposure.data
    exposure = self.exposure.to_region_nd_map(on_region, func=np.sum, weights=mask_e_true)#*(self.exposure.sum_over_axes()).data)#(self.counts.sum_over_axes().data))#*(self.exposure.data/(np.sum(self.exposure.data*mask_e_true.data, axis=(-2,-1)))[:,np.newaxis,np.newaxis]))
    #exposure.data = np.nan_to_num(np.divide(exposure.data,np.multiply(area_factor_e_true.data,np.sum(self.exposure.data, axis=(-2,-1))[:,np.newaxis,np.newaxis])))
    exposure.data = np.nan_to_num(np.divide(exposure.data,area_factor_e_true.data))
    dataset.exposure = exposure
    try:
        counts_off = self.counts_off.get_spectrum(
            on_region, np.sum, weights=self.mask_safe
        )
        counts_off.data = np.nan_to_num(np.divide(counts_off.data,area_factor.data))
        kwargs["counts_off"] = counts_off

        acceptance = self.acceptance.get_spectrum(
            on_region, np.sum, weights=self.mask_safe
        )
        acceptance.data = np.nan_to_num(np.divide(acceptance.data,area_factor.data))
        kwargs["acceptance"] = acceptance
        
        acceptance_off = self.acceptance_off.get_spectrum(on_region, np.sum, weights=self.mask_safe)#kwargs["acceptance"] * kwargs["counts_off"] / norm
        area_factor = self.mask_safe.get_spectrum(on_region, np.sum, weights=self.mask_safe)
        #print(area_factor.data)
        #acceptance_off.data = np.nan_to_num(np.divide(acceptance_off.data,(area_factor.data)))
        acceptance_off.data = np.nan_to_num(np.divide(acceptance_off.data,(area_factor.data*area_factor.data)))
        #np.nan_to_num(acceptance_off.data, copy=False)
        kwargs["acceptance_off"] = acceptance_off

        return SpectrumDatasetOnOff.from_spectrum_dataset(dataset=dataset, **kwargs)
    except:
        return dataset

def wstat_rebin(self,n=5):
        z=[0]
        x=0

        backgr_spec= np.sum(self.counts_off.data,axis=(1,2), where=self.mask_safe.data)/np.sum(self.mask_safe.data,axis=(1,2))
        #np.sum(self.dataset[j].counts_off.data,axis=(1,2))/np.sum(self.dataset[j].mask_safe.data,axis=(1,2))
        #backgr_spec=np.sum([self.backgr_spectra[i].counts.data.flatten() for i in range(len(self.backgr_spectra)) ],axis=0)
        
        for i,val in enumerate(backgr_spec):
            x=x+val
            #if self.stacked.counts_off.data[i,0,0]!=0.0:
                #z.append(i)
            if x>=n:
                z.append(i+1)
                x=0

        edges_new=[]
        edges_old=self.counts.geom.axes["energy"].edges#.value
        for i, val in enumerate(z):
            edges_new.append(edges_old.value[int(val)])
        edges_new.append(edges_old.value[-1])

        energy_axis_new = MapAxis.from_edges(edges_new,name="energy", unit=edges_old.unit, interp="lin")

        return self.resample_energy_axis(energy_axis_new)

def rebin_energy_bounds(self,low_E=0.2*u.keV,high_E=10.0*u.keV):
        edges_new=[]
        edges_old=self.counts.geom.axes["energy"].edges#.value
        for i, val in enumerate(edges_old):
            if val>=low_E and val<=high_E:
                edges_new.append(val.value)
        edges_new.insert(0,edges_old[0].value)
        edges_new.append(edges_old[-1].value)

        energy_axis_new = MapAxis.from_edges(edges_new,name="energy", unit=edges_old.unit, interp="lin")
        return self.resample_energy_axis(energy_axis_new)

def to_lin(self):
    energy_axis_lin = MapAxis.from_edges(edges=self.counts.geom.axes["energy"].edges,unit="keV",name="energy",interp="lin")
    e_true_axis_lin = MapAxis.from_edges(edges=self.exposure.geom.axes["energy_true"].edges,unit="keV",name="energy_true",interp="lin")
    geom_lin=self.counts.geom.copy().replace_axis(energy_axis_lin)

    if self.psf and (type(self) not in [gammapy.datasets.spectrum.SpectrumDataset,gammapy.datasets.spectrum.SpectrumDatasetOnOff]):
        rad_axis=self.psf.psf_map.geom.copy().axes["rad"]
        binsz_irf=self.edisp.edisp_map.geom.pixel_scales[0]
        dataset_lin = type(self).create(geom=geom_lin,energy_axis_true=e_true_axis_lin,rad_axis=rad_axis,binsz_irf=binsz_irf)
    else:
        rad_axis=None
        binsz_irf=None
        dataset_lin = type(self).create(geom=geom_lin,energy_axis_true=e_true_axis_lin)
    
    #dataset_lin = type(self).create(geom=geom_lin,energy_axis_true=e_true_axis_lin,rad_axis=rad_axis,binsz_irf=binsz_irf)

    dataset_lin.counts.data=self.counts.data
    if type(self) in [gammapy.datasets.spectrum.MapDatasetOnOff,gammapy.datasets.spectrum.SpectrumDatasetOnOff]:
        dataset_lin.counts_off.data=self.counts_off.data
        dataset_lin.acceptance.data=self.acceptance.data
        dataset_lin.acceptance_off.data=self.acceptance_off.data
    elif self.background:
        dataset_lin.background.data=self.background.data
    else:
        pass
    try:
        dataset_lin.gti=self.gti.copy()
    except:
        pass
    dataset_lin.exposure.data=self.exposure.data*((self.exposure.unit/dataset_lin.exposure.unit).decompose().scale)
    dataset_lin.edisp.edisp_map.data=self.edisp.edisp_map.data
    dataset_lin.edisp.exposure_map.data=self.edisp.exposure_map.data
    if self.psf and (type(self) not in [gammapy.datasets.spectrum.SpectrumDataset,gammapy.datasets.spectrum.SpectrumDatasetOnOff]):
        dataset_lin.psf.psf_map.data=self.psf.psf_map.data
        dataset_lin.psf.exposure_map.data=self.psf.exposure_map.data
    else:
        dataset_lin.psf=None
    dataset_lin.mask_safe.data=self.mask_safe.data
    
    return dataset_lin

def get_eRASS1_point_src(self,coord,radius,catalog_path):
    erass = Table(fits.open(catalog_path+"eRASS1_Main.v1.1.fits")[1].data)
    catalog=SkyCoord(ra=erass["RA"].data*u.deg,dec=erass["DEC"].data*u.deg)
    c = SkyCoord(ra=[coord.ra.to(u.deg).value]*u.degree, dec=[coord.dec.to(u.deg).value]*u.degree)
    idx1, idx2, sep, dist=catalog.search_around_sky(c,radius*u.deg)
    return Regions([CircleSkyRegion(center=i,radius=30*u.arcsec) for i in catalog[idx2]])

def export_pyxspec_model(outfile,name=""):
    xspec.Plot.xAxis = "keV"
    xspec.Plot("model "+name)
    energy=xspec.Plot.x()
    #xspec.Plot.xAxis = "A"
    #xspec.Plot("eemodel "+name)
    model=xspec.Plot.model()

    f= open(outfile,"w")
    #f.write("# energy (keV), flux (ergs cm$^{-2}$ s$^{-1}$) \n")
    f.write("# energy (keV), flux (1/ (keV cm2 s) ) \n")
    for i in range(len(model)):
        f.write(str(energy[i])+" , "+str(model[i])+"\n")
    f.close()

def import_erosita_background_model(filename, fwc=False, spatial=True, **kwargs):
    if fwc:
        y_unit="1/keV"
    else:
        y_unit="1/ ( keV cm2 s1)"
    
    vals = np.loadtxt(filename,skiprows=1,delimiter=",")
    energy = u.Quantity(vals[:, 0], "keV")
    values = u.Quantity(vals[:, 1], y_unit)

    kwargs.setdefault("interp_kwargs", {"fill_value": None})

    if spatial:
        spatial_model = ConstantSpatialModel()
    else:
        spatial_model=None
#    spectral_model = (
#        TemplateSpectralModel(energy=energy, values=values, **kwargs)
#        * PowerLawNormSpectralModel()
#    )
    template = TemplateSpectralModel(energy=energy, values=values, **kwargs)

    if fwc:
        name="fwc"
        apply_irf={"psf": False, "exposure": False, "edisp": True}
        spectral_model = (
        #TemplateSpectralModel(energy=energy, values=values, **kwargs)
        template
        * PowerLawNormSpectralModel()
    )
    else:
        name="bkg"
        apply_irf={"psf": False, "exposure": True, "edisp": True}
        spectral_model = (
        #TemplateSpectralModel(energy=energy, values=values, **kwargs)
        template
        * PowerLawNormSpectralModel()
    )
    
    return SkyModel(
        spatial_model=spatial_model,
        spectral_model=spectral_model,
        name=name,
        apply_irf=apply_irf,
    )

def make_exclusion_mask(self,include=None, exclude=None, point_src_exclusion=True, catalog_path=None):
    # make exclusion mask
    regions = []
    if point_src_exclusion:
        if not catalog_path:
            raise ValueError('Path to eRASS1 catalog needed, if point source exclusion should be performed.')
        else:
            regions=get_eRASS1_point_src(self,self.counts.geom.center_skydir,np.max(self.counts.geom.width).value*np.sqrt(2),catalog_path=catalog_path)
    if exclude:
        regions.extend(exclude)

    geom = self._geom
    mask_exclude= ~geom.region_mask(regions)

    if include:
    # make inclusion mask
        mask_include= geom.region_mask(include)
        data_and=np.logical_and(mask_exclude.data,mask_include.data)
        mask_exclude.data=data_and

    return mask_exclude

def fit_erosita_background(dataset,tm="stacked"):
    # fit FWC model in 6 to 9 keV range
    fwc=import_erosita_background_model(os.path.dirname(os.path.abspath(__file__))+"/data/FWC_models/backgr_FWC_model_DR1_"+str(tm)+".txt",fwc=True)
    dataset.models=[fwc]
    make_energy_mask(dataset,[6.0,9.0]*u.keV)
    fit = Fit()
    result = fit.run(datasets=[dataset])
    fwc.parameters["norm"].frozen=True

    # define model for diffuse background (to be changed to YAML file)
    area_arcmin=dataset.counts.geom.region.to_pixel(dataset.counts.geom.wcs).area*(dataset.counts.geom.binsz_wcs[0]**2).to(u.arcmin**2)
    from sherpa.astro.xspec import XSTBabs,XSapec,XSpowerlaw,XSParameter, XSgaussian, XSexpfac, XSParameter, XSbkn2pow, XSvnei, XSconstant
    from sherpa.models import parameter
    from sherpa.astro import xspec

    xspec.set_xsabund('wilm')
    # area constant
    const_sherpa = XSconstant()
    const_sherpa.factor = area_arcmin.value
    const_sherpa._use_caching=False
    const_sherpa.freeze()
    # TBabs absorption
    absorption_sherpa = XSTBabs()
    absorption_sherpa.nH = 1.36
    absorption_sherpa._use_caching=False
    absorption_sherpa.freeze()
    # LHB
    lhb_sherpa = XSapec()
    lhb_sherpa._use_caching=False
    lhb_sherpa.pars[0].name="kT_lhb"
    lhb_sherpa.pars[0].val=0.1
    lhb_sherpa.pars[0].min=0.09
    lhb_sherpa.pars[0].max=0.12
    lhb_sherpa.pars[0].freeze()
    lhb_sherpa.pars[1].freeze()
    lhb_sherpa.pars[2].freeze()
    lhb_sherpa.pars[-1].name="norm_lhb"
    lhb_sherpa.pars[-1].val=1e-6
    lhb_sherpa.pars[-1].min=1e-8
    lhb_sherpa.pars[-1].max=1e-4
    #CGM
    cgm_sherpa = XSapec()
    cgm_sherpa._use_caching=False
    cgm_sherpa.pars[0].name="kT_cgm"
    cgm_sherpa.pars[0].val=0.16
    cgm_sherpa.pars[0].min=0.15
    cgm_sherpa.pars[0].max=0.25
    cgm_sherpa.pars[1].name="Abundanc_cgm"
    cgm_sherpa.pars[1].val=0.08
    cgm_sherpa.pars[1].min=0.04
    cgm_sherpa.pars[1].max=2.0
    cgm_sherpa.pars[1].freeze()
    cgm_sherpa.pars[-1].name="norm_cgm"
    cgm_sherpa.pars[-1].val=4e-7
    cgm_sherpa.pars[-1].min=4e-9
    cgm_sherpa.pars[-1].max=4e-2
    # Cor
    cor_sherpa = XSapec()
    cor_sherpa._use_caching=False
    cor_sherpa.pars[0].name="kT_cor"
    cor_sherpa.pars[0].val=0.7
    cor_sherpa.pars[0].freeze()
    cor_sherpa.pars[1].name="Abundanc_cor"
    cor_sherpa.pars[1].val=1
    cor_sherpa.pars[1].freeze()
    cor_sherpa.pars[-1].name="norm_cor"
    cor_sherpa.pars[-1].val=1e-6
    cor_sherpa.pars[-1].min=1e-8
    cor_sherpa.pars[-1].max=1e-4
    # CXB
    cxb_sherpa = XSpowerlaw()
    cxb_sherpa._use_caching=False
    cxb_sherpa.pars[0].val=1.46
    cxb_sherpa.pars[0].freeze()
    cxb_sherpa.pars[-1].name="norm_cxb"
    cxb_sherpa.pars[-1].val=8.88e-7
    cxb_sherpa.pars[-1].freeze()
    
    spectral_model_sherpa=const_sherpa*(lhb_sherpa + absorption_sherpa*(cgm_sherpa+cor_sherpa +cxb_sherpa))
    spectral_model = SherpaSpectralModel(spectral_model_sherpa)
    
    manual_bkg = SkyModel(spectral_model=spectral_model,name="manual_bkg")

    # set both bkg components on dataset
    dataset.models=[manual_bkg,fwc]
    make_energy_mask(dataset,[0.3,9.0]*u.keV,reset=True)

    fwc.parameters["norm"].min=0.1
    
    fit = Fit(backend="sherpa")
    result = fit.run(datasets=[dataset])

    manual_bkg.parameters["Abundanc_cgm"].frozen=False
    manual_bkg.parameters["Abundanc_cgm"].max=5.0
    cgm_sherpa.pars[1].max=5.0

    fit = Fit(backend="sherpa")
    result = fit.run(datasets=[dataset])

    manual_bkg.parameters["norm_cxb"].frozen=False

    fit = Fit(backend="sherpa")
    result = fit.run(datasets=[dataset])

    fwc.parameters["norm"].frozen=False

    fit = Fit(backend="sherpa")
    result = fit.run(datasets=[dataset])

    print(result)



def erosita_background_to_template(model, energy_range=[0.1*u.keV,10.5*u.keV]):
    energies = np.linspace(energy_range[0],energy_range[1],2000)
    values = model.spectral_model(energies)
 
    t = TemplateSpectralModel(energy = energies, values = values)
    bkg = SkyModel(spectral_model=t*PowerLawNormSpectralModel(),spatial_model=ConstantSpatialModel(),name="bkg")
    return bkg












