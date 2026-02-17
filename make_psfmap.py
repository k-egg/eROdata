#import gammapy
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from regions import SkyRegion
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs import utils
#from regions import SkyRegion

# %matplotlib inline
import matplotlib.pyplot as plt
#from IPython.display import display
#from gammapy.data import EventList, Observation, GTI
#from gammapy.datasets import Datasets, MapDataset
#from gammapy.irf import EDispKernelMap, PSFMap,EDispKernel,EffectiveAreaTable2D
#from gammapy.maps import Map, MapAxis, WcsGeom,RegionNDMap,RegionGeom,MapAxes
#from gammapy.makers import MapDatasetMaker
#from gammapy.modeling import Fit
#from gammapy.modeling.models import (
#    Models,
#    PointSpatialModel,
#    PowerLawNormSpectralModel,
#    PowerLawSpectralModel,
#    SkyModel,
#    TemplateSpatialModel,
#    create_fermi_isotropic_diffuse_model,
#)
import numpy as np
from os import path
import os, os.path, time, subprocess
from astropy.wcs import WCS
import argparse

def make_psfmap(evtfile,tm, outfile,downsample=1):
    
    #read table of pointing positions from event file
    #(integrating readPointing into this function for convenience):
    
    hdul=fits.open(evtfile)
    P_table = np.array(hdul["CORRATT"+str(tm)].data.tolist())[:,1:3]
    
    img_dim_x = hdul[0].shape[0]
    img_dim_y = hdul[0].shape[1]
    img_coord = WCS(hdul[0].header)
    
    
    #create array with skycoord for every pixel
    coord_arr=np.zeros((img_dim_y,img_dim_x,2))
    for i in range(img_dim_y):
        for j in range(img_dim_x):
            #skycoor=utils.pixel_to_skycoord(j,i,wcs)
            coord_arr[i,j,0]=j
            coord_arr[i,j,1]=i

    coord_arr=np.array(img_coord.pixel_to_world_values(coord_arr.reshape((img_dim_y**2,2)))).reshape((img_dim_y,img_dim_x,2))
    
    #prepare output array:
    arr_a=np.zeros((img_dim_y,img_dim_x,6))
    
    for i in range(len(P_table)):
        print(i)
        p_coord=(P_table[i,0],P_table[i,1])#SkyCoord(P_table[i,0],P_table[i,1],unit="deg",frame="fk5")
        for j in range(img_dim_y):#len(coord_arr)):
            for k in range(img_dim_x):#len(coord_arr)):
                dist=np.rad2deg(np.arccos( np.sin(np.deg2rad(coord_arr[j,k,1]))*np.sin(np.deg2rad(p_coord[1])) \
                          + np.cos(np.deg2rad(coord_arr[j,k,1]))*np.cos(np.deg2rad(p_coord[1])) \
                          *np.cos(np.deg2rad(coord_arr[j,k,0]- p_coord[0]))))
                #np.sqrt(((coord_arr[j,k,0]-p_coord[0])**2)+ ((coord_arr[j,k,1]-p_coord[1])**2))#p_coord.separation(coord_arr[j,k])
                
                if dist < 3*(1/60):
                    arr_a[j,k,0]+=1
                elif dist < 9*(1/60):
                    arr_a[j,k,1]+=1
                elif dist < 15*(1/60):
                    arr_a[j,k,2]+=1
                elif dist < 21*(1/60):
                    arr_a[j,k,3]+=1
                elif dist < 27*(1/60):
                    arr_a[j,k,4]+=1
                elif dist < 30*(1/60):
                    arr_a[j,k,5]+=1
    
    #open PSF file
    hdul2=fits.open("/home/wecapstor1/caph/mppi147h/pwn_analysis/eROSITA_calib/PSF/PSFrad_TM"+str(tm)+"_100_new.fits")
    x=np.array(hdul2[1].data.tolist()[0][6])
    b=hdul2[1].data.tolist()[0]
    rad_lo=np.array(b[4]*7)
    rad_hi=np.array(b[5]*7)
    rad=np.array([(rad_hi[i]+rad_lo[i])/2 for i in range(len(rad_lo))])
    
    e_lo=np.repeat(b[0],len(b[4]))
    e_hi=np.repeat(b[1],len(b[4]))
    e_true=np.array([(e_hi[i]+e_lo[i])/2 for i in range(len(e_lo))])
    
    
    #define output array
    psf_map_arr=np.zeros((7,100,int(img_dim_y/downsample),int(img_dim_x/downsample)))
    
    #average PSF curves in accordance to the distributions found in arr_a
    #divide arr_a through vignetting factors (outer FoV has less exposure and larger factor)
    
    v_f=[1.01652741, 1.16308784, 1.45869291, 1.84955096, 2.36298895, 2.83500886]
    arr_a=np.divide(arr_a,v_f)
    
    for i in range(0,img_dim_y,downsample):
        for j in range(0,img_dim_x,downsample):
            for e in range(7):
                try:
                    avg_psf=np.average(x[:,:,e],axis=1, weights=np.average(arr_a[i:i+downsample,j:j+downsample],axis=(0,1)))
                    psf_map_arr[e,:,int(i/downsample),int(j/downsample)]=avg_psf
                except:
                    pass
    #creating the output file:
    #empty Primary image:
    
    hdu_empty = fits.PrimaryHDU()
    
    #PSF Map:
    
    hdu_psfmap = fits.ImageHDU(psf_map_arr, header=hdul[0].header)
    hdu_psfmap.header["RADECSYS"]="ICRS"
    hdu_psfmap.header["BUNIT"] = "arcsec-2"
    hdu_psfmap.name='PSF'
    hdu_psfmap.header['NAXIS']=4
    #hdu_psfmap.header['NAXIS1']=890
    #hdu_psfmap.header['NAXIS2']=890
    hdu_psfmap.header['NAXIS3']=100
    hdu_psfmap.header['NAXIS4']=7
    hdu_psfmap.header['BANDSHDU']='PSF_BANDS'
    #hdu_psfmap.header['WCSSHAPE']='(890,890,100,7)'
    hdu_psfmap.header['EXTNAME'] = 'PSF'
    hdu_psfmap.header['INTERP2'] = 'log     '
    hdu_psfmap.header['AXCOLS2'] = 'ENERGY_TRUE_MIN,ENERGY_TRUE_MAX'
    hdu_psfmap.header['INTERP1'] = 'lin     '
    hdu_psfmap.header['AXCOLS1'] = 'RAD_MIN,RAD_MAX'
    
    # create downsampled WCS in header:
    hdu_psfmap.header["CDELT1"]=hdu_psfmap.header["CDELT1"]*downsample
    hdu_psfmap.header["CDELT2"]=hdu_psfmap.header["CDELT2"]*downsample
    hdu_psfmap.header["CRPIX1"]=((hdu_psfmap.header["CRPIX1"]-0.5)/downsample)+0.5
    hdu_psfmap.header["CRPIX2"]=((hdu_psfmap.header["CRPIX2"]-0.5)/downsample)+0.5
    hdu_psfmap.header["NAXIS1"]=hdu_psfmap.header["NAXIS1"]/downsample
    hdu_psfmap.header["NAXIS2"]=hdu_psfmap.header["NAXIS2"]/downsample
    hdu_psfmap.header['WCSSHAPE']='('+str(hdu_psfmap.header["NAXIS1"])+','+str(hdu_psfmap.header["NAXIS2"])+',100,7)'
    
    #PSF Bands extension:
    col1=fits.Column(name='CHANNEL', format='K',array=range(len(rad_lo)))
    col2=fits.Column(name='RAD', format='D',array=rad,unit='arcsec')
    col3=fits.Column(name='RAD_MIN', format='D',array=rad_lo,unit='arcsec')
    col4=fits.Column(name='RAD_MAX', format='D',array=rad_hi,unit='arcsec')
    col5=fits.Column(name='ENERGY_TRUE', format='D',array=e_true,unit='keV')
    col6=fits.Column(name='ENERGY_TRUE_MIN', format='D',array=e_lo,unit='keV')
    col7=fits.Column(name='ENERGY_TRUE_MAX', format='D',array=e_hi,unit='keV')
    hdu2 = fits.BinTableHDU.from_columns([col1, col2, col3,col4,col5,col6,col7])
    hdu2.name = 'PSF_BANDS'
    hdu2.header['EXTNAME']='PSF_BANDS'
    hdu2.header['AXCOLS2'] = 'ENERGY_TRUE_MIN,ENERGY_TRUE_MAX'
    hdu2.header['INTERP1'] = 'lin     '
    hdu2.header['AXCOLS1'] = 'RAD_MIN,RAD_MAX'
    
    
    hdul_psfmap = fits.HDUList([hdu_empty,hdu_psfmap,hdu2])
    hdul_psfmap.writeto(outfile,overwrite=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    #input_evtfile,input_expmap,input_mask=None,lower,upper
    parser.add_argument("--evtfile", required=True, type=str, help="Input calibrated eROSITA eventfile")
    parser.add_argument("--tm", required=True, type=int, help="TM to analyze. Choose 1, 2, 3, 4, or 6.")
    parser.add_argument("--outfile", required=True, type=str, help="Output filepath for PSFMap")

    parser.add_argument("--downsample", required=False, type=int, help="Factor for downsampling of PSFMap")

    
    args = parser.parse_args()

    make_psfmap(args.evtfile,args.tm,args.outfile,args.downsample)    

