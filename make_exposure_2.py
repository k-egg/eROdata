# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import utils

import matplotlib.pyplot as plt

import numpy as np
from os import path
import os, os.path, time, subprocess
from astropy.wcs import WCS
import argparse
import glob

def make_arf_map(input_evtfile,input_expmap,arf_path=None,output_path=None,downsample=2,tm=1,clobber=False):
    # make output path:
    if output_path==None:
        filepath_long=os.path.abspath(__file__)
        filepath=filepath_long[0:filepath_long.rfind("/")+1]
    else:
        filepath=output_path
        if filepath[-1]!="/":
            filepath=filepath+"/"
    
    arf_expmap_file = filepath+"quantized_arf_exp_b"+str(downsample)+"_"+str(tm)+".fits.gz"

    # get WCS information and exposure time data
    hdu = fits.open(input_evtfile)[0]
    img = hdu.data
    wcs = WCS(hdu.header)
    y_len=img.shape[0]
    x_len=img.shape[1]
    
    hdu_exp = fits.open(input_expmap)[0]
    img_exp = hdu_exp.data

    # create downsampled WCS and output arrays:
    x_len_out=int(x_len/downsample)
    y_len_out=int(y_len/downsample)

    # try writing into existing file
    if clobber==False:
        try:
            out_arr_exp= fits.open(arf_expmap_file)[0].data

        except:
            out_arr_exp=np.zeros((1024,y_len_out,x_len_out))
    else:
        out_arr_exp=np.zeros((1024,y_len_out,x_len_out))
    
    
    #find files in ARF folder:
    files=glob.glob(arf_path+"*"+str(tm)+"20_ARF*"+str(tm)+".fits")
    
    for i in files:
        # detemine correct pixel and read data
        d=i[i.rfind("ARF_00001_")+10:i.rfind("_")]
        x=int(d[0:d.find("_")])
        y=int(d[d.find("_")+1:])

        hdu_arf = fits.open(i)[1]
        data_arf = np.array(hdu_arf.data.tolist())[:,2]

        # delete ARF file:
        os.remove(i)

        # put ARF data into correct output map pixel
        out_arr_exp[:,int(y/downsample),int(x/downsample)]=data_arf*np.average(img_exp[y:y+downsample,x:x+downsample])



    energy_low = np.array(hdu_arf.data.tolist())[:,0]
    energy_high = np.array(hdu_arf.data.tolist())[:,1]
    energy = (energy_low + energy_high)/2
    channel = np.array(range(1,len(energy)+1))
    header=hdu_arf.header



    #create output ARF FITS file:
    #file format is designed to be read in as a WcsNDMap

    col1=fits.Column(name='CHANNEL', format='I',array=channel)
    col2=fits.Column(name='ENERGY_TRUE', format='D',array=energy,unit='keV')
    col3=fits.Column(name='ENERGY_TRUE_MIN', format='D',array=energy_low,unit='keV')
    col4=fits.Column(name='ENERGY_TRUE_MAX', format='D',array=energy_high,unit='keV')
    hdu2 = fits.BinTableHDU.from_columns([col1, col2, col3,col4])
    hdu2.name = 'PRIMARY_BANDS'
    hdu2.header['AXCOLS1']='ENERGY_TRUE_MIN,ENERGY_TRUE_MAX'
    

    #change NAXIS in image header
    hdu.header['NAXIS']=3
    hdu.header['NAXIS3']=1024
    #add BANDSHDU for energy axis
    hdu.header['BANDSHDU']='PRIMARY_BANDS'
    #Define WCS
    hdu.header["RADECSYS"]="ICRS"
    
    # create downsampled WCS in header:
    hdu.header["CDELT1"]=hdu.header["CDELT1"]*downsample
    hdu.header["CDELT2"]=hdu.header["CDELT2"]*downsample
    hdu.header["CRPIX1"]=((hdu.header["CRPIX1"]-0.5)/downsample)+0.5
    hdu.header["CRPIX2"]=((hdu.header["CRPIX2"]-0.5)/downsample)+0.5
    hdu.header["NAXIS1"]=hdu.header["NAXIS1"]/downsample
    hdu.header["NAXIS2"]=hdu.header["NAXIS2"]/downsample
    
    #Remove physical coordinate system:
    hdu.header.remove("WCSNAMEP")
    hdu.header.remove("WCSTY1P")
    hdu.header.remove("WCSTY2P")
    hdu.header.remove("CTYPE1P")
    hdu.header.remove("CTYPE2P")
    hdu.header.remove("CRPIX1P")
    hdu.header.remove("CRVAL1P")
    hdu.header.remove("CDELT1P")
    hdu.header.remove("CRPIX2P")
    hdu.header.remove("CRVAL2P")
    hdu.header.remove("CDELT2P")


    #Prepare output file:

    #Exposure 3D
    hdu_arffile_exp = fits.PrimaryHDU(out_arr_exp, header=hdu.header)
    hdu_arffile_exp.header["BUNIT"] = "cm^2*s"
    hdu_arffile_exp.name='PRIMARY'
    hdul_arffile_exp = fits.HDUList([hdu_arffile_exp,hdu2])

    hdul_arffile_exp.writeto(arf_expmap_file,overwrite=True)

    
if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--input_evtfile", required=True, type=str, help="Input calibrated eROSITA eventfile")
    parser.add_argument("--input_expmap", required=True, type=str, help="Input eROSITA expmap")
    parser.add_argument("--arf_path", required=False, type=str, help="Input ARF path")
    parser.add_argument("--output_path", required=False, type=str, help="Filepath for output files")
    parser.add_argument("--downsample", required=False, type=int, help="Bin parameter. 2 means the map is binned in steps of 2x2=4 pixels.")
    parser.add_argument("--tm", required=False, type=int, help="TM to analyze. Choose 1, 2, 3, 4, or 6.")
    parser.add_argument("--clobber", required=False, type=bool, help="TM to analyze. Choose 1, 2, 3, 4, or 6.")

    args = parser.parse_args()

    make_arf_map(args.input_evtfile,args.input_expmap,args.arf_path,args.output_path,args.downsample,args.tm,args.clobber)    
