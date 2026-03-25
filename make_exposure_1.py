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

def make_arfs(input_evtfile,output_path=None,x_lo=None, x_hi=None, y=None,downsample=2,tm=1):

    # make output path:
    if not output_path:
        filepath_long=os.path.abspath(__file__)
        filepath=filepath_long[0:filepath_long.rfind("/")+1]
    else:
        filepath=output_path
        if filepath[-1]!="/":
            filepath=filepath+"/"

    # open evtfile and get image size parameters
    hdu = fits.open(input_evtfile)[0]
    img = hdu.data
    wcs = WCS(hdu.header)
    img_mask = np.zeros(img.shape)
    arf_mask = np.zeros(img.shape)
    y_len=img.shape[0]
    x_len=img.shape[1]
    pix_num=x_len*y_len
    pixel_rad=abs(hdu.header["CDELT1"])*downsample*np.sqrt(2) # calculate diagonal of exposure map pixel

    # adjust parameters, so all arf files are created if the function is not called as part of snakemake:
    if not x_lo:
        x_lo=0
    if not x_hi:
        x_hi=x_len-downsample
    if y is None:
        y=range(0,y_len - downsample,downsample)
    else:
        y=[y]

  
    # create downsampled WCS:
    x_len_out=int(x_len/downsample)
    y_len_out=int(y_len/downsample)
        
    # calculate array with coords for every pixel:
    coord_arr=np.zeros((y_len,x_len,2))
    for i in range(y_len):
        for j in range(x_len):
            coord_arr[i,j,0]=j
            coord_arr[i,j,1]=i

    coord_arr=np.array(wcs.pixel_to_world_values(coord_arr.reshape((y_len**2,2)))).reshape((y_len,x_len,2))
    
    mask=np.ones(img.shape)

    for i in y:
        for j in range(x_lo,x_hi,downsample):
            
            #produce one-pixel-mask:
            img_mask = np.zeros(img.shape)

            img_mask[i:i+downsample,j:j+downsample]=mask[i:i+downsample,j:j+downsample]

            hdu_out = fits.PrimaryHDU(img_mask, header=hdu.header)
            hdul = fits.HDUList([hdu_out])
            hdul.writeto(filepath+"mask_temporary_"+str(j)+"_"+str(i)+"_"+str(tm)+".fits.gz",overwrite=True)

            # extract one-pixel-ARF:
            cmd=["srctool", "eventfiles="+input_evtfile, "exttype=tophat", "srccoord=fk5;circle("+str(coord_arr[int(i+downsample/2),int(j+downsample/2),0])+","+str(coord_arr[int(i+downsample/2),int(j+downsample/2),1])+","+str(pixel_rad)+")","srcreg= mask "+filepath+"mask_temporary_"+str(j)+"_"+str(i)+"_"+str(tm)+".fits.gz", "todo= ARF", "insts="+str(tm), "writeinsts="+str(tm),"psftype=NONE", "prefix="+filepath, "xgrid=0.5", "clobber=yes", "extpars="+str(3600*pixel_rad),"suffix= _"+str(j)+"_"+str(i)+"_"+str(tm)]
                #xgrid originally 0.5
            #print(cmd)
            subprocess.run(cmd)

            # remove temporary mask
            os.remove(filepath+"mask_temporary_"+str(j)+"_"+str(i)+"_"+str(tm)+".fits.gz")
            

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_evtfile", required=True, type=str, help="Input calibrated eROSITA eventfile")
    parser.add_argument("--output_path", required=False, type=str, help="Filepath for output files")
    parser.add_argument("--x_lo", required=False, type=int, help="Lower bound of cols analyzed in this step")
    parser.add_argument("--x_hi", required=False, type=int, help="Lower bound of cols analyzed in this step")
    parser.add_argument("--y", required=False, type=int, help="Row analyzed in this step")
    parser.add_argument("--downsample", required=False, type=int, help="Bin parameter. 2 means the map is binned in steps of 2x2=4 pixels.")
    parser.add_argument("--tm", required=False, type=int, help="TM to analyze. Choose 1, 2, 3, 4, or 6.")

    args = parser.parse_args()

    make_arfs(args.input_evtfile,args.output_path,args.x_lo,args.x_hi,args.y,args.downsample,args.tm)

