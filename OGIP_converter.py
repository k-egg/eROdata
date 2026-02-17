import gammapy
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from regions import SkyRegion

# %matplotlib inline
import matplotlib.pyplot as plt
from IPython.display import display
from gammapy.data import EventList
from gammapy.datasets import Datasets, MapDataset
from gammapy.irf import EDispKernelMap, PSFMap,EDispKernel
from gammapy.maps import Map, MapAxis, WcsGeom,RegionNDMap,RegionGeom,MapAxes
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    Models,
    PointSpatialModel,
    PowerLawNormSpectralModel,
    PowerLawSpectralModel,
    SkyModel,
    TemplateSpatialModel,
    create_fermi_isotropic_diffuse_model,
)
import numpy as np
from os import path
import os, os.path, time, subprocess


def OGIP_converter(pha_file,suffix=None):
    if suffix==None:
        suffix=""

    # open counts spectrum, get filepath, get ARF & RMF filenames
    pha_hdulist = fits.open(pha_file)
    pha1=pha_hdulist[1]
    outfile_pha=pha_file.replace(".fits","")+suffix+".fits"
    filepath=pha_file[0:pha_file.rfind("/")]
    
    rmf_file=pha1.header["RESPFILE"]
    arf_file=pha1.header["ANCRFILE"]

    # define backgorund file output path if necessary
    if "BACKFILE" in pha1.header: 
        try:
            bkg_file=pha1.header["BACKFILE"]
            bkg_hdulist = fits.open(bkg_file)
            bkg1=bkg_hdulist[1]
            outfile_bkg=filepath+bkg_file[bkg_file.rfind("/"):].replace(".fits","")+suffix+".fits"
            pha1.header["BACKFILE"]=outfile_bkg
        except:
            pha1.header.remove("BACKFILE")
    
    rmf_hdulist = fits.open(rmf_file)
    
    rmf1=rmf_hdulist[1]
    rmf_ebounds=rmf_hdulist["EBOUNDS"]
    
    ebounds=rmf_hdulist["EBOUNDS"].data
    outfile_rmf=filepath+rmf_file[rmf_file.rfind("/"):].replace(".fits","")+suffix+".fits"
    pha1.header["RESPFILE"]=outfile_rmf
    #pha1.header["ANCRFILE"]=arf_file
    
    # Adjust counts spectrum file:
    
    #BACKSCAL keyword to column
    if "BACKSCAL" in pha1.columns.names:
        pass
    elif "BACKSCAL" in pha1.header:
        backscal_array = [pha1.header["BACKSCAL"]*pha1.header["EXPOSURE"] for i in range(len(pha1.data))]
        backscal = fits.Column(name="BACKSCAL",format="1D",array=backscal_array)
        pha1.columns.add_col(backscal)
    else:
        raise KeyError("no BACKSCAL column or keyword in pha file")

    # get quality and grouping information if present
    if "QUALITY" in pha1.columns.names:
        if "GROUPING" in pha1.columns.names:
        #Grouping
            gr=pha1.columns["GROUPING"].array
            qu=pha1.columns["QUALITY"].array

            bounds=[]
            for i, val in enumerate(gr):
                if val==1 and qu[i]==0:
                    bounds.append(rmf_ebounds.data[i][1])

                if val==-1 and qu[i]==2:
                    bounds.append(rmf_ebounds.data[i-1][2])
                    break
    else:
        qual_array = np.zeros(len(pha1.data))
        qual = fits.Column(name="QUALITY",format="1D",array=qual_array)
        pha1.columns.add_col(qual)
    
    # Add EBOUNDS extension from RMF file to counts spectrum and write output:
    pha_hdulist.append(rmf_ebounds)
    pha_hdulist.writeto(outfile_pha,overwrite=True)
    
    # Adjust background file:
    
    # Change BACKSCAL keyword to column & adjust for weighting with exposure time:
    if "BACKFILE" in pha1.header:
        if "BACKSCAL" in bkg1.columns.names:
            pass
        elif "BACKSCAL" in bkg1.header:
            backscal_array = [bkg1.header["BACKSCAL"]*bkg1.header["EXPOSURE"] for i in range(len(bkg1.data))]
            backscal = fits.Column(name="BACKSCAL",format="1D",array=backscal_array)
            bkg1.columns.add_col(backscal)
        else:
            raise Exception("no BACKSCAL column or keyword in bkg file")
    
    # Add EBOUNDS extension from RMF file to background spectrum and write output:
        bkg_hdulist.append(rmf_ebounds)
        bkg_hdulist.writeto(outfile_bkg,overwrite=True)
    
    
    # RMF file:
    # reduce F_CHAN column (in eROSITA data) from 1 to 0 to ensure proper indexing
    rmf_matrix_data=rmf1.data
    for i in range(len(rmf_matrix_data["F_CHAN"])):
        k=rmf_matrix_data["F_CHAN"][i]
        rmf_matrix_data["F_CHAN"][i]=[k[0]-1]

    # write RMF output
    fits.writeto(outfile_rmf,data=rmf_matrix_data,overwrite=True,header=rmf1.header)
    fits.append(outfile_rmf,data=ebounds,header=rmf_ebounds.header)
    
    print("Files successfully converted!")
    
    #Definition of rebinned energy axis:
    if bounds:
        axiis = MapAxis.from_edges(bounds,name="energy",unit="keV")
        return axiis
    else:
        pass
