evtfile_list= #path to txt file list of eventfiles
source_reg="fk5;circle(228.5338154,-59.1935136,0.7)" #ds9 region encompassing both fg and bkg data

fg_reg="fk5;circle(228.5338154,-59.1935136,0.4895308)" #source region
bg_reg="fk5;circle(228.5338154,-59.1935136,0.7)" #bkg region

tms=[1,2,3,4,6] #TMs to use [1,2,3,4,6] is recommended (5 and 7 are affected by optical light leak)
pointed=False #True if pointed EDR data, False if DR1 or scan data
on_off=True #choose if an OnOff bkg is created
wstat_rebin_n=5 #rebin to minimum number of off counts per bin for correct wstat treatment (recommended: 5)

output_stacked=True #to get stacked dataset
output_tms=False #True to output a dataset for each TM


#define paths
out_path= #output path for data products and datasets
scratch_path= #path for temporary files
eROdata_path= #path where eROdata is located
catalog_path= #path to eRASS1_Main catalog for point source exclusion

#bin size parameters
#massively affect file size and computation time, choose wisely and with respect to region size
evt_bin=16 #arcsec, ideally multiple of 4
exp_bin=48 #arcsec, use multiple (at least x2) of evt_bin
psf_bin=48 #arcsec, use multiple of evt_bin
psf_downsample_time=1 #downsample factor for the time axis over which the PSFMap is created
srctool_like=False #only necessary for exact agreement with srctool


#more detailed options:

bkg_ext_mask=None #External mask to exclude from background?
bkg_exclude_reg="fk5;circle(228.3789225,-59.0181796,0.0832210)" #Additional region to exclude from background?
bkg_weighting_exp=True #Bkg weighted by exposure time (True) or full exposure (False)? Important for pointed data.
binsize_irf=3 #arcmin, Change PSF & EDisp binning in dataset


#cluster options:
runtime_make_evtfiles = 20 #min
mem_mb_make_evtfiles = 5000 #MB

runtime_make_arfs = 120 #min
mem_mb_make_arfs = 5000 #MB

runtime_make_arf_map = 120 #min
mem_mb_make_arf_map = 20000 #MB

runtime_make_psf_map = 180 #min
mem_mb_make_psf_map = 5000 #MB

runtime_make_datasets = 120 #min
mem_mb_make_datasets = 200000 #MB
