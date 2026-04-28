evtfile_list=#path to txt file list of eventfiles
source_reg=#ds9 region encompassing both fg and bkg data

fg_reg=#source region
bg_reg=#bkg region

tms=[1,2,3,4,6] #TMs to use [1,2,3,4,6] is recommended (5 and 7 are affected by optical light leak)
#dr1=True #True if data is from the DR1 or a scanned observation, False if pointed data
pointed=False #True if pointed EDR data, False if DR1 or scan data
on_off=True#choose if an OnOff bkg is created
wstat_rebin_n=None#rebin to minimum number of off counts per bin for correct wstat treatment (recommended: 5)

output_stacked=True #to get stacked dataset
output_tms=False #True to output a dataset for each TM


#define paths
out_path=
scratch_path=#path for temporary files"/home/wecapstor1/caph/mppi147h/pwn_analysis/DR1_eROSITA_data/MSH15-56/datasets/scratch/"
eROdata_path='../'#path where eROdata is located
catalog_path=#path to eRASS1_Main catalog for point source exclusion

#bin size parameters
#massively affect file size and computation time, choose wisely and with respect to region size
evt_bin=8 #arcsec, ideally multiple of 4
exp_bin=120 #arcsec, use multiple (at least x2) of evt_bin
psf_bin=160 #arcsec, use multiple of evt_bin
srctool_like=False#only necessary for exact agreement with srctool


#more detailed options:

bkg_ext_mask=None#Additional mask to exclude from background?
bkg_exclude_reg=None#Additional region to exclude from background?
bkg_weighting_exp=True#Bkg weighted by exposure time (True) or full exposure (False)? Important for pointed data.
#e_reco_axis="original"#e_reco_axis to use
#e_true_axis="original"#e_true_axis to use
binsize_irf=10#arcmin Change PSF & EDisp binning in dataset