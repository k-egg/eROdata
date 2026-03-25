evtfile_list=#path to txt file list of eventfiles
source_reg=#ds9 region encompassing both fg and bkg data

fg_reg=#source region
bg_reg=#bkg region

tms=[1,2,3,4,6] #TMs to use [1,2,3,4,6] is recommended (5 and 7 are affected by optical light leak)
stacked=True #to get stacked dataset
dr1=True #True if data is from the DR1 or a scanned observation, False if pointed data
output_tms=False #True to output a dataset for each TM

#define paths
out_path=
scratch_path=#path for temporary files"/home/wecapstor1/caph/mppi147h/pwn_analysis/DR1_eROSITA_data/MSH15-56/datasets/scratch/"
eROdata_path='../'#path where eROdata is located
catalog_path=#path to eRASS1_Main catalog for point source exclusion

#bin size parameters
#massively affect file size and computation time, choose wisely and with respect to region size
evt_bin=8 #arcsec
exp_bin=120 #arcsec
psf_bin=160 #arcsec
srctool_like=False