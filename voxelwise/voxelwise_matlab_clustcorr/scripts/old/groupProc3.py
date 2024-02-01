import nibabel as nib
from nilearn.image import math_img,load_img
import pandas as pd 
import numpy as np
import os
import sys


script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
outfold=script_directory+"/"


icv_mask_loc = "/nobackup/h_vmac/user/robbwh/voxelwise/fdg_pet/mask_ICV_resampled.nii.gz"
icv_mask = load_img(icv_mask_loc)
 
data_fold = "/nobackup/h_vmac/user/robbwh/voxelwise/fdg_pet/data/"
subjsCsv = "/nobackup/h_vmac/user/robbwh/voxelwise/fdg_pet/data/newcovariates.csv"

subjs = pd.read_csv(subjsCsv)



subjs['pre_post'] = [0, 1] * (len(subjs) // 2) + [0] * (len(subjs) % 2)

print(subjs)


pt_ids = subjs["Patient"].to_list()


file_locs = []

for i,pt in enumerate(pt_ids):

  if i%2==0:
    subjs.at[i,'DeltaPost'] = 0
    suff="sSUV_pre.nii"
    outloc = data_fold+pt+"/crop/"+suff
  else:
    suff="sSUV_post.nii"
    outloc = data_fold+pt+"/crop/"+suff

  file_locs.append(outloc)


#print(file_locs)



petImgsComb = nib.concat_images(file_locs)

###normalization
petImgsData = petImgsComb.get_fdata()
imgdims = petImgsData.shape[3]


#masksSubjs = np.zeros(petImgsData.shape)
petSubjs_SUVR = np.zeros(petImgsData.shape)

#http://spect.yale.edu/analysis_details.html
#trying this now for the scaling
for imgdim in range(imgdims):
  print(imgdim)
  #get global mean
  img_mean = np.nanmean(petImgsData[:,:,:,imgdim])
  #get regions where intensity is above mean value divided by 8
  img_excl = np.where(petImgsData[:,:,:,imgdim]>(img_mean/8),petImgsData[:,:,:,imgdim],np.nan)
  #now get new global mean
  img_excl_mean = np.nanmean(img_excl)
  #now threshold by (0.8*img_excl_mean)
  img_thr = np.where(petImgsData[:,:,:,imgdim]>(img_excl_mean*0.8),petImgsData[:,:,:,imgdim],np.nan)
  #scale by the new global mean
  petSubjs_SUVR[:,:,:,imgdim] = img_thr/img_excl_mean




##remask zeros using nan for stats (which need to exclude nans..)
petSubjs_SUVR = nib.Nifti1Image(petSubjs_SUVR,petImgsComb.affine,petImgsComb.header)
#mask with nan after image creation.. (looks like image creation will reintroduce 0s where nans were due to header issues)
petSubjs_SUVR = math_img("np.where(petSubjs_SUVR>0.000001,petSubjs_SUVR,np.nan)",petSubjs_SUVR=petSubjs_SUVR)
petSubjs_SUVR.to_filename(outfold+"pet_imgs_4d.nii.gz")


##binarize the images to make a mask
petSubjs_bin = math_img("np.where(np.isnan(petSubjs_SUVR),0,1)",petSubjs_SUVR=petSubjs_SUVR)
petSubjs_bin = math_img("np.mean(petSubjs_bin,axis=3)",petSubjs_bin=petSubjs_bin)

#multiply this by the icv mask
petSubjs_bin = math_img("petSubjs_bin*icv_mask",petSubjs_bin=petSubjs_bin,icv_mask=icv_mask)

#keep areas where at least 98% of people have data within analysis mask
petSubjs_bin = math_img("np.where(petSubjs_bin>0.98,1,0)",petSubjs_bin=petSubjs_bin)

petSubjs_bin.to_filename(outfold+"explicit_mask.nii.gz")

#petImgs_mean = math_img("np.mean(img,axis=3)",img=petImgsComb)
#petImgs_mean.to_filename(outfold+"pet_imgs_mean.nii.gz")


subjs.to_csv(outfold+"sample.csv",index=False)

