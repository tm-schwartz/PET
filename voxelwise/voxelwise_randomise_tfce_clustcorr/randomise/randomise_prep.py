import nibabel as nib
from nilearn.image import math_img
import pandas as pd 
import os
import sys


#script made to setup files along with mask for regions to be excluded on a per-subject basis
#the "setup_masks" part may not be necessary if there are no missing voxels within analysis mask
#or voxels you want to exclude for a subject

inFold= sys.argv[1]
analysis = sys.argv[2]


fsl="/scratch/Hudson/resource/fsl_6.0.5.1.sif"
inFold="/scratch/Hudson/ASL/vwASL/"


outFold=inFold+analysis

#cmd = "mkdir -p "+outFold
#os.system(cmd)

subjsCsv = inFold+analysis+"/finaldesign"+".csv"

subjs = pd.read_csv(subjsCsv,dtype=str)


mapids = subjs["map.id"].to_list()
sesss = subjs["asl.session.id"].to_list()
print(mapids)

foldBase = inFold+"sessions"

subjFolds = []

for i in range(len(mapids)):
  

  mapid = mapids[i]
  sess = sesss[i]

  #subjFolds.append(foldBase+"/"+sess+"/"+mapid)
  subjFolds.append(foldBase+"/"+sess)


#cbfImgsPVC = [x+"/cbf_pvc_inpainted_template_masked.nii.gz" for x in subjFolds]

cbfImgs = [x+"/cbf_nonpvc_template_masked.nii.gz" for x in subjFolds]

print("need different template for pvc")

gmSegsUnthresh = [x+"/gmSegtemplateNonPVC.nii.gz" for x in subjFolds]

gmSegsThresh = [x+"/gmSegtemplateThreshNonPVC.nii.gz" for x in subjFolds]

gmSegsThreshInv = [x+"/gmSegtemplateThreshNonPVC_inverted.nii.gz" for x in subjFolds]

cmd = "singularity exec "+fsl+" setup_masks "+outFold+"/fsl_design.mat "+outFold+"/fsl_contrast.con "+outFold+ "/desMask "+" ".join(gmSegsThreshInv)
print(cmd)
os.system(cmd)



print("here1")

#for img in gmSegsUnthresh:
#  print(img)
#  tempImg = load_img(img)


gmSegsUnthreshComb = nib.concat_images(gmSegsUnthresh)
print("here2")
gmSegsUnthreshComb = math_img("np.mean(img,axis=3)",img=gmSegsUnthreshComb)
gmSegsUnthreshComb = math_img("np.where(img>0.6,1,0)",img=gmSegsUnthreshComb)
gmSegsUnthreshComb.to_filename(outFold+"/gmMaskThresh.nii.gz")


gmSegsThreshComb = nib.concat_images(gmSegsThresh)
gmSegsThreshComb.to_filename(outFold+"/gmMasks.nii.gz")



#cbfImgsPVC = nib.concat_images(cbfImgsPVC)
#cbfImgsPVC.to_filename(outFold+"/cbfImgsPVC.nii.gz")

#cbfMean = math_img("np.mean(img,axis=3)",img=cbfImgsComb)
#cbfMean.to_filename(outFold+"/cbfMean.nii.gz")

#cbf images have extra 4th dim (maybe only pvc images do)
cbfImgsComb = nib.concat_images(cbfImgs)
cbfImgsComb.to_filename(outFold+"/cbfImgs.nii.gz")


#cbfMean = math_img("np.mean(img,axis=3)",img=cbfImgsComb)

#cbfMean.to_filename(outFold+"/cbfMeanNonPVC.nii.gz")



