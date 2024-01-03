#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=18:00:00
#SBATCH --mem=40G
#SBATCH --output=/scratch/Hudson/ASL/vwASL/ab42RandomiseSCoV0_5/logs/randomise_jobarray_%a.out
#SBATCH --job-name=randomise_ab42

#use this with a lower number of permutations to check that code is working
#not sure if the other parallelized script in this folder is working

inFold=/scratch/Hudson/ASL/vwASL/ab42RandomiseSCoV0_5
gmMask="$inFold"/gmMaskThresh.nii.gz
inputCBF="$inFold"/cbfImgs.nii.gz
outFold="$inFold"/outRandomise

mkdir -p "$outFold"

outputFile="$outFold"/out500

#--vxl=-10 will depend on output of previous step (setup_masks in python script) it will say on command line
#will only need --vxl and --vxf if masking out voxels on per-subject basis
#500 perm for results, 5000 for publication
singularity exec /scratch/Hudson/resource/fsl_6.0.5.1.sif randomise -i "$inputCBF" -o "$outputFile" -m "$gmMask" -d "$inFold"/desMask.mat -t "$inFold"/desMask.con --vxl=-10 --vxf="$inFold"/desMask -n 500 -T

