#!/bin/bash

#Copy scripts for voxelwise analysis

#call using "./prepareAnalysis.sh ."

#The single arguement to script is the desired project directory
inPath=$1


in_scriptsDir=/nobackup/h_vmac/user/robbwh/resource/gitS/PET/voxelwise/voxelwise_matlab_clustcorr/scripts

#Put needed scripts in scripts directory
scriptsDir="$inPath"/scripts
mkdir "$scriptsDir" 

cp $in_scriptsDir/* $scriptsDir

cp $scriptsDir/groupProc3.py "$inPath"

#Create directory structure for R processing
procInRDir="$inPath"/processingInR
mkdir "$procInRDir"
mkdir "$procInRDir"/merged_chunks
mkdir "$procInRDir"/vwReg_chunks
mkdir "$procInRDir"/model_chunks
mkdir "$procInRDir"/reconstructed_images
mkdir "$procInRDir"/restructured_images
mkdir "$procInRDir"/spm_results
mkdir "$procInRDir"/logs
