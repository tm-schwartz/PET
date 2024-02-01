#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=03:30:00
#SBATCH --mem=70G
#SBATCH --output=/projDir/processingInR/logs/Chunk4DNifti_vwReg_job_%A.out

#Specify project directory
#Project Directory Path
projDir=/scratch/Hudson/ASL/vwASL/vwProcessing
#Path to scripts folder
scriptsDir="$projDir"/scripts
#Path to R processing directory
processingInRDir="$projDir"/processingInR

#Variables to specify
#Specify your Merged 4D filename below. Only replace the filename in caps below.
inputFile="$processingInRDir"/gmMasks.nii.gz

#Paths to other variables and directories. Do not modify anything here and below.
#Output path
outPath="$processingInRDir"/vwReg_chunks
#Chunk Thickness
chunkThickness=4

#Script call
singularity exec /scratch/vmac/vwASL_Long.simg Rscript --no-save "$scriptsDir"/chunk4DNifti.R "$inputFile" "$chunkThickness" "$outPath"
