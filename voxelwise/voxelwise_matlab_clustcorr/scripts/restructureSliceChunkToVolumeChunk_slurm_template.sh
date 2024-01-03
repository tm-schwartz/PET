#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:30:00
#SBATCH --mem=64G
#SBATCH --output=/projDir/logs/restructureSliceChunkToVolumeChunk_%A.out

#Specify project directory
#Project Directory Path
ml GCC/8.2.0  OpenMPI/3.1.4
ml R/3.6.0

projDir=/nobackup/h_vmac/user/robbwh/voxelwise/fdg_pet/pre_post
#Path to scripts folder
scriptsDir="$projDir"/scripts
#Path to R processing directory
processingInRDir="$projDir"/processingInR

#Variables to specify
#Specify your Merged 4D filename below. Only replace the filename in caps below.
inputFile="$processingInRDir"/pet_imgs_4d.nii.gz

#Paths to other variables and directories. Do not modify anything here and below.
chunkPath="$processingInRDir"/reconstructed_images
outPath="$processingInRDir"/restructured_images

#Number of standardized residual chunks (same as the number of model chunks and number of chunks your initial 4D images were broken into)
nChunks=31
#Number of volumes in each output chunk
nVolumesPerOutput=25

#singularity exec /scratch/vmac/vwASL_Long.simg Rscript --no-save "$scriptsDir"/restructureSliceChunkToVolumeChunk.R "$inputFile" "$chunkPath" "$nChunks" "$nVolumesPerOutput" "$outPath"
Rscript --no-save "$scriptsDir"/restructureSliceChunkToVolumeChunk.R "$inputFile" "$chunkPath" "$nChunks" "$nVolumesPerOutput" "$outPath"




