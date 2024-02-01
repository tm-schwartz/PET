#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:30:00
#SBATCH --mem=16G
#SBATCH --output=/projDir/logs/reconstructImgFromModelChunksjob_%A.out

ml GCC/8.2.0  OpenMPI/3.1.4 R/3.6.0

#Project Directory Path
projDir=/nobackup/h_vmac/user/robbwh/voxelwise/fdg_pet/pre_post
scriptsDir="$projDir"/scripts
processingInRDir="$projDir"/processingInR

#Variables to specify
#Specify your Merged 4D filename below. Only replace the filename in caps below.
inputFile="$processingInRDir"/pet_imgs_4d.nii.gz
#Residual (error) degrees of freedom of your desired model term
dof=66

#Paths to other variables and directories. Do not modify anything here and below.
chunkPath="$processingInRDir"/model_chunks
outPath="$processingInRDir"/reconstructed_images
#Number of total model chunks (same as the number of chunks your 4D file was broken into)
nChunks=31

#singularity exec /scratch/vmac/r_pandoc_rstudio.simg Rscript --no-save "$scriptsDir"/reconstructImgFromModelChunks.R "$inputFile" "$chunkPath" "$nChunks" "$dof" "$outPath"
Rscript --no-save "$scriptsDir"/reconstructImgFromModelChunks.R "$inputFile" "$chunkPath" "$nChunks" "$dof" "$outPath"
