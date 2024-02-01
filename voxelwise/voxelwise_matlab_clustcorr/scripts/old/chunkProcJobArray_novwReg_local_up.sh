#!/bin/bash
#SBATCH --mail-user=william.h.robb@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --array=0-23
#SBATCH --output=/scratch/Hudson/voxelwise/wmhBurdenDelay/processingInR/logs/job_chunkProc__%A_%a_.out

#Specify project directory
projDir=/nobackup/h_vmac/user/robbwh/voxelwise/fdg_pet/cog_status
scriptsDir="$projDir"/scripts
processingInRDir="$projDir"/processingInR

#Variables to specify
#Longitudinal Dataframe CSV File. Only change the csv filename in caps below.
covFile="$projDir"/sample.csv
#Specify your Merged 4D filename below. Only replace the filename in caps below.
inputFile="$processingInRDir"/pet_imgs_4d.nii.gz
vwRegFile="$processingInRDir"/Merged_lagmasks.nii.gz
#Mask img
maskImg="$processingInRDir"/pet_img_mean_mask_b4.nii


#Other needed variables, directories, and files. Do not modify anything here and below
#Number of total chunks your 4D files were broken into
nChunks=31
inChunkPath="$processingInRDir"/merged_chunks
vwRegChunkPath="$processingInRDir"/vwReg_chunks
outPath="$processingInRDir"/model_chunks

#Uncomment the loop to run locally on badger
for SLURM_ARRAY_TASK_ID in `seq 20 -1 15`; do
  echo "PROCESSING CHUNK ${SLURM_ARRAY_TASK_ID}"
singularity exec /nobackup/h_vmac/vwASL_Long.simg Rscript --no-save "$scriptsDir"/chunkProc_novwReg.R "$inputFile" "$inChunkPath" ${SLURM_ARRAY_TASK_ID} "$nChunks" "$maskImg" "$covFile" "$outPath"
done
