#!/bin/bash
#SBATCH --mail-user=william.h.robb@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --array=0-23
#SBATCH --output=/scratch/Hudson/ASL/vwASL/vwProcessing/processingInR/logs/job_chunkProc__%A_%a_.out

#Specify project directory
projDir=/scratch/Hudson/ASL/vwASL/vwProcessing
scriptsDir="$projDir"/scripts
processingInRDir="$projDir"/processingInR

#Variables to specify
#Longitudinal Dataframe CSV File. Only change the csv filename in caps below.
covFile="$projDir"/sample.csv
#Specify your Merged 4D filename below. Only replace the filename in caps below.
inputFile="$processingInRDir"/cbfImgs.nii.gz
vwRegFile="$processingInRDir"/gmMasks.nii.gz
#Mask img
maskImg="$processingInRDir"/gmMaskThresh.nii


#Other needed variables, directories, and files. Do not modify anything here and below
#Number of total chunks your 4D files were broken into
nChunks=24
inChunkPath="$processingInRDir"/merged_chunks
vwRegChunkPath="$processingInRDir"/vwReg_chunks
outPath="$processingInRDir"/model_chunks

echo "here"

#Uncomment the loop to run locally on badger
#for SLURM_ARRAY_TASK_ID in `seq 0 51`; do
 # echo "PROCESSING CHUNK ${SLURM_ARRAY_TASK_ID}"
singularity exec /gpfs23/scratch/vmac/vwASL_Long.simg Rscript --no-save "$scriptsDir"/chunkProc_vwReg.R "$inputFile" "$inChunkPath" "$vwRegFile" "$vwRegChunkPath" ${SLURM_ARRAY_TASK_ID} "$nChunks" "$maskImg" "$covFile" "$outPath"
#done
