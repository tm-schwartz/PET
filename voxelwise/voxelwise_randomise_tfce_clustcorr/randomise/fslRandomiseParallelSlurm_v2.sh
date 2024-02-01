#!/bin/sh

function Usage {

cat <<USAGE
This script takes a design and contrast file and runs randomise parallel on accre for
the vwASL pipeline.

Usage: ./fslRadomiseSlurm.sh -i inputImage -l logDir -o outPrefix -d designMatFile -c contrastFile -m maskFile -n nPerm
  -i: inputImages: Image to be
  -l: LogDirectory: Directory where slurm job files will be kept
  -o: outputImage: output file prefix or path/prefix
  -d designMatFile: design matrix mat file
  -c contrastFile: contrast matrix con file
  -m maskFile: image to be used as mask
  -n nPerm: number of permutations
USAGE
  exit 1
}

if [[ $# -ne 18 ]];
  then
  Usage >&2
fi

#Set flags
while getopts ":i:l:o:d:c:m:n:v:f:" opt; do
  case $opt in
    i)
      inputImage=$OPTARG
      ;;
    l)
      logDir=$OPTARG
      ;;
    o)
      outPrefix=$OPTARG
      ;;
    d)
      desFile=$OPTARG
      ;;
    c)
      conFile=$OPTARG
      ;;
    m)
      maskFile=$OPTARG
      ;;
    n)
      nPerm=$OPTARG
      ;;
    v)
      vxl=$OPTARG
      ;;
    f)
      vxf=$OPTARG
      ;;
  esac
done

#Check if any of the flags are empty
if [[ -z $inputImage ]] || [[ -z $outPrefix ]] || [[ -z $logDir ]] || [[ -z $desFile ]] || [[ -z $conFile ]] || [[ -z $maskFile ]] || [[ -z $nPerm ]] || [[ -z $vxl ]] || [[ -z $vxf ]];
  then
  Usage >&2
fi

time_start=`date +%s`

#specify slurm script
timestamp=`date +%m-%d-%Y-%H-%M-%S`
mkdir -p ${logDir}
qscript=${logDir}/job_fslRandomise_${timestamp}.slurm
echo '#!/bin/sh' > $qscript
echo '#SBATCH --nodes=1' >> $qscript
echo '#SBATCH --ntasks=1' >> $qscript
echo '#SBATCH --cpus-per-task=8' >> $qscript
echo '#SBATCH --time=18:00:00' >> $qscript
echo '#SBATCH --mem-per-cpu=25G' >> $qscript
echo "#SBATCH --output=${logDir}/fslRandomise_%j.out" >> $qscript
echo 'fsl6=/scratch/Hudson/resource/fsl_6.0.5.1.sif' >> $qscript
echo 'export FSLPARALLEL=8' >> $qscript
echo 'echo "Running FSL randomise"' >> $qscript
echo 'singularity exec $fsl6 randomise_parallel -i '"${inputImage} -o ${outPrefix} -d ${desFile} -t ${conFile} -m ${maskFile} -n ${nPerm} --vxl=-${vxl} --vxf=${vxf} -T" >> $qscript
echo 'if [[ $? -ne 0 ]];' >> $qscript
echo '  then' >> $qscript
echo '    echo "Error running randomise_parallel"' >> $qscript
echo '    exit 1' >> $qscript
echo 'fi' >> $qscript
echo 'echo "Finished fsl randomise"' >> $qscript
echo 'exit 0' >> $qscript

echo
echo "--------------------------------------------------------------------------------------"
echo "Submitting fslRandomise job will report back when no job remain"
echo "--------------------------------------------------------------------------------------"

sbatch --wait ${qscript}

echo
echo "--------------------------------------------------------------------------------------"
echo "No job remains."
echo "--------------------------------------------------------------------------------------"

if [[ ! -f ${outPrefix}_tstat1.nii.gz ]] || [[ ! -f ${outPrefix}_tstat2.nii.gz ]] || [[ ! -f ${outPrefix}_tfce_corrp_tstat1.nii.gz ]] || [[ ! -f ${outPrefix}_tfce_corrp_tstat2.nii.gz ]];
  then
    echo "One or more randomise outputs not found. Check success of SLURM Job"
    exit 1
fi

time_end=`date +%s`
time_elapsed=$((time_end - time_start))
echo
echo "--------------------------------------------------------------------------------------"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

exit 0
