#!/bin/sh
#SBATCH --job-name=CUTNRUN_pkcalling
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=CutnrunPkCalling_%A-%a.txt

# ========== VARIABLES ==========
# put in a file call samples.txt the name of the variables

N=$(( $(wc -l < samples.txt) - 1 ))
describer=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

module load MACS2/2.2.5-foss-2021b-Python-3.8.5

source ./config.sh
for dir in "${path_macs2}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done

for num in 1 05
do

if ["${narrow}"=="yes" ];then

  if ["${paired}"=="no" ]; then
    macs2 callpeak --format BAM -t ${path_bam}/${describer}_clean.bam -c ${path_bam}/${input}_clean.bam -g hs -n ${describer}_0${num} -q 0.${num} --outdir ${path_macs2}
  else
    macs2 callpeak --format BAMPE -t ${path_bam}/${describer}_clean.bam -c ${path_bam}/${input}_clean.bam -g hs -n ${describer}_0${num} -q 0.${num} --outdir ${path_macs2}
  fi


else

  if ["${paired}"=="no" ]; then
    macs2 callpeak --format BAM -t ${path_bam}/${describer}_clean.bam -c ${path_bam}/${input}_clean.bam -g hs -n ${describer}_0${num} --broad -B -q 0.${num} --outdir ${path_macs2}
  else 
    macs2 callpeak --format BAMPE -t ${path_bam}/${describer}_clean.bam -c ${path_bam}/${input}_clean.bam -g hs -n ${describer}_0${num} --broad -B -q 0.${num} --outdir ${path_macs2}
  fi


fi

done
