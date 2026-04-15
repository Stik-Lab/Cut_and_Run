#!/bin/sh
#SBATCH --job-name=CUTNRUN
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=CutnrunPIPELINE_%A-%a.log

# ========== VARIABLES ==========
# put in a file call samples.txt the name of the variables

N=$(wc -l < samples.txt)
describer=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

source ./config.sh
for dir in "${path_bam}" "${path_temp}" "${path_bw}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done

# ========== MODULES ==========
module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m
module load Trim_Galore/0.6.6-foss-2021b-Python-3.8.5
module load Bowtie2/2.4.4.1-GCC-11.2.0
module load samtools-1.12-gcc-11.2.0-n7fo7p2
module load picard/2.26.3-Java-11
module load bedtools2-2.30.0-gcc-11.2.0-q7z4zez
module load deepTools/3.5.1-foss-2021b

# ========== FASTQC ==========
echo "................................................................ 1. START_FASTQC ${describer} ................................................................"

fastqc ${path_fq}/${describer}_*.fastq.gz -o ${path_fq}

echo "................................................................ 1. END_FASTQC ${describer} ................................................................"


# ========== TRIMMING ==========
echo   "................................................................     2. START_TRIMMING ${describer}   ................................................................"

if [ "${paired}" == "no" ]; then

    R1=$(ls ${path_fq}/${describer}*.f*q.gz | head -n 1)
    trim_galore --fastqc --output_dir ${path_fq} "${R1}"

else

    R1=$(ls ${path_fq}/${describer}*{_1,_R1}*.f*q.gz | head -n 1)
    R2=$(ls ${path_fq}/${describer}*{_2,_R2}*.f*q.gz | head -n 1)
    trim_galore --fastqc --output_dir ${path_fq} --paired "${R1}" "${R2}"

fi


echo    "................................................................    2. END_TRIMMING ${describer}   ................................................................"


# ========== ALIGNMENT ==========
echo "................................................................ 3. START_ALIGNMENT ${describer} ................................................................"

if [ "${paired}" == "no" ]; then

T1=$(ls ${path_fq}/${describer}*_trimmed.fq.gz 2>/dev/null || ls ${path_fq}/${describer}*.fq.gz | head -n 1)

    bowtie2 -x ${indexgenome} \
           -U "${T1}" \
           --very-sensitive-local --no-unal -k 2 --phred33 -I 10 -X 700 -p 8 \
           -S ${path_temp}/${describer}.sam


else

    T1=$(ls ${path_fq}/${describer}*val_1.fq.gz)
    T2=$(ls ${path_fq}/${describer}*val_2.fq.gz)

    bowtie2 -x ${indexgenome} \
           -1 "${T1}" \
           -2 "${T2}" \
           --very-sensitive-local --no-unal --no-mixed --no-discordant -k 2 --phred33 -I 10 -X 700 --dovetail -p 8 \
           -S ${path_temp}/${describer}.sam

fi

echo "................................................................ 3. END_ALGINMENT ${describer} ................................................................"


echo "................................................................ 4. START_SAMtoBAM ${describer} ................................................................"

samtools view -@ 8 -bS ${path_temp}/${describer}.sam > ${path_bam}/${describer}.bam

echo "................................................................ 4. END_SAMtoBAM ${describer} ................................................................"

# ========== FILTERING ==========

echo "................................................................ 5. START_remove_mtDNA_READS ${describer} ................................................................"

samtools view -h ${path_bam}/${describer}.bam | grep -v chrM | samtools sort -O bam -o ${path_temp}/${describer}.rmChrM.bam -T ${path_temp}

echo "................................................................ 5. END_remove_mtDNA_READS ${describer} ................................................................"

echo "................................................................ 6. START_filtering_protperly_reads ${describer} ................................................................"

samtools view -F 2304 -b -q 10  ${path_temp}/${describer}.rmChrM.bam > ${path_temp}/${describer}.qual.bam

# -b: sortida en bam
# -q 10: filtra les lectures si es menor a 10 fora

echo "................................................................ 6. END_filtering_protperly_reads ${describer} ................................................................"

echo "................................................................ 7. START_mark_duplicates ${describer} ................................................................"

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=${path_temp}/${describer}.qual.bam \
O=${path_temp}/${describer}_removed_duplicates.bam \
M=${path_temp}/${describer}_marked_dup_metrics.txt \
REMOVE_DUPLICATES=true ASSUME_SORTED=true VERBOSITY=WARNING

echo "................................................................ 7. END_mark_duplicates ${describer} ................................................................"

echo "................................................................ 8. START_blacklist regions + index ${describer} ................................................................"

bedtools intersect -nonamecheck -v -abam ${path_temp}/${describer}_removed_duplicates.bam -b ${blacklist} > ${path_bam}/${describer}_clean.bam
samtools index ${path_bam}/${describer}_clean.bam

echo "................................................................ 8. END_blacklist + index ${describer} ................................................................"

# ========== COVERAGE ==========

echo "................................................................ 9. START_bigwig ${describer} ................................................................"

bamCoverage --bam ${path_bam}/${describer}_clean.bam --outFileName ${path_bw}/${describer}.bw --effectiveGenomeSize ${effectiveGenomeSize} --outFileFormat bigwig --binSize 10 --normalizeUsing RPGC > ${path_bw}/${describer}.log

echo "................................................................ 9. END_bigwig ${describer} ................................................................"

if [ -f "${path_bw}/${describer}.bw" ]; then
    echo "job successful"
else
    echo "job failed: bigwig not found for ${describer}"
    exit 1
fi

# ==========  LAUNCH ANALYSIS SCRIPTS ==========

if [ "$(grep 'job successful' CutnrunPIPELINE_*.log | wc -l)" -eq "${N}" ]; then
    sbatch --array=1-${N} pipeline/peak_calling.sh
else
    echo "Number of completed jobs: $(grep 'job successful' CutnrunPIPELINE_*.log | wc -l)"
fi
