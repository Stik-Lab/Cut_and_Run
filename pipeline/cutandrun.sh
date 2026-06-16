#!/bin/sh
#SBATCH --job-name=CUTNRUN
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=CutnrunPIPELINE_%A-%a.log

# ========== VARIABLES ==========
# put in a file call samples.txt the name of the variables

describer=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
source ./config.sh

SAM=${path_temp}/${describer}.sam
BAM_CLEAN=${path_bam}/${describer}_clean.bam
BW=${path_bw}/${describer}.bw

for dir in "${path_bam}" "${path_temp}" "${path_bw}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done

N=$(( $(wc -l < samples.txt) - 1 ))

# ========== MODULES ==========
module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m
module load Trim_Galore/0.6.6-foss-2021b-Python-3.8.5
module load Bowtie2/2.4.4.1-GCC-11.2.0
module load samtools-1.12-gcc-11.2.0-n7fo7p2
module load picard/2.26.3-Java-11
module load bedtools2-2.30.0-gcc-11.2.0-q7z4zez
module load deepTools/3.5.1-foss-2021b

# ========== STEP 1: FASTQC ==========
echo "................................................................ 1. START_FASTQC ${describer} ................................................................"

FQC_DONE=$(ls ${path_fq}/${describer}*_fastqc.zip 2>/dev/null | head -n 1)
if [ -n "${FQC_DONE}" ]; then
    echo "SKIP FastQC: reports already exist for ${describer}"
else
  fastqc ${path_fq}/${describer}_*.fastq.gz -o ${path_fq}
fi

echo "................................................................ 1. END_FASTQC ${describer} ................................................................"


# ========== STEP 2: TRIMMING ==========
echo   "................................................................     2. START_TRIMMING ${describer}   ................................................................"

if [ "${paired}" == "no" ]; then
    R1_TRIM=$(ls ${path_fq}/${describer}*_trimmed.fq.gz 2>/dev/null | head -n 1)
    if [ -s "${R1_TRIM}" ]; then
        echo "SKIP Trim Galore: trimmed file already exists for ${describer}"
    else
      R1=$(ls ${path_fq}/${describer}*.f*q.gz | head -n 1)
      trim_galore --fastqc --output_dir ${path_fq} "${R1}"
    fi
else
    R1_TRIM=$(ls ${path_fq}/${describer}*val_1.fq.gz 2>/dev/null | head -n 1)
    R2_TRIM=$(ls ${path_fq}/${describer}*val_2.fq.gz 2>/dev/null | head -n 1)
    if [ -s "${R1_TRIM}" ] && [ -s "${R2_TRIM}" ]; then
        echo "SKIP Trim Galore: trimmed files already exist for ${describer}"
    else
      R1=$(ls ${path_fq}/${describer}*{_1,_R1}*.f*q.gz | head -n 1)
      R2=$(ls ${path_fq}/${describer}*{_2,_R2}*.f*q.gz | head -n 1)
      trim_galore --fastqc --output_dir ${path_fq} --paired "${R1}" "${R2}"
    fi
fi

echo    "................................................................    2. END_TRIMMING ${describer}   ................................................................"


# ========== STEP 3: ALIGNMENT ==========
echo "................................................................ 3. START_ALIGNMENT ${describer} ................................................................"

if [ -s "${SAM}" ]; then
    echo "SKIP Bowtie2: SAM already exists for ${describer}"
else
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
fi

echo "................................................................ 3. END_ALIGNMENT ${describer} ................................................................"

# ========== STEP 4: SAM TO BAM ==========
echo "................................................................ 4. START_SAMtoBAM ${describer} ................................................................"
if [ -s "${path_bam}/${describer}.bam" ]; then
    echo "SKIP SAM to BAM: BAM already exists for ${describer}"
else
  samtools view -@ 8 -bS ${path_temp}/${describer}.sam > ${path_bam}/${describer}.bam
fi
echo "................................................................ 4. END_SAMtoBAM ${describer} ................................................................"


# ========== STEP 5: REMOVE chrM ==========
echo "................................................................ 5. START_remove_mtDNA_READS ${describer} ................................................................"
if [ -s "${path_temp}/${describer}.rmChrM.bam" ]; then
    echo "SKIP chrM removal: file already exists for ${describer}"
else
  samtools view -h ${path_bam}/${describer}.bam | grep -v chrM | samtools sort -O bam -o ${path_temp}/${describer}.rmChrM.bam -T ${path_temp}
fi
echo "................................................................ 5. END_remove_mtDNA_READS ${describer} ................................................................"

# ========== STEP 6: FILTER LOW-QUALITY READS ==========
echo "................................................................ 6. START_filtering_protperly_reads ${describer} ................................................................"
if [ -s "${path_temp}/${describer}.qual.bam" ]; then
    echo "SKIP quality filter: file already exists for ${describer}"
else
  samtools view -F 2304 -b -q 10  ${path_temp}/${describer}.rmChrM.bam > ${path_temp}/${describer}.qual.bam
fi

# -b: sortida en bam
# -q 10: filtra les lectures si es menor a 10 fora

echo "................................................................ 6. END_filtering_protperly_reads ${describer} ................................................................"

# ========== STEP 7: MARK DUPLICATES ==========
echo "................................................................ 7. START_mark_duplicates ${describer} ................................................................"
if [ -s "${path_temp}/${describer}_removed_duplicates.bam" ]; then
    echo "SKIP MarkDuplicates: file already exists for ${describer}"
else
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=${path_temp}/${describer}.qual.bam \
    O=${path_temp}/${describer}_removed_duplicates.bam \
    M=${path_temp}/${describer}_marked_dup_metrics.txt \
    REMOVE_DUPLICATES=true ASSUME_SORTED=true VERBOSITY=WARNING
fi
echo "................................................................ 7. END_mark_duplicates ${describer} ................................................................"

# ========== STEP 8: BLACKLIST + INDEX ==========
echo "................................................................ 8. START_blacklist regions + index ${describer} ................................................................"
if [ -s "${BAM_CLEAN}" ] && [ -s "${BAM_CLEAN}.bai" ]; then
    echo "SKIP blacklist removal: file already exists for ${describer}"
else
  bedtools intersect -nonamecheck -v -abam ${path_temp}/${describer}_removed_duplicates.bam -b ${blacklist} > ${path_bam}/${describer}_clean.bam
  samtools index ${path_bam}/${describer}_clean.bam
fi
echo "................................................................ 8. END_blacklist + index ${describer} ................................................................"

# ========== STEP 9: BIGWIG ==========

echo "................................................................ 9. START_bigwig ${describer} ................................................................"
if [ -s "${BW}" ]; then
    echo "SKIP bamCoverage: bigwig already exists for ${describer}"
else
    bamCoverage --bam ${path_bam}/${describer}_clean.bam --outFileName ${path_bw}/${describer}.bw --effectiveGenomeSize ${effectiveGenomeSize} --outFileFormat bigwig --binSize 10 --normalizeUsing RPGC > ${path_bw}/${describer}.log
fi
echo "................................................................ 9. END_bigwig ${describer} ................................................................"

if [ -f "${path_bw}/${describer}.bw" ]; then
    echo "job successful"
else
    echo "job failed: bigwig not found for ${describer}"
    exit 1
fi

# ==========  LAUNCH ANALYSIS SCRIPTS ==========

if [ "$(grep 'job successful' CutnrunPIPELINE_${SLURM_ARRAY_JOB_ID}-*.log | wc -l)" -eq "${N}" ]; then
    if [ "${wo_input}" == "no" ]; then
        sbatch --array=1-${N} pipeline/peak_calling.sh
    else
        sbatch --array=1-$((N+1)) pipeline/peak_calling_woinput.sh
    fi
else
    echo "Number of completed jobs: $(grep 'job successful' CutnrunPIPELINE_${SLURM_ARRAY_JOB_ID}-*.log | wc -l)"
fi
