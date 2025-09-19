#  Folder paths 
path_fq='path_to_fastq'
path_bam='path_to_bamfiles'
path_temp='path_to_temporaryfiles'

# Files and programs 
indexgenome='/mnt/beegfs/public/references/index/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as'   # requirements
refgenome='/mnt/beegfs/public/references/genome/human/GRCh38.primary_assembly.genome.fa' 
blacklist='/mnt/beegfs/eferre/bin/files/black_list'

# Parameters 
N=$(wc -l < samples.txt)
paired='yes/no'
narrow='yes/no'
genome='hg38'
input=$(tail -n 1 samples.txt)
effectiveGenomeSize=2913022398

N=$(( $(wc -l < samples.txt) - 1 ))
