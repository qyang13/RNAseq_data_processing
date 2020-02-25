#!/bin/bash

###############################################################################
# Setting the parameters
###############################################################################
# Set your email to receive notification about the job status
email=jane.doe@colorado.edu

# Specify where the new sequencing data is stored
gz=/Shares/dtra_collab/<path_to_new_sequencing_data>

###############################################################################
# Copy the entire script to Fiji terminal and run :)
###############################################################################

export PATH=~:$PATH
mkdir ~/RNA_proc/
mkdir ~/RNA_proc/fastq
mkdir ~/RNA_proc/output
mkdir ~/RNA_proc/bam
WD=~/RNA_proc
module load bbmap fastqc samtools hisat2 subread
ref=/scratch/Shares/sawyer/ref/hg38
bbmap_adapters=/scratch/Shares/sawyer/ref/adapters.fa

cp_id=$(sbatch --mem=10g -p short --time=23:50:00 --ntasks=8 --nodes=1 --output=./%x_%j.out --error=./%x_%j.err --mail-type=ALL --job-name=cp \
 --mail-user=${email} \
 --wrap="cp ${gz}/*.gz ${WD}/fastq")

cp_id=${cp_id##* }

for name in $(ls ${WD}/fastq/*_R1.fastq.gz)
do
 name="${name%_R1.*.*}"
 proc_id=$(sbatch --dependency=afterok:${cp_id} --mem=20g -p short --time=23:50:00 --ntasks=32 \
  --nodes=1 --output=./%x_%j.out --error=./%x_%j.err \
  --mail-type=ALL --mail-user=${email} --job-name=RNAseq_proc \
  --wrap="gunzip ${name}*;
         bbduk.sh -Xmx20g \
           t=32 \
           in=${name}_R1.fastq \
           in2=${name}_R2.fastq \
           out=${name}_R1.trim.fastq \
           out2=${name}_R2.trim.fastq \
           ref=${bbmap_adapters} \
           ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
           maq=10 minlen=25 \
           tpe tbo \
           literal=AAAAAAAAAAAAAAAAAAAAAAA \
           stats=${name}.trimstats.txt \
           refstats=${name}.refstats.txt \
           ehist=${name}.ehist.txt;
         fastqc ${name}_R1.trim.fastq;
         fastqc ${name}_R2.trim.fastq;
         mv ${name}*.html ${WD}/output;
         hisat2  -p 32 \
           --very-sensitive \
           -x ${ref}/HISAT2/hg38_ercc \
           --pen-noncansplice 14 \
           --rna-strandness RF \
           --mp 1,0 \
           --sp 2,1 \
           -1 ${name}_R1.trim.fastq \
           -2 ${name}_R2.trim.fastq \
           --new-summary \
           > ${name}.sam \
           2> ${name}.hisat2_mapstats.txt;
         mv ${name}.*.txt ${WD}/output;
         samtools view -@ 32 -bS -o ${name}.bam ${name}.sam;
         samtools sort -@ 32 ${name}.bam > ${name}.sorted.bam;
         samtools index ${name}.sorted.bam ${name}.sorted.bam.bai;
         mv ${name}.sorted.bam* ${WD}/bam;
         rm *trim* *.zip *.sam *.bam")
    proc_id=${proc_id##* }
    proc_ids=${proc_ids}:${proc_id}
done

sbatch --dependency=afterok${proc_ids} --mem=20g -p short --time=23:50:00 --ntasks=32 \
 --nodes=1 --output=./%x_%j.out --error=./%x_%j.err \
 --mail-type=ALL --mail-user=${email} --job-name=featureCounts \
 --wrap="featureCounts -T 16 \
         -a ${ref}/hg38_ercc.refseq.gtf \
         -o ${WD}/output/hg38_ercc_readcount.txt \
         -s 2 \
         ${WD}/bam/*.sorted.bam"
