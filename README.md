# RNA-seq data processing
This tutorial introduces steps to download and process **human RNA-seq data** generated using Ilumina NovaSeq platform at CU-Anschutz.

**NOTE: Unless noted otherwise, all of the following steps are to be executed on Fiji.**

## Table of Contents:
* [Step 1: Retrieve sequencing data and setup](#step1)
* [Step 2: Trim reads, QC, Align & Read Count](#step2)
* [Step 3: Download the read counting results](#step3)
* [Step 4: Download the processed results](#step4)

<H2 id="step1">Step 1: Retrieve sequencing data and setup</H2>
Typically, once Halley finishes her sample preparation, she will be responsible to send the samples to CU-Anschutz for sequencing. The sequencing generally takes up to 3 weeks, with library preparation taking additional 1-2 weeks. The new sequencing data will be stored under Halley's account on CU-Anschutz FTP server, so she is responsible for transferring the new data to Fiji.

**NOTE: If the sequencing is not paid via DTRA grant, make sure to add a different speedtype and PI name in the sample submission form. Also, to avoid delay, add someone from Sawyer Lab in the sample submission form so we are notified as soon as the sequencing is complete.**

Once Halley notifies us that the data has been transferred to Fiji, prepare a copy of the files in your home directory. Don't forget to change the `<path to new sequencing data>` and `<your email>` below:
```
# Set your email to receive notification about the job status
email=jane.doe@colorado.edu

# Specify where the new sequencing data is stored
gz=/Shares/dtra_collab/<path to new sequencing data>

# Set up the environment:
export PATH=~:$PATH

# Create a few new directories to store all processing data
mkdir ~/RNA_proc/
mkdir ~/RNA_proc/fastq
mkdir ~/RNA_proc/output
mkdir ~/RNA_proc/bam

# Set the working directory
WD=~/RNA_proc

# Copy the raw fastq files to your own working directory, this might take a few minutes
sbatch --mem=10g -p short --time=23:50:00 --ntasks=8 --nodes=1 --output=./%x_%j.out --error=./%x_%j.err --mail-type=ALL --job-name=cp \
 --mail-user=${email} \
 --wrap="cp ${gz}/*.gz ${WD}/fastq"
```

Once the job is complete, change into your working directory to make sure the files are transferred correctly.

<H2 id="step2">Step 2: Trim reads, QC, Align & Read Count</H2>
Once the sequencing data is copied successfully, we will first need to tell Fiji where everything is, this includes the sequencing data, the reference genome, and where all softwares are. Once everything is setup, simply copy the following code snippet into Fiji terminal, and wait until the job completes :)

```
module load bbmap fastqc samtools hisat2 subread
ref=/scratch/Shares/sawyer/ref/hg38
bbmap_adapters=/scratch/Shares/sawyer/ref/adapters.fa

for name in $(ls ${WD}/fastq/*_R1.fastq.gz)
do
    name="${name%_R1.*.*}"
    sbatch --mem=20g -p short --time=23:50:00 --ntasks=32 \
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
            rm *trim* *.zip *.sam *.bam"
done
```
The whole process should take around 1-2 hours to complete, once the job is completed, it's always a good idea to check the output directory to make sure all the files, with suffixes `.html`, `.txt` and `.sorted.bam` are there, and its filesize is not zero. This can be done using the command:
```
ls -lsh ${WD}/output
```

<H2 id="step3">Step 3: Count reads for each gene</H2>
By this step, all the sequencing reads have gone through quality filtering, trimming, QC, mapping and format conversion. To finally figure out how many reads are mapped onto each gene, we will need to run one more command. Since this command finishes within 2 minutes, we will simply run it through Fiji head node. No matter which directory you are in, simply copy and paste the following code into the terminal and run:

```
ref=/scratch/Shares/sawyer/ref/hg38
WD=~/RNA_proc
module load subread

featureCounts -T 16 \
-a ${ref}/hg38_ercc.refseq.gtf \
-o ${WD}/output/hg38_ercc_readcount.txt \
-s 2 \
${WD}/bam/*.sorted.bam
```
Once the process completes, there should be a `hg38_ercc_readcount.txt` file containing all the read counts as well as a `hg38_ercc_readcount.txt.summary` file summarizing the total number of reads mapped/discarded.

<H2 id="step4">Step 4: Download the processed results</H2>
Once the program finishes, all the files needed to evaluate the success of the processing and for downstream analysis are all going to be stored in `~/RNA_proc/output`. The files should be small in size (< 50 MB total) so it can be easily downloaded and transferred. To download the files, the fastest way would require you first **logout from Fiji**, then type in the following command in your local terminal (replace your username below):
```
rsync -r <Your_Username>@fiji.colorado.edu:~/RNA_proc/output ~/Desktop
```
Once the command is complete, the output folder should show up on your Desktop.
