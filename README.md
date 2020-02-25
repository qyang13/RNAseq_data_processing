# RNA-seq data processing
This tutorial introduces steps to download and process **human RNA-seq data** generated using Ilumina NovaSeq platform at CU-Anschutz.

### Table of Contents:
* [Step 1: Retrieve sequencing data from CU-Anschutz](#step1)
* [Step 2: Trim reads, QC, Align & Read Count](#step2)
* [Step 3: Download the read counting results](#step3)

### Step 1: Retrieve sequencing data from CU-Anschutz
Typically, once Halley finishes her sample preparation, she will be responsible to send the samples to CU-Anschutz for sequencing. The sequencing generally takes up to 3 weeks, with library preparation taking additional 1-2 weeks. The new sequencing data will be stored under Halley's account on CU-Anschutz FTP server, so she is responsible for transferring the new data to Fiji.
**NOTE: If the sequencing is not paid via DTRA grant, make sure to add a different speedtype and PI name in the sample submission form. Also, to avoid delay, add someone from Sawyer Lab in the sample submission form so we are notified as soon as the sequencing is complete.**

Once Halley notifies us that the data has been transferred to Fiji, prepare a copy of the files in your home directory. Don't forget to change the `<path to new sequencing data>` and `<your email>` below:
```
# Set up the environment:
export PATH=~:$PATH

# Create a new directory to store all processing data
mkdir ~/RNA_proc/

# Set the working directory
WD=~/RNA_proc/

# Specify where the new sequencing data is stored
gz=/Shares/dtra_collab/<path to new sequencing data>

# Copy the raw fastq files to your own working directory, this might take a few minutes
sbatch --mem=10g -p short --time=23:50:00 --ntasks=8 --nodes=1 --output=./%x_%j.out --error=./%x_%j.err --mail-type=ALL --job-name=cp \
 --mail-user=<your email> \
 --wrap="cp ${gz}/*.gz ~/RNA_proc/"
```
Once the job is complete, change into your working directory to make sure the files are transferred correctly.

### Step 2: Trim reads, QC, Align & Read Count



### Step 3: Download the read counting results
