#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name="trimming" # Job name
#SBATCH --array=1-64
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=4 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=20000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=01-00:00:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=trimming_%A_%a.out # STDOUT
#SBATCH --error=trimming_%A_%a.err # STDERR
#SBATCH --mail-user=rzlli@ucdavis.edu # does not work yet
#SBATCH --mail-type=ALL # does not work yet
# This will be run once for a single process

/bin/hostname

start=`date +%s`
module load trimmomatic/0.33 

cd /share/malooflab/Fastqs/Chunmei/RNAseq

# trimming 
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list_trimming`

echo ${sample}
outdir=/share/malooflab/Chunmei/trimmed

trimmomatic SE -threads 4 ${sample}_L007_R1_001.fastq.gz $outdir/${sample}_trimmed.fq.gz ILLUMINACLIP:Bradseq_adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

end=`date +%s`
runtime=$((end-start))
echo $runtime seconds to completion
