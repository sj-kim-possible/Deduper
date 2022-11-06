#!/bin/bash

#SBATCH --partition=bgmp            ### Partition
#SBATCH --job-name=dedup%j          ### Job Name
#SBATCH --output=dedup%j.out        ### File in which to store job output
#SBATCH --error=dedup%j.err         ### File in which to store job error messages
#SBATCH --nodes=1                   ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1         ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1           ### Number of cpus per node
#SBATCH --account=bgmp              ### Account used for job submission

deduper="/projects/bgmp/skim6/bioinfo/Bi624/Deduper-sj-kim-possible/kim_deduper.py"
sorted_sam_input="/projects/bgmp/skim6/bioinfo/Bi624/Deduper-sj-kim-possible/C1_SE_uniqAlign.sorted.sam"
umis="/projects/bgmp/skim6/bioinfo/Bi624/Deduper-sj-kim-possible/STL96.txt"
output_filename="/projects/bgmp/skim6/bioinfo/Bi624/Deduper-sj-kim-possible/C1_SE_uniqAligned_deduped.sorted.sam"

/usr/bin/time -v $deduper -f $sorted_sam_input -u $umis -o $output_filename