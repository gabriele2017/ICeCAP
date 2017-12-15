#!/bin/bash
#$ -cwd
#$ -q long.qc 
#$ -P jknight.prjc
#$ -N GM_rep1
#$ -o /well/jknight/gabriele/ICeCap/to
#$ -e /well/jknight/gabriele/ICeCap/te
#$ -pe shmem 6

/well/jknight/gabriele/ICeCap/ICeCAP -J 12 -E AAGCTT -C 1 -P ./JURKATRAW -R /well/jknight/gabriele/ICeCap/reference_bowtie -l JURKATrep1 -T 100 -B /well/jknight/gabriele/ICeCap/jurkatbaits.bed -f /well/jknight/gabriele/ICeCap/
wait
