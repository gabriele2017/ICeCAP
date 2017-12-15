#!/bin/bash
#BSUB -W 168:00 # Run max
#BSUB -e ICeCap_II.e # stderr, %J is the job ID
#BSUB -J ICeCap_II   # stdout, %J is the job ID
#BSUB -o ICeCAP_II # stdout, %J is the job ID
#BSUB -R "span[hosts=1]"
#BSUB -n 12

START=1
END=3
STEP=2

## save $START, just in case if we need it later ##
i=$START
while [[ $i -le $END ]]
do
      j=$(($i+$STEP))
/well/jknight/gabriele/ICeCap/ICeCAP -N stage2 -l GMrep1  -B /well/jknight/gabriele/ICeCap/baits.bed  -Q /well/jknight/gabriele/ICeCap/reference_map/ -s $i -e $j -S ziw -G 1 -Z 5 -M 200 -R /well/jknight/gabriele/ICeCap/reference_bowtie/ -f /well/jknight/gabriele/ICeCap/ -I /well/jknight/gabriele/ICeCap/GMrep1/data/tmp/  -O /well/jknight/gabriele/ICeCap/GMrep1/data/tmp/ -X yes &
((i = i + $STEP+1))
done

wait
