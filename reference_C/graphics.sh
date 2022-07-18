
for FILE in `ls $HM/*tad `; do

echo $FILE

cat $FILE |gawk '{if($1=="chr"'$CHROM'){print $0}}' > $GRAPHICS"/TADS_"$CHROM

done


# Bing Ren's Lab Domain caller path required here, including .fai files : default is the home folder.
################################################
~/domaincall_software/perl_scripts/DI_from_matrix.pl $GRAPHICS"/TADS_"$CHROM  $RES 2000000 ~/domaincall_software/chr$CHROM.fa.fai > ~/DI
octave -qf ~/domaincall_software/HMM_calls.m
gawk '{if(NF==5){print $0}}' ~/DI.output > ~/DI2.output
~/domaincall_software/perl_scripts/file_ends_cleaner.pl ~/DI2.output ~/DI|perl ~/domaincall_software/perl_scripts/converter_7col.pl > ~/DI_$CHROM.output
~/domaincall_software/perl_scripts/hmm_probablity_correcter.pl ~/DI_$CHROM.output 2 0.99 $RES|perl ~/domaincall_software/perl_scripts/hmm-state_caller.pl ~/domaincall_software/chr$CHROM.fai $CHROM|perl ~/domaincall_software/perl_scripts/hmm-state_domains.pl > $GRAPHICS"/temp"


gawk '{print $1 "\t" $2 "\t" $3 "\t" "TAD_"NR}' $GRAPHICS"/temp" > $GRAPHICS"/TADS_CHR"$CHROM
gawk '{print $1 "\t" $2 "\t" $3 "\t" "TAD_"NR}' $GRAPHICS"/temp" >> $GRAPHICS"/TADS.bed"


