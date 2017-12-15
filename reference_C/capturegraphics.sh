
for IN in $(seq $INT); do

for FILE in `ls $HM"/"frag_*_"$CHROM"_"$IN"_*.in|gawk '{print $NF}' FS="/"|sort -nk2,2 -t "_" `; do

cat $HM"/"$FILE >> $GRAPHICS"/HMRAW_"$CHROM"_"$IN
echo "" >> $GRAPHICS"/HMRAW_"$CHROM"_"$IN

gawk '{if(FNR==1){n=split(FILENAME,a,"_");if(NR>1){printf("\n")};printf("chr%s\t%s\t%s\t",a[n-3],a[n-1],a[n-1]+'$RES');};printf("%f\t",$5)}' $HM"/"$FILE >> $GRAPHICS"/RAWTADS_"$CHROM"_"$IN
done

CB=`gawk '{if($5!=0){c=c+1;a=a+$5;print a/c}}' $GRAPHICS"/HMRAW_"$CHROM"_"$IN|tail -n 1`
/usr/bin/gnuplot <<EOF

set title "Enriched_Region_"$IN
set xlabel "$CHROM"
set ylabel "$CHROM"
set pm3d map
set size 0.8,1
set palette model CMY rgbformulae 7,5,15
set cbrange [0:$CB]
set terminal postscript enhanced color
set output '$GRAPHICS/HMRAW_${CHROM}_${IN}.ps'
splot '$GRAPHICS/HMRAW_${CHROM}_${IN}' u 2:4:5
save '$GRAPHICS/HMRAW_${CHROM}_${IN}.gnu'
EOF

# ps2epsi ${GRAPHICS}/HMRAW_${CHROM}.ps ${GRAPHICS}/HMRAW_${CHROM}.epsi
# epstopdf ${GRAPHICS}/HMRAW_${CHROM}.epsi

for FILE in `ls $HM"/"frag_*_"$CHROM"_"$IN"_*.fin|gawk '{print $NF}' FS="/"|sort -nk2,2 -t "_" `; do

cat $HM"/"$FILE >> $GRAPHICS"/HM_"$CHROM"_"$IN
echo "" >> $GRAPHICS"/HM_"$CHROM"_"$IN

gawk '{if(FNR==1){n=split(FILENAME,a,"_");if(NR>1){printf("\n")};printf("chr%s\t%s\t%s\t",a[n-3],a[n-1],a[n-1]+'$RES');};printf("%f\t",$5)}' $HM"/"$FILE  >> $GRAPHICS"/TADS_"$CHROM"_"$IN
echo "" >> $GRAPHICS"/TADS_"$CHROM"_"$IN
done

# Bing Ren's Lab Domain caller path required here, including .fai files : default is the home folder.
################################################
#~/domaincall_software/perl_scripts/DI_from_matrix.pl $GRAPHICS"/TADS_"$CHROM"_"$IN  $RES 2000000 ~/domaincall_software/chr$CHROM.fa.fai > ~/DI
#octave -qf ~/domaincall_software/HMM_calls.m
#gawk '{if(NF==5){print $0}}' ~/DI.output > ~/DI2.output
#~/domaincall_software/perl_scripts/file_ends_cleaner.pl ~/DI2.output ~/DI|perl ~/domaincall_software/perl_scripts/converter_7col.pl > ~/DI_$CHROM_$IN.output
#~/domaincall_software/perl_scripts/hmm_probablity_correcter.pl ~/DI_$CHROM_$IN.output 2 0.99 $RES|perl ~/domaincall_software/perl_scripts/hmm-state_caller.pl ~/domaincall_software/chr$CHROM.fai $CHROM|perl ~/domaincall_software/perl_scripts/hmm-state_domains.pl > $GRAPHICS"/temp"

#gawk '{print $2 "\t" $2 "\t" 200;print $2 "\t" $3 "\t" 200;print $3 "\t"$3 "\t" 200}' $GRAPHICS"/temp" |gawk '{if(NF==3){print $0}}' >>$GRAPHICS"/TADS_CHR"$CHROM"_"$IN
#gawk '{print $2 "\t" $2 "\t" 200;print $3 "\t" $2 "\t" 200;print $3 "\t"$3 "\t" 200}' $GRAPHICS"/temp" |gawk '{if(NF==3){print $0}}' >>$GRAPHICS"/TADS_CHR"$CHROM"_"$IN


#/usr/bin/gnuplot <<EOF
#CB=`gawk '{if($5!=0){c=c+1;a=a+$5;print a/c}}' $GRAPHICS"/HM_"$CHROM"_"$IN|tail -n 1`
#set title "Enriched Region $IN"
#set xlabel "Chrom$CHROM"
#set ylabel "Chrom$CHROM"
#set pm3d map
#set size 0.8,1
#set palette model CMY rgbformulae 7,5,15
#set cbrange [0:$CB]
#set terminal postscript enhanced color
#set output '$GRAPHICS/HM_${CHROM}_${IN}.ps'
#splot '$GRAPHICS/HM_${CHROM}_${IN}' u 2:4:5
#rep '$GRAPHICS/TADS_CHR${CHROM}_${IN}' w l
#save '$GRAPHICS/HM_${CHROM}_${IN}.gnu'
#EOF
#done
# ps2epsi ${GRAPHICS}/HM_${CHROM}.ps ${GRAPHICS}/HM_${CHROM}.epsi
# epstopdf ${GRAPHICS}/HM_${CHROM}.epsi


