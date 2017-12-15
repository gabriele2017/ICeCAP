for FILE in `ls $HM"/"frag_*_"$CHROM"_*.in|gawk '{print $NF}' FS="/"|sort -nk2,2 -t "_" `; do
NNR=$((NNR+1))
gawk '{if(FNR==1){n=split(FILENAME,a,"_")};if(NR>1){printf("%d\t%d\t%d\t%d\t%f\n",a[n-2],'$NNR'*'$RES',a[n-2],NR*'$RES',$5)}}' $HM"/"$FILE |gawk '{if($1==""'$CHROM'){print $0}}' >> $GRAPHICS"/PromoterHeatMapRAW_"$CHROM
echo "" >> $GRAPHICS"/PromoterHeatMapRAW_"$CHROM
gawk '{if(FNR==1){n=split(FILENAME,a,"_");if(NR>1){printf("\n")};printf("chr%s\t%d\t%d\t",a[n-2],a[n-3]*'$RES',(a[n-3]+1)*'$RES');};printf("%f\t",$5)}' $HM"/"$FILE |gawk '{if($1=="chr"'$CHROM'){print $0}}' >> $GRAPHICS"/RawPromoterFlocks_"$CHROM

done
CB=`gawk '{if($5!=0){c=c+1;a=a+$5;print a/c}}' $GRAPHICS"/PromoterHeatMapRAW_"$CHROM|tail -n 1`

/usr/bin/gnuplot <<EOF
set xlabel "Chrom$CHROM"
set ylabel "Chrom$CHROM"
set pm3d map
set size 0.8,1
set palette model CMY rgbformulae 7,5,15
set cbrange [0:$CB]
set terminal postscript enhanced color
set output '$GRAPHICS/PromoterHeatMapRAW_${CHROM}.ps'
splot '$GRAPHICS/PromoterHeatMapRAW_${CHROM}' u 2:4:5
save '$GRAPHICS/PromoterHeatMapRAW_${CHROM}.gnu'
EOF

ps2epsi $GRAPHICS/PromoterHeatMapRAW_${CHROM}.ps  $GRAPHICS/PromoterHeatMapRAW_${CHROM}.epsi
epstopdf $GRAPHICS/PromoterHeatMapRAW_${CHROM}.epsi 

for FILE in `ls $HM"/"frag_*_"$CHROM"_*.fin|gawk '{print $NF}' FS="/"|sort -nk2,2 -t "_" `; do

gawk '{if(FNR==1){n=split(FILENAME,a,"_")};if(NR>1){printf("%d\t%d\t%d\t%d\t%f\n",a[n-2],a[n-3]*'$RES',a[n-2],NR*'$RES',$5)}}' $HM"/"$FILE |gawk '{if($1==""'$CHROM'){print $0}}' >> $GRAPHICS"/PromoterHeatMap_"$CHROM
echo "" >> $GRAPHICS"/PromoterHeatMap_"$CHROM
gawk '{if(FNR==1){n=split(FILENAME,a,"_");if(NR>1){printf("\n")};printf("chr%s\t%d\t%d\t",a[n-2],a[n-3]*'$RES',(a[n-3]+1)*'$RES');};printf("%f\t",$5)}' $HM"/"$FILE |gawk '{if($1=="chr"'$CHROM'){print $0}}' >> $GRAPHICS"/PromoterFlocks_"$CHROM

done

# Bing Ren's Lab Domain caller path required here, including .fai files : default is the home folder.
################################################
~/domaincall_software/perl_scripts/DI_from_matrix.pl $GRAPHICS"/PromoterFlocks_"$CHROM  $RES 50 ~/domaincall_software/chr$CHROM.fa.fai > ~/DI
~/domaincall_software/perl_scripts/DI_from_matrix.pl $GRAPHICS"/PromoterFlocks_"$CHROM  $RES 50 ~/domaincall_software/chr$CHROM.fa.fai >> $GRAPHICS"/Directionality_Index"

octave -qf ~/domaincall_software/HMM_calls.m
gawk '{if(NF==5){print $0}}' ~/DI.output > ~/DI2.output
~/domaincall_software/perl_scripts/file_ends_cleaner.pl ~/DI2.output ~/DI|perl ~/domaincall_software/perl_scripts/converter_7col.pl > ~/DI_$CHROM.output
~/domaincall_software/perl_scripts/hmm_probablity_correcter.pl ~/DI_$CHROM.output 2 0.99 $RES|perl ~/domaincall_software/perl_scripts/hmm-state_caller.pl ~/domaincall_software/chr$CHROM.fai $CHROM|perl ~/domaincall_software/perl_scripts/hmm-state_domains.pl > $GRAPHICS"/temp"

gawk '{print $2 "\t" $2 "\t" 200;print $2 "\t" $3 "\t" 200;print $3 "\t"$3 "\t" 200}' $GRAPHICS"/temp" |gawk '{if(NF==3){print $0}}' > $GRAPHICS"/PromoterFlocks_CHR"$CHROM
gawk '{print $2 "\t" $2 "\t" 200;print $3 "\t" $2 "\t" 200;print $3 "\t"$3 "\t" 200}' $GRAPHICS"/temp" |gawk '{if(NF==3){print $0}}'  >> $GRAPHICS"/PromoterFlocks_CHR"$CHROM

CB=`gawk '{if($5!=0){c=c+1;a=a+$5;print a/c}}' $GRAPHICS"/PromoterHeatMap_"$CHROM|tail -n 1`

/usr/bin/gnuplot <<EOF
set xlabel "Chrom$CHROM"
set ylabel "Chrom$CHROM"
set pm3d map
set size 0.8,1
set palette model CMY rgbformulae 7,5,15
set cbrange [0:$CB]
set terminal postscript enhanced color
set out '$GRAPHICS/PromoterHeatMap_${CHROM}.ps'
splot  '$GRAPHICS/PromoterHeatMap_${CHROM}' u 2:4:5
rep '$GRAPHICS/PromoterFlocks_CHR${CHROM}' w l
save '$GRAPHICS/PromoterFlocks_${CHROM}.gnu'
EOF

ps2epsi $GRAPHICS/PromoterHeatMap_${CHROM}.ps  $GRAPHICS/PromoterHeatMap_${CHROM}.epsi
epstopdf $GRAPHICS/PromoterHeatMap_${CHROM}.epsi 

for FILE in `ls $HM"/"frag_*_"$CHROM"_*.in|gawk '{print $NF}' FS="/"|sort -nk2,2 -t "_" `; do
A=`gawk '{if(FNR==1){n=split(FILENAME,a,"_");printf("%d",a[n-3])}}' $HM"/"$FILE`
B=`gawk '{if(FNR==1){n=split(FILENAME,a,"_");printf("%d",a[n-1])}}' $HM"/"$FILE`
C[$A]=$B
done




while read -a line; do echo -e "${C[${line[1]}]} \t ${C[${line[1]}]} \t 200 \n ${C[${line[1]}]} \t ${C[${line[2]}]} \t 200 \n ${C[${line[2]}]} \t ${C[${line[2]}]} \t 200" >> $GRAPHICS"/HMFlocksposition_CHR"$CHROM;
 done < $GRAPHICS"/temp"
while read -a line; do
echo -e "${C[${line[1]}]} \t ${C[${line[1]}]} \t 200 \n ${C[${line[2]}]} \t ${C[${line[1]}]} \t 200 \n ${C[${line[2]}]} \t ${C[${line[2]}]} \t 200" >> $GRAPHICS"/HMFlocksposition_CHR"$CHROM;
 done < $GRAPHICS"/temp"


echo "Gene Flocks have been populated. Good bye.\n";
