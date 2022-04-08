#TITLE
## setting the run parameters ##
################################

########################
##########################################################################################################

exec_cmd()
{
    echo $*
    if [ -z "$DRY_RUN" ]; then
        eval "$@" ##|| die 'Error'
    fi
}

LANE="FLOWCELL"
REF="REFERENCE"
TRIM=TRIMLENGTH
ZOM="RES"
# scratch folder for temporary files
if [ -n "$WORK_PATH" ]; then
    work_path=$WORK_PATH"/"
else
    work_path=$home_path$LANE/BOWTIE/$TRIM/
fi

if [ -n "$HOME_PATH" ]; then
    ref_path=$HOME_PATH"/"
else
    ref_path=REFBOWTIE
fi

if [ -n "$OUT_PATH" ]; then
   out_path=$OUT_PATH"/"
else
   out_path=$home_path$LANE/data/
#  echo "No output folder specified; quitting."
#    exit
fi

# reference source folder for source etc
if [ -n "$HOME_PATH" ]; then
    ref_pathC=$HOME_PATH"/"
else
    ref_pathC=REFERENCEFOLDEC
fi

gawk '{print $18+$17}' $out_path"/"*pairs*bfide.hic > $out_path"/"GSIZE
gawk '{print $18+$17}' $out_path"/"*pairs*pruned.hic > $out_path"/"BSIZE

file=$LANE

valid=`cat $out_path"/"*bfide.hic|wc -l`

invalid=`cat $out_path"/"*pruned.hic|wc -l`

echo -e "Valid Pairs \t $valid" >> $out_path$file.COUNTS.txt
echo -e "Invalid Pairs \t $invalid" >> $out_path$file.COUNTS.txt

self=`cat $out_path"/"*pairs*pruned.hic|gawk '{if($21==1){count++;print count}}'|wc -l`
echo -e "Invalid Pair: Self-ligations\t$self" >> $out_path$file.COUNTS.txt

same_internal=`cat $out_path"/"*pairs*pruned.hic|gawk '{if($21==2){count++;print count}}'|wc -l`
echo -e "Invalid Pair: Same_internal\t$same_internal" >> $out_path$file.COUNTS.txt

artefacts=`cat $out_path"/"*pairs*pruned.hic|gawk '{if($21==3){count++;print count}}'|wc -l`
echo -e "Invalid Pair: Artefacts\t$artefacts" >> $out_path$file.COUNTS.txt

next_frag=`cat $out_path"/"*pairs*pruned.hic|gawk '{if($21==4){count++;print count}}'|wc -l`
echo -e "Invalid Pair: Next-fragment\t$next_frag" >> $out_path$file.COUNTS.txt

religation=`cat $out_path"/"*pairs*pruned.hic|gawk '{if($21==5){count++;print count}}'|wc -l`
echo -e "Invalid Pair: Religation\t$religation" >> $out_path$file.COUNTS.txt

next_frag2=`cat $out_path"/"*pairs*pruned.hic|gawk '{if($21==6){count++;print count}}'|wc -l`
echo -e "Invalid Pair: Next-fragment\t$next_frag2" >> $out_path$file.COUNTS.txt

Total_pairs_processed=$((valid+invalid))

Self_Ligations=$((invalid-same_internal-artefacts-religation-next_frag-next_frag2))

Same_internal=$((invalid-self-artefacts-religation-next_frag-next_frag2))

Artefacts=$((invalid-self-same_internal-religation-next_frag-next_frag2))

Next_fragment=$((invalid-self-same_internal-religation-artefacts-next_frag2))

Religation=$((invalid-self-same_internal-next_frag-artefacts-next_frag2))

Next_fragment2=$((invalid-self-same_internal-religation-artefacts-next_frag))

frac=`gawk 'BEGIN{print '"$valid"'/'"$Total_pairs_processed"*100'}'`
frac1=`gawk 'BEGIN{print '"$Self_Ligations"'/'"$Total_pairs_processed"*100'}'`
frac2=`gawk 'BEGIN{print '"$Same_internal"'/'"$Total_pairs_processed"*100'}'`
frac3=`gawk 'BEGIN{print '"$Artefacts"'/'"$Total_pairs_processed"*100'}'`
frac4=`gawk 'BEGIN{print '"$Next_fragment"'/'"$Total_pairs_processed"*100'}'`
frac5=`gawk 'BEGIN{print '"$Religation"'/'"$Total_pairs_processed"*100'}'`
frac6=`gawk 'BEGIN{print '"$Next_fragment2"'/'"$Total_pairs_processed"*100'}'`

echo -e "Total_pairs_processed\t$Total_pairs_processed\t100" > $out_path/$LANE.vpairs.pairstat
echo -e "Valid\t$valid\t$frac"  >> $out_path/$LANE.vpairs.pairstat
echo -e "Self_Ligations\t$Self_Ligations\t$frac1"  >> $out_path/$LANE.vpairs.pairstat
echo -e "Same_internal\t$Same_internal\t$frac2" >> $out_path/$LANE.vpairs.pairstat
echo -e "Artefacts\t$Artefacts\t$frac3"  >> $out_path/$LANE.vpairs.pairstat
echo -e "Next_fragment\t$Next_fragment\t$frac4"  >> $out_path/$LANE.vpairs.pairstat
echo -e "Religation\t$Religation\t$frac5"  >> $out_path/$LANE.vpairs.pairstat
echo -e "Next_fragment2\t$Next_fragment2\t$frac6"  >> $out_path/$LANE.vpairs.pairstat
echo -e "Reported_pairs\t$valid\t$frac"  >> $out_path/$LANE.vpairs.pairstat

cmd="R CMD BATCH --no-save --no-restore \"--args picDir='$out_path' bwtDir='$out_path/' sampleName='$LANE' rmMulti='1' rmSingle='1'\" $ref_pathC/plot_pairing_portion2.R $out_path/plot_pairing_portion.Rout"

exec_cmd " $cmd"

cat $out_path"/"*pairs*bfide.hic |grep -v random|gawk '{if(($2==$6&&$3<=$7)||$2<$6){print $1 " " "chr"$2 " " $3 " " $4 " " $5 " " "chr"$6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13 " " $14 " " $15 " " $16}else{print $5 " " "chr"$6 " " $7 " " $8 " " $1 " " "chr"$2 " " $3 " " $4 " "  $12 " " $13 " " $14 " "  $9 " " $10 " " $11 " " $16 " " $15}}' |sort -T $work_path -V -k2,2 -k6,6 -k3,3d -k7,7d > $out_path"/"$file.merged_nodups_pre_format.txt

cat $out_path"/"*pairs*ontarget.hic |grep -v random|gawk '{if(($2==$6&&$3<=$7)||$2<$6){print $1 " " "chr"$2 " " $3 " " $4 " " $5 " " "chr"$6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13 " " $14 " " $15 " " $16}else{print $5 " " "chr"$6 " " $7 " " $8 " " $1 " " "chr"$2 " " $3 " " $4 " "  $12 " " $13 " " $14 " "  $9 " " $10 " " $11 " " $16 " " $15}}' |sort -T $work_path -V -k2,2 -k6,6 -k3,3d -k7,7d > $out_path"/"$file.merged_nodups_pre_format_ontarget.txt

#juicer_tools pre $out_path"/"$file.merged_nodups_pre_format.txt  $out_path"/"$file.merged_nodups_pre_format.hic $ref_path"/"$REF".chrom.sizes"
#juicer_tools pre $out_path"/"$file.merged_nodups_pre_format.txt  $out_path"/"$file.merged_nodups_pre_format_ontarget.hic $ref_path"/"$REF".chrom.sizes"

#cooler makebins /well/jknight/gabriele/ICeCAP/reference_bowtie/hg19.chrom2.sizes  20000 > bins.bed
#cooler load -f coo --one-based ./bins.bed  /well/jknight/gabriele/ICeCAP/HIC_CD14_R1/analysis20/PKY/HIC_CD14_R1_Z20b.pkY ./out.20000
#cooler zoomify out.20000

#clodius aggregate bedpe --chromsizes-filename hg19.chrom3.sizes --chr1-col 1 --chr2-col 1 --from1-col 2 --to1-col 3 --from2-col 2 --to2-col 3 --no-header TADS_CD14 

#higlass-manage ingest --filetype bed2ddb --datatype 2d-rectangle-domains  TADS_CD14.multires.db


for FILE in  $out_path"/"*GSIZE ; do
echo $FILE
/usr/bin/gnuplot <<EOF

set xlabel "d1+d2"
set ylabel "Frequency"
n=1000 #number of intervals
max=3000. #max value
min=0. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5 #fillstyle
set tics out nomirror
set xlabel "x"
set ylabel "Frequency"
#count and plot
set terminal png
set output "${out_path}/${file}_ditag_distribution.png"
plot "${FILE}" u (hist(\$1,width)):(1.0) smooth freq w boxes lc rgb"green" notitle

EOF
done

