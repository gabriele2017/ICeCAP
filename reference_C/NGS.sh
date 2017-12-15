#TITLE
## setting the run parameters ##
################################

########################
# zcat
if [ -n "$ZCAT_PATH" ]; then
    zcat=$ZCAT_PATH"/"
else
    zcat=/bin/zcat
fi

# gawk
if [ -n "$GAWK_PATH" ]; then
    gawk=$GAWK_PATH"/"
else
    gawk=/usr/bin/gawk
fi

# split
if [ -n "$SPLIT_PATH" ]; then
    split=$SPLIT_PATH"/"
else
    split=/usr/bin/split
fi

## input and output folders and file prefixes
#############################################

# reference folder for indexes, fragment lists etc
if [ -n "$HOME_PATH" ]; then
    ref_path=$HOME_PATH"/"
else
    ref_path=REFBOWTIE
fi

if [ -n "$HOME_PATH" ]; then
    ref_path=$HOME_PATH"/"
else
ref_path2=REFERENCEFOLDER
fi

# reference source folder for source etc
if [ -n "$HOME_PATH" ]; then
    ref_pathC=$HOME_PATH"/"
else
    ref_pathC=REFERENCEFOLDEC
fi

# home folder for tools
if [ -n "$HOME_PATH" ]; then
    home_path=$HOME_PATH"/"
else
    home_path=HOMEFOLDER
fi

ALLSPEC="ALLELE"
LANE="FLOWCELL"
STRING="DITAG"
SAVE="STARTPATHL"
SAVE2="STARTPATHR"
OLIGOS="BAITS"
ENZYME="CUTTER"
REF="REFERENCE"
CORE=NODES
TRIM=TRIMLENGTH
CHREND=CHROMEND
# pruning
if [ -n "$PRUNE" ]; then
    prune=$PRUNE
else
    prune=$ref_pathC"prune2.pl"
fi

#############################################

if [ -n "$GENOME_REF_FASTA" ]; then
    genome_ref_fasta=$GENOME_REF_FASTA
else
    genome_ref_fasta=$ref_path$REF
fi

file=$LANE
trim=$TRIM
mkdir -p $home_path$LANE
mkdir -p $home_path$LANE/BOWTIE/
mkdir -p $home_path$LANE/BOWTIE/$TRIM
mkdir -p $home_path$LANE/data/

if [ -n "$SHIFTREV_PATH" ]; then
shiftrev=$SHIFTREV_PATH"/"
else
shiftrev=$ref_pathC"/"shiftrev
fi

# scratch folder for temporary files
if [ -n "$WORK_PATH" ]; then
    work_path=$WORK_PATH"/"
else
    work_path=$home_path$LANE/BOWTIE/$TRIM/
fi

# output folder

if [ -n "$OUT_PATH" ]; then
   out_path=$OUT_PATH"/"
else
   out_path=$home_path$LANE/data/
#  echo "No output folder specified; quitting."
#    exit
fi

## other run parameters & files
###############################

# input folder
if [ -n "$IN_PATH" ]; then
    in_path=$IN_PATH"/"
else

    in_path=$SAVE
#    echo "No input folder specified; quitting."
#    exit
fi

# input folder
if [ -n "$IN_PATH" ]; then
    in_path2=$IN_PATH"/"
else

    in_path2=$SAVE2
#    echo "No input folder specified; quitting."
#    exit
fi

# input cell_line/replicate name prefix
if [ -n "$NAME" ]; then
    name=$NAME"/"
else
    name=$LANE
#    echo "No input file prefix specified; quitting."
#    exit
fi

# print all input parameters:

echo -e "Input parameters:"
echo -e "\tHOME_PATH = < $home_path >"
echo -e "\tWORK_PATH = < $work_path >"
echo -e "\tOUT_PATH = < $out_path >"
echo -e "\tIN_PATH = < $in_path >"
echo -e "\tIN_PATH = < $in_path2 >"
echo -e "\tGENOME_REF_FASTA = < $genome_ref_fasta >"
echo -e ""

echo -e "Input parameters:" >> $work_path"HiCinC.txt"
echo -e "\tHOME_PATH = < $home_path >"  >> $work_path"HiCinC.txt"
echo -e "\tSCRATCH = < $work_path >"  >> $work_path"HiCinC.txt"
echo -e "\tOUT_PATH = < $out_path >"  >> $work_path"HiCinC.txt"
echo -e "\tIN_PATH = < $in_path >"  >> $work_path"HiCinC.txt"
echo -e "\tIN_PATH = < $in_path2 >"  >> $work_path"HiCinC.txt"
echo -e "\tGENOME_REF_FASTA = < $genome_ref_fasta>"  >> $work_path"HiCinC.txt"
echo -e ""  >> $work_path"HiCinC.txt"

#####################################################

######################
## START PROCESSING ##
######################

echo -e "Starting the data processing."

for path in $in_path $in_path2; do

#####################################################
## STEP 0: .gz to .fastq conversion (if required) ###
#####################################################

################################
## STEP 1: alignment          ##
################################
IFS=,
tmp=($path)
for key in "${!tmp[@]}"; do
FILE=`echo "${tmp[$key]}"`;
NAME=`echo "${tmp[$key]}" | gawk -F'/' '{print $NF}'| sed s/\.fastq//| sed s/\.fq//| sed s/\.gz//`;

echo "$FILE";
echo "..Processing: $NAME";

if [ ${FILE: -3} == ".gz" ]; then
           echo -e "\t....gunzip: < gunzip $FILE >"

gunzip $FILE
FILE=`echo "$FILE" | sed s/\.gz//`;

fi

$split -dl 24000000 $FILE $work_path"/"$NAME"_split_"

echo -e "\n Raw reads: $NAME \t" >> $out_path$LANE.COUNTS.txt
gawk '{if(NR%4==2){print $0}}' $FILE |wc -l|tail -n 1 >> $out_path$LANE.COUNTS.txt

done
done

for filename in $work_path"/"*_split_*;
do
cat $filename | gawk '{if(index($0,"'$STRING'")!=0){flag=index($0,"'$STRING'")-1; count++; len+=length($0)+length("'$ENZYME'");};if((NR%4==2)&&flag>20){print substr($0,1,flag)"'$ENZYME'"}; if((NR%4==0)&&flag>20){  print substr($0,1,flag+length("'$ENZYME'"))};if((NR%4==0||NR%4==2)&&flag==0){print $0}else if((NR%4==0||NR%4==2)&&flag>0&&flag<=20){print substr($0,1,1)}else if(NR%4!=0&&NR%4!=2){print $0};if(NR%4==0){flag=0};}' > $filename".trimmed.fq"
rm $filename
done

for filename in $work_path"/"*trimmed.fq; do
                echo -e "\t..running BOWTIE on: < $filename >";
                file3=`echo $filename | awk -F"/" '{print $NF}'`
                file2=`echo $file3 | sed s/\.fq//g`;
                file=`echo $file2 | sed s/\_split//`;

bowtie2 -p $CORE --very-sensitive -x $genome_ref_fasta -U $filename -S $work_path$file.bowtiesam 2> $work_path"/"$file"_bowtie_out.txt"

rm $filename
done
head -n 10000 $work_path*00.trimmed.bowtiesam | grep  '^[[:space:]]*@' > $work_path$"/"header

for filename in $work_path$"/"*.bowtiesam; do
file2=`echo $filename | awk -F"/" '{print $NF}'`
file=`echo $file2 | sed s/\.bowtiesam//g`;
$shiftrev $filename $work_path$file $ALLSPEC
rm $filename
done

#######################################################
## STEP 2:  (pair-matching and de-duplication)       ##
#######################################################

for filename in $work_path*.sam; do
file2=`echo $filename | awk -F"/" '{print $NF}' | sed s/\.trimmed.sam//g`
file3=`echo $file2 | awk -F"_" '{printf("%s_%s_",$1,$3);}'`
file=$file3"all"
grep -v '^[[:space:]]*@' $filename >> $work_path$file.sam 
done

file=$LANE

for filename in $work_path*all.sam; do
cat $filename | sort -k1,1 -T $work_path >> $work_path$file.sorted.sam
done

rm $work_path*all.sam

echo -e "\n Aligned reads: \t" >> $out_path$file.COUNTS.txt
wc -l $work_path$file.sorted.sam  >> $out_path$file.COUNTS.txt

gawk '{for(i=1;i<=NF;i++){vec[i]=$i;}; if(vec[1]==prev[1]&&(vec[2]+prev[2]+1<200)){if(vec[5]<prev[5]){prev[5]=vec[5];};printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t",vec[1],vec[2]+2*prev[2]+1,vec[3],vec[4],prev[5],vec[6],prev[3],prev[4]);for(i=9;i<=NF;i++){printf("%s\t",vec[i]);}; printf("\n"); printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t",prev[1],prev[2]+2*vec[2]+1,prev[3],prev[4],prev[5],prev[6],vec[3],vec[4]);for(i=9;i<=NF;i++){printf("%s\t",prev[i]);}; printf("\n");}; for(i=1;i<=NF;i++){prev[i]=$i;}; }' FS='\t'  $work_path$file.sorted.sam > $work_path$file.sorted.matefixed.sam

################################################################################################################################################

echo -e "\n Sorted&Matefixed Reads \t" >> $out_path$file.COUNTS.txt
wc -l $work_path$file.sorted.matefixed.sam  >> $out_path$file.COUNTS.txt

rm $work_path$file.sorted.sam

##############################################################################################################################################

cat $work_path$file.sorted.matefixed.sam | gawk '{if($5>29&&$4!=0&&$8!=0){print $0}}'>$work_path$file.sorted.matefixed.MQ30.sam

rm $work_path$file.sorted.matefixed.sam
echo "QC reads (MQ>=30)" >> $out_path$file.COUNTS.txt

wc -l $work_path$file.sorted.matefixed.MQ30.sam >> $out_path$file.COUNTS.txt

gawk '{for(i=1;i<=NF;i++){vec[i]=$i;}; if(vec[1]==prev[1]){for(i=1;i<=NF;i++){printf("%s\t",vec[i]);};printf("\n");for(i=1;i<=NF;i++){printf("%s\t",prev[i]);}; printf("\n");}; for(i=1;i<=NF;i++){prev[i]=$i;}}' FS='\t' $work_path$file.sorted.matefixed.MQ30.sam > $work_path$file.sorted.matefixed.MQ30.paired.sam
rm $work_path$file.sorted.matefixed.MQ30.sam
echo "Re-paired" >> $out_path$file.COUNTS.txt

wc -l $work_path$file.sorted.matefixed.MQ30.paired.sam >> $out_path$file.COUNTS.txt

########################################################################################################

for i in $(seq $((CHREND-3)) ); do

cat $work_path$file.sorted.matefixed.MQ30.paired.sam |gawk '{if($3=="chr'$i'"){b=$4; c=$7; d=$8; print b"_"c"_"d "\t" $0  }}' | sort -k1,1n -T $work_path|gawk '{for(i=2;i<=NF;i++){printf("%s\t",$i);}printf("\n");}' >  $out_path$file.pairs.chr_$i.sam

done

for i in X Y M ; do 
cat $work_path$file.sorted.matefixed.MQ30.paired.sam |gawk '{if($3=="chr'$i'"){b=$4; c=$7; d=$8; print b"_"c"_"d "\t" $0  }}' | sort -k1,1n -T $work_path|gawk '{for(i=2;i<=NF;i++){printf("%s\t",$i);}printf("\n");}' >  $out_path$file.pairs.chr_$i.sam

done

rm $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.sam 
mv $out_path$file.pairs.chr_X.sam $out_path$file.pairs.chr_"$((CHREND-2))".sam
mv $out_path$file.pairs.chr_Y.sam $out_path$file.pairs.chr_"$((CHREND-1))".sam
mv $out_path$file.pairs.chr_M.sam $out_path$file.pairs.chr_"$((CHREND))".sam

##########################################################################################################


# print all input parameters:

# lane identifiers
if [ -n "$CHR_A" ]; then
chr_a=$CHR_A
else
chr_a=1
#    echo "Input lanes not specified; quitting."
#    exit
fi

if [ -n "$CHR_B" ]; then
    lane_b=$CHR_B
else
chr_b=$((CHREND))
fi

for (( i=$chr_a; i <=$chr_b; i=i+1 )); do
for filename in $out_path$"/"*.pairs.chr_$i.sam; do
echo -e "$prune $ref_path2""chr $filename $ref_path2$ENZYME $OLIGOS $i $CHREND" >> $out_path"pruning_chr"$i
$prune $ref_path2"chr" $filename $ref_path2$ENZYME $OLIGOS $i $CHREND
done
wait
done

cat $out_path"/"*.pairs*bfide | sort -k1,1 -T $work_path > $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.sam 
cat $out_path"/"*.pairs*bfide_ontarget | sort -k1,1 -T $work_path > $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.targeted.sam

rm $out_path"/"*.pairs*bfide 
rm $out_path"/"*.pairs*bfide_ontarget

gawk '{for(i=1;i<=NF;i++){vec[i]=$i;}; if(vec[1]==prev[1]){for(i=1;i<=NF;i++){printf("%s\t",vec[i]);};printf("\n");for(i=1;i<=NF;i++){printf("%s\t",prev[i]);}; printf("\n");}; for(i=1;i<=NF;i++){prev[i]=$i;}}' FS='\t' $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.sam > $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.repaired.sam
gawk '{for(i=1;i<=NF;i++){vec[i]=$i;}; if(vec[1]==prev[1]){for(i=1;i<=NF;i++){printf("%s\t",vec[i]);};printf("\n");for(i=1;i<=NF;i++){printf("%s\t",prev[i]);}; printf("\n");}; for(i=1;i<=NF;i++){prev[i]=$i;}}' FS='\t' $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.targeted.sam > $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.targeted.repaired.sam

rm $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.sam
rm $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.targeted.sam

for i in $(seq "$CHREND");
do
cat $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.repaired.sam |gawk '{if($3=="'$i'"){b=$4; print b "\t" $0 }}' | sort -k1,1n -T $work_path |gawk '{for(i=2;i<=NF;i++){printf("%s\t",$i);}printf("\n");}' >  $out_path$file.pairs.chr_$i.sam.bfide
cat $work_path$file.sorted.matefixed.MQ30.dedupl.sorted.paired.pruned.targeted.repaired.sam |gawk '{if($3=="'$i'"){b=$4; print b "\t" $0 }}' | sort -k1,1n -T $work_path |gawk '{for(i=2;i<=NF;i++){printf("%s\t",$i);}printf("\n");}' >  $out_path$file.pairs.chr_$i.sam.bfide_ontarget

echo "Deduplicated Valid Paired reads per chromosome: chr'$i'" >> $out_path$file.COUNTS.txt
wc -l $out_path$file.pairs.chr_$i.sam >> $out_path$file.COUNTS.txt
done

echo -e "Valid Pairs" >> $out_path$file.COUNTS.txt
wc -l $out_path"/"*.pairs*bfide|tail -1 >> $out_path$file.COUNTS.txt

gawk '{print $12+$13}' $out_path"/"*pairs*bfide > $out_path"/"SIZES

echo -e "Invalid Pairs" >> $out_path$file.COUNTS.txt
wc -l $out_path"/"*pruned|tail -1 >> $out_path$file.COUNTS.txt

echo -e "Invalid Pair: Self-ligations" >> $out_path$file.COUNTS.txt
gawk '{if($14==1){count++;print count}}' $out_path"/"*.pairs*pruned|tail -1 >> $out_path$file.COUNTS.txt

echo -e "Invalid Pair: Same Internal/Same Dangling Ends" >> $out_path$file.COUNTS.txt
gawk '{if($14==2){count++;print count}}' $out_path"/"*.pairs*pruned|tail -1 >> $out_path$file.COUNTS.txt

echo -e "Invalid Pair: Artefacts" >> $out_path$file.COUNTS.txt
gawk '{if($14==3){count++;print count}}' $out_path"/"*.pairs*pruned|tail -1 >> $out_path$file.COUNTS.txt

echo -e "Invalid Pair: Next fragment" >> $out_path$file.COUNTS.txt
gawk '{if($14==4){count++;print count}}' $out_path"/"*.pairs*pruned|tail -1 >> $out_path$file.COUNTS.txt

echo -e "Invalid Pair: Re-ligation/uncut" >> $out_path$file.COUNTS.txt
gawk '{if($14==5){count++;print count}}' $out_path"/"*.pairs*pruned|tail -1 >> $out_path$file.COUNTS.txt

echo -e "Invalid Pair: Next Fragment" >> $out_path$file.COUNTS.txt
gawk '{if($14==6){count++;print count}}' $out_path"/"*.pairs*pruned|tail -1 >> $out_path$file.COUNTS.txt

cat $out_path"/"*badfrags >> $out_path"/BSIZE"
cat $out_path"/"*goodfrags >> $out_path"/GSIZE"
rm $out_path"/"*badfrags
rm $out_path"/"*goodfrags

for FILE in  $out_path"/"*SIZE ; do
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

gawk '/Pair\:/{print}' $out_path$file.COUNTS.txt > $out_path"/"names.txt
gawk '/Pair\:/{getline; print}' $out_path$file.COUNTS.txt > $out_path"/"values.txt
paste $out_path"/"names.txt $out_path"/"values.txt > $out_path"/"pietable.txt
rm $out_path"/"names.txt $out_path"/"values.txt
gawk -v OFS=, '{print ($1, $2)}' FS='\t' $out_path"/"pietable.txt > $out_path"/"piedata.csv

gawk '/Aligned reads/{getline; print ("Aligned", $1)}' $out_path$file.COUNTS.txt > $out_path"/"values.txt
gawk '/Sorted/{getline; print ("Sorted",$1)}' $out_path$file.COUNTS.txt >> $out_path"/"values.txt
gawk '/QC reads /{getline; print ("QC",$1)}' $out_path$file.COUNTS.txt >> $out_path"/"values.txt
gawk '{if($1=="Deduplicated"&&$2!="Valid"&&$3!="resorted"){getline; print ("Deduplicated",$1)}}' $out_path$file.COUNTS.txt >> $out_path"/"values.txt
gawk '/Resorted/{getline; print("Resorted", $1)}' $out_path$file.COUNTS.txt >> $out_path"/"values.txt

seq 0 5 > $out_path"/"order.txt
paste $out_path"/"order.txt $out_path"/"values.txt > $out_path"/"values.dat

for FILE3 in  $out_path"/"*.dat ; do
/usr/bin/gnuplot <<EOF

set terminal png
set output '${out_path}/${file}_NGS_histogram.png'
unset key
set boxwidth 0.5
set format y "%.0s*10^%T"
set ylabel 'Number of Reads'
set yrange [0:]
set style fill solid
plot "${FILE3}" using 1:3:xtic(2) with boxes
EOF

done

for FILE2 in  $out_path"/"*csv ; do
/usr/bin/gnuplot <<EOF
set terminal png
set output '${out_path}/${file}_Invalid_ditags_chart.png'
set datafile separator ','
stats '${FILE2}' u 2 noout
ang(x)=x*360.0/STATS_sum
perc(x)=x*100.0/STATS_sum
set size square
set xrange [-1:1.5]
set yrange [-1.25:1.25]
set style fill solid 1
unset border
unset tics
unset key
Ai = 0.0; Bi = 0.0;
mid = 0.0;
i = 0; j = 0;
yi  = 0.0; yi2 = 0.0;
set output '${out_path}/${file}_Invalid_ditags_chart.png'
plot "${FILE2}" u (0):(0):(1):(Ai):(Ai=Ai+ang(\$2)):(i=i+1) with circle linecolor var,"${FILE2}"  u (1.3):(yi2=yi2+0.5/STATS_records):(j=j+1) w p pt 5 ps 2 linecolor var,"${FILE2}" u (mid=Bi+ang(\$2)*pi/360.0, Bi=2.0*mid-Bi, 0.5*cos(mid)):(0.5*sin(mid)):(sprintf('%.0f (%.1f\%)', \$2, perc(\$2))) w labels
EOF
done

