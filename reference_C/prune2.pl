#!/usr/bin/perl
use File::Basename;
my $file=$ARGV[0];
my $file2=$ARGV[1];
my $file6=$ARGV[2];
my $file9=$ARGV[3];
my $ints=$ARGV[4];
my $ints2=$ARGV[5];

for my $i($ints..$ints){

$file1=$file.$i;

open FILEIN, $file1||die("Could not open samplelist\n");
open FILEIN2, $file2||die("Could not open samplelist\n");
($file3,$dir,$suffix) = fileparse($file2);
$file3b=~s/chr_$i//;
$file4=$dir.$file3.".pruned";
$file5=$dir.$file3.".bfide";
$file7=$dir.$file3b."badfrags";
$file8=$dir.$file3b."goodfrags";
$file10=$dir.$file3.".bfide_ontarget";

print "open: $file1\n";
print "open: $file2\n";
print "write to: $file5\n";
open SAVE,">> $file5";
open SAVE2,">> $file4";
open SAVE3,">> $file7";
open SAVE4,">> $file8";
open SAVE5,">> $file10";

my $rev1=0;
my $rev2=0;
my $nlines = 0;
my $hd =0;
my $nchrom = 0;
my %list =();
my %remove =();
my @restr_fr=();
my @count_restr_fr=();
my @baits=();
my @count_baits=();
my $d1=0;
my $d2=0;
my $d1prime=0;
my $d2prime=0;

my @line=();
my $count=0;
my $chrom=1;

while(<FILEIN>)
{
print $_;
 $count=$count+1;
       chomp($_); # gets rid of newlines
        $_=~s/^\s+//; # gets rid of leading whitespace
        $_=~s/\s+$//; # gets rid of trailing whitespace
        @line = split(/\t/,$_);
        $line[0]=~s/chr//;
 $list{$count}=[@line];
}

my @name=();
my $counter=0;
my $cell1=0;
my $cell2=0;
my $cell3=0;
my $cell4=0;

my $counter_h=0;

close FILEIN;
$file6b=$file6;

print "open: $file6b\n";
open FILEIN3, $file6b||die("Could not open samplelist\n");
my $max=0;

while(<FILEIN3>)
{
   chomp($_); # gets rid of newlines
        $_=~s/^\s+//; # gets rid of leading whitespace
        $_=~s/\s+$//; # gets rid of trailing whitespace
   @line = split(/\t/,$_);

$line[0]=~s/chr//;
$count_restr_fr[$line[0]]+=1;
$restr_fr[$line[0]][$count_restr_fr[$line[0]]]=$line[3];
}

close FILEIN3;

open FILEIN4, $file9||die("Could not open samplelist\n");
print "open: $file9\n";

while(<FILEIN4>)
{
chomp($_); # gets rid of newlines
        $_=~s/^\s+//; # gets rid of leading whitespace
        $_=~s/\s+$//; # gets rid of trailing whitespace
@line = split(/\t/,$_);
$line[0]=~s/chr//;
if($line[0] eq "X"){$line[0]=$ints2-2;}
if($line[0] eq "Y"){$line[0]=$ints2-1;}
if($line[0] eq "M"){$line[0]=$ints2;}
$count_baits[$line[0]]+=1;
$baits[$line[0]][$count_baits[$line[0]]]=int(($line[1]+$line[2])/2);
}
close FILEIN4;

my $save_a=0;
my $save_b=0;
my $save_aprime=0;
my $save_bprime=0;

my $dist1=0.;
my $dist2=0.;

open FILEIN2, $file2||die("Could not open samplelist\n");
print "open: $file2\n";

while(<FILEIN2>)
{
$counter=$counter+1;

chomp($_); # gets rid of newlines
$_=~s/^\s+//; # gets rid of leading whitespace
$_=~s/\s+$//; # gets rid of trailing whitespace

@name=split(/\t/,$_);
$name[2]=~s/chr//;
$name[6]=~s/chr//;


for my $i(1..$count){
my @listed = @{$list{$i}};

$cell1=grep(/=/,$name[6]);
$cell2=grep(/X/,$name[6]);
$cell3=grep(/Y/,$name[6]);
$cell4=grep(/M/,$name[6]);

$cell6=grep(/X/,$name[2]);
$cell7=grep(/Y/,$name[2]);
$cell8=grep(/M/,$name[2]);

if($cell2==1){$name[6]=$ints2-2;}
if($cell3==1){$name[6]=$ints2-1;}
if($cell4==1){$name[6]=$ints2;}
if($cell6==1){$name[2]=$ints2-2;}
if($cell7==1){$name[2]=$ints2-1;}
if($cell8==1){$name[2]=$ints2;}
if($cell1==1){$name[6]=$name[2];}

if($name[2]==$listed[0]){
$tmp=sprintf("%i",$name[6]);

if($tmp<25&&$tmp>0&&($name[3]>=$listed[1]&&$name[3]<=$listed[2])){
    if($name[1]==1||$name[1]==17||$name[1]==33||$name[1]==49){
    if($name[1]==1){$rev1=1;$rev2=1;}
    if($name[1]==17){$rev1=0;$rev2=1;}
    if($name[1]==33){$rev1=1;$rev2=0;}
    if($name[1]==49){$rev1=0;$rev2=0;}

my $flag=0;


my $up=$count_restr_fr[$name[2]];
my $down=1;
my $flag3=0;
my $middle=int(($down+$up)/2);

while ($flag3<1){
if($name[3]>$restr_fr[$name[2]][$middle]){$down=$middle; }
if($name[3]<$restr_fr[$name[2]][$middle]){$up=$middle; }
if($name[3]==$restr_fr[$name[2]][$middle]){$up=$middle; $flag3=1; $save_a=$middle; }
$middle=int(($up+$down)/2);
if($up==($down+1)){$flag3=2; $save_a=$down;}
#print "enzyme: $down $up $middle $name[3] $restr_fr[$name[2]][$save_a+$rev1] \n";
}

if($up==($down+2)){$flag3=1; $save_a=$middle;}
$d1=abs($name[3]-$restr_fr[$name[2]][$save_a+$rev1]);
$d1prime=0;
if($count_baits[$name[2]]>0){
my $up=$count_baits[$name[2]];
my $down=1;
my $flag3=0;
my $middle=int(($down+$up)/2);

#print "Bait: $name[2] $count_baits[$name[2]] \n";
while ($flag3<1){
if($name[3]>$baits[$name[2]][$middle]){$down=$middle; }
if($name[3]<$baits[$name[2]][$middle]){$up=$middle; }
if($name[3]==$baits[$name[2]][$middle]){$up=$middle; $flag3=1; $save_aprime=$middle; }
$middle=int(($up+$down)/2);
if($up==($down+1)){$flag3=2; $save_aprime=$down;}
#print "Bait: $down $up $middle $name[3] $baits[$name[2]][$middle+$rev1] \n";
}

if($up==($down+2)){$flag3=1; $save_aprime=$middle;}
$d1prime=abs($name[3]-$baits[$name[2]][$save_aprime+$rev1]);
}
my $up=$count_restr_fr[$name[6]];
my $down=1;
my $flag3=0;
my $middle=int(($down+$up)/2);

#print "$name[0] \t $name[1] \t  $name[2] \t $name[3] \t $name[4] \t $name[5] \t $name[6] \t $name[7]\t $name[8] \n" ;
while ($flag3<1){
if($name[7]>$restr_fr[$name[6]][$middle]){$down=$middle; }
if($name[7]<$restr_fr[$name[6]][$middle]){$up=$middle; }
if($name[7]==$restr_fr[$name[6]][$middle]){$up=$middle; $flag3=1; $save_b=$middle; }
$middle=int(($up+$down)/2);
if($up==($down+1)){$flag3=2; $save_b=$down;}
#print "enzyme: $down $up $middle $name[7] $restr_fr[$name[6]][$save_b+$rev2] \n";
}
if($up==($down+2)){$flag3=1; $save_b=$middle;}
$d2=abs($name[7]-$restr_fr[$name[6]][$save_b+$rev2]);

$d2prime=0;
if($count_baits[$name[6]]>0){
my $up=$count_baits[$name[6]];
my $down=1;
my $flag3=0;
my $middle=int(($down+$up)/2);
while ($flag3<1){
if($name[7]>$baits[$name[6]][$middle]){$down=$middle; }
if($name[7]<$baits[$name[6]][$middle]){$up=$middle; }
if($name[7]==$baits[$name[6]][$middle]){$up=$middle; $flag3=1; $save_bprime=$middle; }
$middle=int(($up+$down)/2);
if($up==($down+1)){$flag3=2; $save_bprime=$down;}
#print "Baits: * $down * $up * $middle * $name[6] * $baits[$name[6]][$save_bprime+$rev2] \n";
}
if($up==($down+2)){$flag3=1; $save_bprime=$middle;}
$d2prime=abs($name[7]-$baits[$name[6]][$save_bprime+$rev2]);
}
if(($save_a==$save_b&&$name[2]==$name[6])){
if($name[1]==17&&$name[3]<$name[7]){$flag=1;}
if($name[1]==17&&$name[3]>=$name[7]){$flag=2;}
if($name[1]==33&&$name[3]<$name[7]){$flag=2;}
if($name[1]==33&&$name[3]>=$name[7]){$flag=1;}

if($name[1]==1||$name[1]==49){$flag=3;}
}

if(($save_a==($save_b+1)&&$name[2]==$name[6])||($save_b==$save_a+1&&$name[2]==$name[6])){
if($name[1]==17&&$name[3]<$name[7]){$flag=4;}
if($name[1]==17&&$name[3]>=$name[7]){$flag=5;}
if($name[1]==33&&$name[3]<$name[7]){$flag=5;}
if($name[1]==33&&$name[3]>=$name[7]){$flag=4;}

if($name[1]==1||$name[1]==49){$flag=6;}
}

my $size=$d1+$d2;
if($flag>0){
printf (SAVE2 "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t%f\t%f\n",$name[0],$name[1],$name[2],$name[3],$name[4],$name[5],$name[6],$name[7],$name[8],$name[9],$name[10],$d1,$d2,$flag,$d1prime,$d2prime);

printf (SAVE3 "%f\n",$d1+$d2);
}else{
printf (SAVE "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t%f\t%f\t%d\t%d\n",$name[0],$name[1],$name[2],$name[3],$name[4],$name[5],$name[6],$name[7],$name[8],$name[9],$name[10],$d1,$d2,0,$d1prime,$d2prime,$save_aprime,$save_bprime);
if(($d1>$d1prime||$d2>$d2prime)&&($d1prime>0&&$d2prime>0)){
printf (SAVE5 "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t%f\t%f\n",$name[0],$name[1],$name[2],$name[3],$name[4],$name[5],$name[6],$name[7],$name[8],$name[9],$name[10],$d1,$d2,0,$d1prime,$d2prime);}
printf (SAVE4 "%f\n",$d1+$d2);
}
}else{
$flag=7;
printf (SAVE2 "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t%f\t%f\n",$name[0],$name[1],$name[2],$name[3],$name[4],$name[5],$name[6],$name[7],$name[8],$name[9],$name[10],$d1,$d2,$flag,$d1prime,$d2prime);}
}
}
}
}

print "$file5 is really finished now! \n";
close FILEIN2;
close SAVE;
close SAVE2;
}

