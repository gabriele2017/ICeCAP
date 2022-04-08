#!/usr/bin/perl

use warnings;
use Getopt::Long;
use POSIX;

if (@ARGV < 2){
    die "need more arguments
  arg 0 = sam output from bwa
  arg 1 = list of chrs and lengths
  arg 2 = window size
  arg 3 = chromosome we're working on
  arg 4 = read length
  arg 5 = ouput directory
";
}

my $windSize = $ARGV[2];

#first hash reference lengths (entrypoints)
my %chrLengths = ();

my $file=$ARGV[1];
open FILEIN, $file or die "Cannot open ($!).";
while(<FILEIN>){

my $line=$_;
    chomp($line);
    my @fields = split("\t",$line);
    $fields[0]=~s/chr//;
    $chrLengths{$fields[0]} = $fields[2];
}

close FILEIN;

my @posArr = ();
my @gcArr = ();
my @nArr = ();

my $chr = $ARGV[3];

 $file=$ARGV[0];
open FILEIN, $file or die "Cannot open ($!).";
while(<FILEIN>){
my $line=$_;
    chomp ($line);
    #skip header lines
    unless ($line =~ /^@/)
    {
        #assumes the file is named with suffix .fa or .fasta
	#if ($line =~ /^([^\.]+)\.fa(sta)?:(\d+)-\d+/)
	#{
	    my @fields = split("\t",$line);
            $fields[2]=~s/chr//;
	    my $start = $fields[3];
            my $thisChr = $fields[2];
            unless ($fields[2] eq "*"){

		#print STDERR "$chr == $fields[2]\n";
		#print STDERR "$start == $fields[3]\n";	       
		#print STDERR "$chr == $thisChr\n";

		if ( ($fields[2] eq $chr) && ($fields[3] == $start) 
		     && $chr eq $thisChr)
		{ 

		    #calculate appropriate window
		    my $key = floor($start/$windSize);
		     #print STDERR $key . "\n";
		    if (!(exists($posArr[$key]))){
			$posArr[$key] = 1;
			$gcArr[$key] = ($fields[9] =~ tr/GCgc//);
			$nArr[$key] = ($fields[9] =~ tr/nN//);
		    } else {
			$posArr[$key]++;
			$gcArr[$key] = $gcArr[$key] + ($fields[9] =~ tr/gcGC//);
			$nArr[$key] = $nArr[$key] +($fields[9] =~ tr/nN//);

		    }
		}	
	    #}	
	}
    }
}

my $outFh = open (my $mapfile, ">$ARGV[5]/$ARGV[3].map") || die "Can't open output file.\n";
my $outFh2 = open (my $gcfile, ">$ARGV[5]/$ARGV[3].gc") || die "Can't open output file.\n";

my $len = $chrLengths{$chr};
print $chr;
for(my $i=0; $i<($len/$windSize); $i++){
    if (!(exists($posArr[$i]))){
#	print $mapfile "0\n";
#	print $gcfile "NA\n";
    } else {
        printf $mapfile "%d\t%d\t%d\t%f\n",$chr,$i,$windSize,$posArr[$i]/$windSize;
	my $basecount =  ($posArr[$i] * $ARGV[4]);
#	print ($gcArr[$i] / ($basecount - $nArr[$i])) . "\n";
    }
}

