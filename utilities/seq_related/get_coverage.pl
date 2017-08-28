#!/usr/bin/env perl
use File::Basename;
use Getopt::Std;
#get coverage plot for whole genome (optionally for a region)
my $usage = "Usage: $0 [options] <bam1> [bam2...]\n".
"	   -i <GENOME> treat reads as continuous intervals (CIGAR ignored, non-primary considered). requires GENOME file for bedtools\n".
"          -r <BED> bed file for specifying region\n".
"	   -v output per-base depth to STDERR\n";
die $usage unless @ARGV >= 1;

my %opts;
getopts('i:r:v', \%opts);
my $cmd;
if ($opts{i}) {
    $cmd = "samtools view @ARGV | "."perl -ane '".'@n=$F[5]=~/(\d+)[M=D]/g;$l=0;map{$l+=$_}@n;print join("\t",$F[2],$F[3]-1,$F[3]-1+$l),"\n"'."'";
    if ($opts{r}) {
	$cmd .= " | bedtools intersect -a - -b $opts{r} ";
    }
    $cmd .= " | bedtools genomecov -i - -g $opts{i} -d ";
} else {
    $cmd = "samtools depth -m 9999999 -aa ".($opts{r}? " -b $opts{r}" : "")." @ARGV";
}

open DEPTH,"-|", $cmd or die "failed to read from samtools or bedtools: $!\n";
my %len;
while (<DEPTH>) {
    print STDERR if $opts{v};
    s/\s+$//;
    my @F = split /\t/;
    if( !defined($len{$F[2]}) ) {
	$len{$F[2]} = 1;
    } else {
	$len{$F[2]}++;
    }
}
close DEPTH;
my @ints = sort { $a <=> $b } keys(%len);
my %bp_count = (5 => 0, 10 => 0, 20 => 0, 30 => 0);
for my $i(@ints) {
    for $j(keys %bp_count) {
	if ($i <= $j) {
	    $bp_count{$j} += $len{$i}
	}
    }
    print "$i\t$len{$i}\n" if $len{$i};
}

for my $i(sort {$a <=> $b} keys %bp_count) {
    print "#bp covered <= $i reads: $bp_count{$i}\n";
}
