#!/usr/bin/env perl
use File::Basename;
use Getopt::Std;
#get coverage plot for whole genome (optionally for a region)
my $usage = "Usage: $0 [-r <BED>] <bam1> [bam2...]\n".
"          -r BED bed file for specifying region\n";
die $usage unless @ARGV >= 1;

my %opts;
getopts('r:', \%opts);
my $cmd = "samtools depth -m 9999999 -aa ".($opts{r}? " -b $opts{r}" : "")." @ARGV";

open DEPTH,"-|", $cmd or die "failed to read fastqs: $!\n";
my %len;
while (<DEPTH>) {
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
