#!/usr/bin/perl -w
# output sequence length and count and summary stats
# 
# AUTHOR: Yunfei Guo (modified from script by Joseph Fass)
# LAST REVISED: March 2015
use strict;
use POSIX;
use Getopt::Std;

my $usage = "\n\nusage: $0 <fasta file(s)>\n";
our($opt_i);# histogram interval
$opt_i = 1;
getopts('i:') or die $usage;
if (!defined($opt_i) or !($opt_i =~ m/^[0-9]+$/)) {$opt_i = 100;}

if( ( $#ARGV + 1 ) < 1 ) {
    die $usage;
}

# Read in sequences from one or more fasta files
my @data_files = @ARGV;
my $debug = 0;
my $Id;
# Count the number of sequences in the file and create a histogram of the distribution
my $n = 0;
my $individual_len = 0;
my $int;
my $totalLength = 0;
my $gcCount = 0;
my %len= ();
my @seqLengths;
foreach my $file (@data_files){
    warn "NOTICE: reading $file...\n";
    if($file =~ /\.gz$/) {
	open(FASTA, "gunzip -c $file |") or die"Can't open file $file\n";
    } else {
	open(FASTA, $file) or die"Can't open file $file\n";
    }
    while (<FASTA>) {
	if (/^>(.*)$/)  {
	    print if $debug;
	    $n++;
	    if(defined $Id) {
		push @seqLengths, $individual_len; # record length for N50 calc's
		if( !defined($len{$individual_len}) ) {
		    $len{$individual_len} = 1;
		} else {
		    $len{$individual_len}++;
		}
		$Id = undef;
		$individual_len = 0;
	    }
	    $Id = $1; 
	} elsif (/^(\S+)[\r\n]*$/) 	{
	    print if $debug;
	    chomp $_;
	    if(defined $Id){ #prevent counting when there are orphan sequences at the begining
		my $tmp_len += length($_);
		$individual_len += $tmp_len;
		$totalLength += $tmp_len;
		$gcCount += ($_ =~ tr/gGcC/gGcC/);
	    }
	}
    }
    close (FASTA);
}
if(defined $Id) {
    push @seqLengths, $individual_len; # record length for N50 calc's
    if( !defined($len{$individual_len}) ) {
	$len{$individual_len} = 1;
    } else {
	$len{$individual_len}++;
    }
    $Id = undef;
	$individual_len = 0;
}

# Calculate N25, N50, and N75 and counts
my $N25; my $N50; my $N75;
my $N25count=0; my $N50count=0; my $N75count=0;
my $frac_covered = $totalLength;
@seqLengths = reverse sort { $a <=> $b } @seqLengths;
my @allSortedLen = @seqLengths;
$N25 = $seqLengths[0];
while ($frac_covered > $totalLength*3/4) {
    $N25 = shift(@seqLengths);
    $N25count++; $N50count++; $N75count++;
    $frac_covered -= $N25;
}
$N50 = $N25;
while ($frac_covered > $totalLength/2) {
    $N50 = shift(@seqLengths);
    $N50count++; $N75count++;
    $frac_covered -= $N50;
}
$N75 = $N50;
while ($frac_covered > $totalLength/4) {
    $N75 = shift(@seqLengths);
    $N75count++;
    $frac_covered -= $N75;
}

# Print out the results
print "\n";
my @ints = sort { $a <=> $b } keys(%len);
for my $i(@ints) {
    print "$i\t$len{$i}\n" if $len{$i};
}
print "\n";
print "Total length of sequences:\t$totalLength bp\n";
print "Total number of sequences:\t$n\n";
printf "Average length:\t\t\t%.2f\n",($n==0? 0:$totalLength/$n);
printf "Median length:\t\t\t%.2f\n",(&getMedian(@allSortedLen));
# not sure if these right wrt N25 and N75 ..
print "N25 stats:\t\t\t25% of total sequence length is contained in the ".$N25count." sequences >= ".$N25." bp\n" if defined $N25count && defined $N25;
print "N50 stats:\t\t\t50% of total sequence length is contained in the ".$N50count." sequences >= ".$N50." bp\n" if defined $N50count && defined $N50;
print "N75 stats:\t\t\t75% of total sequence length is contained in the ".$N75count." sequences >= ".$N75." bp\n" if defined $N75count && defined $N75;
print "Total GC count:\t\t\t$gcCount bp\n" if defined $gcCount;
printf "GC %%:\t\t\t\t%.2f %%\n", ($gcCount/$totalLength * 100) if $totalLength != 0;
print "\n";

sub getMedian {
    return if @_ == 0;
    if(@_ % 2 == 0) {
	my $middle = @_/2;
	return(($_[$middle-1]+$_[$middle])/2);
    } else {
	my $middle = (@_-1)/2;
	return($_[$middle]);
    }
}
