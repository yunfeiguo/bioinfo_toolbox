#!/usr/bin/env perl
use FindBin qw/$RealBin/;
use lib File::Spec->catdir($RealBin,"..","lib","perl");
use YG::Utils;
use Getopt::Std;
my %opts;
getopts('m:', \%opts);

die "Usage: $0 [-m INT] <fa> <min length INT>\n".
    "	-m INT	optional max length\n" unless @ARGV>=2;
$fa = $ARGV[0];
$fai="$fa.fai";
$l=$ARGV[1];
if (not -f $fai) {
    !system("samtools faidx $fa") or die "$!";
}
open IN,$fai or die "$!";
my @filtered_regions;
my $total = 0;
my $filtered_count = 0;
while(<IN>) {
    @f=split;
    $total++;
    my $pass = $f[1]>$l;
    if (defined $opts{m}) {
	$pass = $pass and $f[1]<=$opts{m};
    }
    if ($pass) {
	push @filtered_regions,"$f[0]:1-$f[1]";
    } else {
	$filtered_count++;
    }
}
close IN;
&YG::Utils::get_fasta_seq($fa, @filtered_regions);
warn "NOTICE: total reads: $total, filtered reads (<= $l bp".
     (defined $opts{m}? ", > $opts{m} bp":"").
     "): $filtered_count\n";
