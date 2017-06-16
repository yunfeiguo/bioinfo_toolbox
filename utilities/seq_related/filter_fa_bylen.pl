#!/usr/bin/env perl
use FindBin qw/$RealBin/;
use lib File::Spec->catdir($RealBin,"..","lib","perl");
use YG::Utils;

die "Usage: $0 <fa> <min length INT>" unless @ARGV==2;
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
    if($f[1]>$l){
	push @filtered_regions,"$f[0]:1-$f[1]";
    } else {
	$filtered_count++;
    }
}
close IN;
&YG::Utils::get_fasta_seq($fa, @filtered_regions);
warn "NOTICE: total reads: $total, filtered reads (<= $l bp): $filtered_count\n";
