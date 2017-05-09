#!/usr/bin/env perl
die "Usage: $0 <fa> <min length INT>" unless @ARGV==2;
$fa = $ARGV[0];
$fai="$fa.fai";
$l=$ARGV[1];
if (not -f $fai) {
    !system("samtools faidx $fa") or die "$!";
}
open IN,$fai or die "$!";
while(<IN>) {
    @f=split;
    if($f[1]>$l){
	!system("samtools faidx $fa $f[0]") or die "$!";
    }
}
