#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <.fa>\n" unless @ARGV == 1;
my $fa=shift @ARGV;
my $fai = "$fa.fai";
!system("samtools faidx $fa") or die "faidx: $!\n" unless -e $fai;
my @a=`cat $fai`;
for(@a){
    my $i=$_;
    $i=~s/^(\S+).*/$1/s;
    warn "processing $i\n";
    !system("samtools faidx $fa $i > $i.fa") or die "$!\n";
}
