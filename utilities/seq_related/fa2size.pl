#!/usr/bin/env perl
use strict;
use warnings;
die "Usage: $0 <fasta>\n" unless @ARGV==1;
my $file = shift @ARGV;
my $fai = "$file.fai";
unless(-e $fai or -l $fai) {
    warn "Call SAMtools to create FASTA index first\n";
    !system("samtools faidx $file") or die "samtools $file indexing failed: $!\n";
}
die "$fai missing!\n" unless -e $fai or -l $fai;
!system("cut -f 1,2 $fai | sort -k 1,1 -k2,2n") or die "cut $fai: $!\n";
warn "All done\n";
