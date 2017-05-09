#!/usr/bin/env perl
use strict;
use warnings;

for my $bam(@ARGV) {
  !system("samtools fasta $bam") or die "$!";
}
