#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: <indexed bams>\n" unless @ARGV;
for my $i(@ARGV) {
    print("$i: ");
    !system("samtools idxstats $i | perl -ane '\$m+=\$F[2];\$u+=\$F[3];END{print \$u/(\$u+\$m),\"\\n\"}'") or die "samtools: $!\n";
}
