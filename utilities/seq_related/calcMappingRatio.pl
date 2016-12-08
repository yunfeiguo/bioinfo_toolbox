#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <1.bam 2.bam ...>\n" unless @ARGV;
for my $i(@ARGV) {
    if (not -e "$i.bai") {
	!system("samtools index $i") or die "samtools index: $!\n";
    }
    print("Mapping ratio for $i:\n");
    !system("samtools idxstats $i | perl -ane '\$mapped += \$F[2]; \$unmapped = \$F[3] if \$F[0] eq \"*\"; END{\$total=\$mapped+\$unmapped;print \"\$unmapped/\$total (\".(\$unmapped/\$total).\") unmapped reads\\n\"}'") or die "samtools idxstats: $!\n";
}

