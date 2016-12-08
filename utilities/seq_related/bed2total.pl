#!/usr/bin/env perl
die "Usage: $0 <1.bed 2.bed ...>\n" unless @ARGV;
for $i(@ARGV) {
    print "Processing $i:\n";
    $t = 0;
    open IN,"<",$i or die "open: $!\n";
    while (<IN>) {
	@F = split;
	$t += $F[2]-$F[1];
    }
    print $t,"\n";
}
