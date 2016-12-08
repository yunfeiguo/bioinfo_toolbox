#!/usr/bin/env perl
use strict;
use warnings;
die "Usage: $0 <.fai>\n" unless @ARGV==1;
my $in = shift @ARGV;
open IN,'<',$in or die "open($in): $!\n";
while(<IN>) {
    my @f=split;
    print join("\t",$f[0],0,$f[1]),"\n";
}
close IN;
