#!/usr/bin/env perl
use warnings;
use strict;
use FindBin qw/$RealBin/;
use lib File::Spec->catdir($RealBin,"..","lib","perl");
use YG::Utils;

die "Usage: $0 <kmer size> <fa>" unless @ARGV==2;
my $k = $ARGV[0];
my $fa = $ARGV[1];

warn "counting all ".$k."mers in $fa\n";
my %kmers = ();
my $fa_records = &YG::Utils::read_fasta($fa);

for my $j(0..$#{$fa_records}) {
    my ($id, $seq) = @{$fa_records->[$j]};
    my $l = length $seq;
    for (my $i = 0; $i + $k <= $l; $i++) {
	$kmers{substr($seq, $i, $k)}++;
    }
}

print "KMER\tcount\n";
for my $kmer(sort keys %kmers) {
    print uc($kmer),"\t$kmers{$kmer}\n";
}
warn "all done.\n";
