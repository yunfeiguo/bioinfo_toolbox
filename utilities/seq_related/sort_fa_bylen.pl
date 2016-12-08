#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

die "Usage: $0 <fasta>\n" unless @ARGV==1;
my $fa = shift;
my $seqio_obj = Bio::SeqIO->new(-file=>$fa,-format=>'fasta');
my @seq_obj;
while(my $i = $seqio_obj->next_seq) {
    push @seq_obj,$i;
}
for my $i(sort {$a->length <=> $b->length} @seq_obj) {
    print ">".$i->id."\n";
    print $i->seq."\n";
}
