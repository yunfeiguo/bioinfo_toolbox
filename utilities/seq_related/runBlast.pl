#!/usr/bin/env perl
use strict;
use warnings;
my $desc="create blast index if not exist, then run blast";
die "Usage: $0 <db> <query>\n" unless @ARGV == 2;
my $db = shift @ARGV;
my $query = shift @ARGV;
my $nproc = 12;
unless(&existIdx($db)){
    &createIdx($db);
}
&runBlast($db,$query);

sub existIdx {
    my $db = shift;
    return(-f "$db.nhr" && -f "$db.nin" && -f "$db.nsq" && -f "$db.shd");
#000850F.fa
#000850F.fa.00.idx
#000850F.fa.fai
#000850F.fa.nhr
#000850F.fa.nin
#000850F.fa.nsq
#000850F.fa.shd
}
sub createIdx {
    my $db = shift;
    !system("makeblastdb -in $db -dbtype nucl 1>&2") or die "makeblastdb: $!\n";
    !system("makembindex -input $db -iformat blastdb 1>&2") or die "makembindex: $!\n";
}
sub runBlast {
    my $db = shift;
    my $q = shift;
    my $out = "$q.blastout";
    !system("blastn -db $db -query $q -num_threads $nproc -task megablast -max_target_seqs 1 -outfmt '10 qseqid qlen sseqid slen qstart qend sstart send sstrand evalue pident qcovs' > $out") or die "blast failed: $!\n";
    warn "NOTICE: output written to $out\n";
}
