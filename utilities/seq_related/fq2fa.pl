#!/usr/bin/env perl

die "Usage: $0 <1.fq[.gz]> [2.fq.[.gz]] ...>" unless @ARGV >= 1;
for $fq(@ARGV) {
	die "unrecognized suffix\n" unless $fq=~/\.fq$|\.fastq$|\.fastq\.gz$|\.fq\.gz$/;
	$read_cmd = $fq =~ /\.gz$/? "gunzip -c ":"cat ";
	!system( "$read_cmd $fq |paste - - - - | perl -ane '".'$F[0]=~s/^@/>/;print $F[0],"\n",$F[1],"\n"'."'" ) or die "$!\n";
}
warn "converted ",scalar(@ARGV), " files to fasta.\n";
