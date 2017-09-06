#!/usr/bin/env perl
use File::Basename;
#simple wrapper for bwa mem 0.7.15 and SAMtools 1.3.1
die "Usage: $0 <ref prefix> <output prefix> <fq1> [fq2...]\n" unless @ARGV >= 3;

$ref = shift @ARGV;
$prefix = shift @ARGV;
@fq = grep {/\.(fastq|fq|fa|fasta|fa\.gz|fasta\.gz)$/} @ARGV;
@bam = grep {/\.bam$/} @ARGV;
chomp(my $cpu_count = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
#$cpu_count *= 0.5; #use only 50% of CPU cores available
#$cpu_count++; #at least 1 CPU core

open ALIGN,"|-", "bwa mem -R '\@RG\\tID:".(basename $prefix).".ID\\tSM:$prefix.SM' -x pacbio -t $cpu_count $ref /dev/stdin | samtools sort -o $prefix.sort.bam --output-fmt BAM -" or die "bwa or samtools failed: $!\n";
if (@fq) {
    open FQ,"-|","cat @fq" or die "failed to read fastqs: $!\n";
    while(<FQ>) {
	print ALIGN;
    }
    close FQ;
}
for my $bam(@bam) {
    open BAM,"-|","bamToFastq -i $bam -fq /dev/stdout" or die "failed to read fastqs: $!\n";
    while(<BAM>) {
	print ALIGN;
    }
    close BAM;
}
close ALIGN;

!system("samtools index $prefix.sort.bam") or die "samtools index: $!";

##!system("cat @fq | bwa mem -R '\@RG\\tID:$prefix.ID\\tSM:$prefix.SM' -x pacbio -t $cpu_count $ref /dev/stdin | samtools sort -o $prefix.sort.bam --output-fmt BAM -") or die "bwa or samtools failed: $!\n";
