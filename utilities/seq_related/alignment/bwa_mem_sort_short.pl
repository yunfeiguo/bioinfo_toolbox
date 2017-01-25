#!/usr/bin/env perl
#simple wrapper for bwa mem 0.7.15 and SAMtools 1.3.1
die "Usage: $0 <ref prefix> <output prefix> <fq1> [fq2...]\n" unless @ARGV >= 3;

$ref = shift @ARGV;
$prefix = shift @ARGV;
@fq = @ARGV;
chomp(my $cpu_count = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
$cpu_count *= 0.5; #use only 50% of CPU cores available
$cpu_count++; #at least 1 CPU core

!system("bwa mem -R '\@RG\\tID:$prefix.ID\\tSM:$prefix.SM' -t $cpu_count $ref @fq | samtools sort -o $prefix.sort.bam --output-fmt BAM -") or die "bwa or samtools failed: $!\n";
