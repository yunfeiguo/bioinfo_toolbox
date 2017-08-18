#!/usr/bin/env perl
use File::Basename;
use Getopt::Std;
#simple wrapper for bwa mem 0.7.15 and SAMtools 1.3.1
my $usage = "Usage: $0 [-ol] <ref prefix> <output prefix> <fq1> [fq2...]\n".
    "          -o STRING additional parameters to bwa mem (other than -R, -t)\n".
    "	       -n use ngmlr instead of bwa-mem\n".
    "          -l use '-x pacbio', mutually exclusive with -o\n";
die $usage unless @ARGV >= 3;

my %opts;
getopts('o:ln', \%opts);
if (defined $opts{o} and defined $opts{l}) {
    die $usage;
}
if (defined $opts{o} and $opts{o}=~/-R|-t/) {
    die "-R, -t will be automatically added.\n";
}
if (defined $opts{n} and defined $opts{o} and $opts{o}=~/-R|-t/) {
    die "-q, -r, -t will be automatically added.\n";
}
$ref = shift @ARGV;
$prefix = shift @ARGV;
@fq = grep {/\.(fastq|fq|fa|fasta)$/} @ARGV;
@fqgz = grep {/\.(fa\.gz|fasta\.gz|fastq\.gz|fq\.gz)$/} @ARGV;
@bam = grep {/\.bam$/} @ARGV;
chomp(my $cpu_count = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
#$cpu_count *= 0.5; #use only 50% of CPU cores available
#$cpu_count++; #at least 1 CPU core

my $sample = basename $prefix;
my $cmd;
if ($opts{n}) {
    $cmd = "ngmlr ".
    (defined $opts{l}? " -x pacbio ":"").
    (defined $opts{o}? $opts{o}:"").
    " -t $cpu_count -r $ref -q /dev/stdin | samtools sort -o $prefix.sort.bam --output-fmt BAM -" or die "bwa or ngmlr or samtools failed: $!\n";
} else {
    $cmd = "bwa mem -R '\@RG\\tID:$sample.ID\\tSM:$sample.SM' ".
    (defined $opts{l}? " -x pacbio ":"").
    (defined $opts{o}? $opts{o}:"").
    " -t $cpu_count $ref /dev/stdin | samtools sort -o $prefix.sort.bam --output-fmt BAM -" or die "bwa or ngmlr or samtools failed: $!\n";
}
warn "Running: $cmd\n";		
open ALIGN,"|-", $cmd or die "bwa or ngmlr or samtools failed: $!\n";
if (@fq) {
    open FQ,"-|","cat @fq" or die "failed to read fastqs: $!\n";
    while(<FQ>) {
	print ALIGN;
    }
    close FQ;
}
if (@fqgz) {
    open FQ,"-|","zcat @fqgz" or die "failed to read gzipped fastqs: $!\n";
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

#ngmlr -t 24 -r /net/hippo/volumes/wadi/shared/prj/guoy28/Downloads/databases/assemblies/ecoli/k12/U00096.3.fasta -q output.contigs.fasta | samtools sort -OBAM -o aug_rnn_asm_on_k12.ngmlr.sort.bam -
##!system("cat @fq | bwa mem -R '\@RG\\tID:$prefix.ID\\tSM:$prefix.SM' -x pacbio -t $cpu_count $ref /dev/stdin | samtools sort -o $prefix.sort.bam --output-fmt BAM -") or die "bwa or samtools failed: $!\n";

