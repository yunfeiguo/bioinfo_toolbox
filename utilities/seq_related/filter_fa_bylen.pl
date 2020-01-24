#!/usr/bin/env perl
use FindBin qw/$RealBin/;
use lib File::Spec->catdir($RealBin,"..","lib","perl");
use YG::Utils;
use Getopt::Std;
my %opts;

die "Usage: $0 [-m INT] <.fa|.fasta|.fq|.fastq> <min length INT, inclusive>\n".
    "	-m INT	optional max length, exclusive\n" unless @ARGV>=2;
getopts('m:', \%opts);
$fa = $ARGV[0];
$fai="$fa.fai";
$l=$ARGV[1];
my $total = 0;
my $filtered_count = 0;
if ($fa =~ /\.(fa|fasta)$/) {
	if (not -f $fai) {
		!system("samtools faidx $fa") or die "samtools faidx error: $!";
	}
	open IN,$fai or die "$!";
	my @filtered_regions;
	while(<IN>) {
		@f=split;
		$total++;
		my $pass = $f[1]>=$l;
		if (defined $opts{m}) {
			$pass = ($pass and $f[1]<$opts{m});
		}
		if ($pass) {
			push @filtered_regions,"$f[0]:1-$f[1]";
		} else {
			$filtered_count++;
		}
	}
	close IN;
	&YG::Utils::get_fasta_seq($fa, @filtered_regions);
} elsif ($fa =~ /\.(fastq|fq)$/) {
	open IN,$fa or die "$!";
	my $line_number = 0;
	my @current_read;
	while(<IN>) {
		chomp;
		push @current_read, $_;
		if ($line_number % 4 == 3) {
			$total++;
			my $pass = length($current_read[1]) >= $l;
			if (defined $opts{m}) {
				$pass = ($pass and (length($current_read[1]) < $opts{m}));
			}
			if ($pass) {
				print join("\n", @current_read),"\n";
			} else {
				$filtered_count++;
			}
			@current_read = ();
		}
		$line_number++;
	}
}
warn "NOTICE: total reads: $total, filtered reads (< $l bp".
     (defined $opts{m}? ", >= $opts{m} bp":"").
     "): $filtered_count\n";
