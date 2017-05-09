#!/usr/bin/env perl

use strict;
use warnings;

die "Usage: <seq1> <seq2>\n" unless @ARGV==2;

my ($seq1,$seq2)=@ARGV;

my $MATCH=10;
my $MISMATCH=-2;
my $GAP=-5;

#initialization
my @matrix;
$matrix[0][0]{score}=0;
$matrix[0][0]{pointer}="none";
for(my $j=1;$j<=length $seq1;$j++)
{
	$matrix[0][$j]{score}=$GAP*$j;
	$matrix[0][$j]{pointer}="left";
}
for (my $i=1;$i<=length $seq2;$i++)
{
	$matrix[$i][0]{score}=$GAP*$i;
	$matrix[$i][0]{pointer}="up";
}

#fill
for(my $i=1;$i<=length $seq2;$i++)
{
	for(my $j=1;$j<=length $seq1;$j++)
	{
		my($diagonal_score,$left_score,$up_score);

		#calculate match score
		my $letter1=substr($seq1,$j-1,1);
		my $letter2=substr($seq2,$i-1,1);
		if($letter1 eq $letter2)
		{
			$diagonal_score=$matrix[$i-1][$j-1]{score}+$MATCH;
		}else
		{
			$diagonal_score=$matrix[$i-1][$j-1]{score}+$MISMATCH;
		}
		#calculate gap score
		$up_score=$matrix[$i-1][$j]{score}+$GAP;
		$left_score=$matrix[$i][$j-1]{score}+$GAP;
		#choose max score
		if ($diagonal_score>=$up_score)
		{
			if ($diagonal_score>=$left_score)
			{
				$matrix[$i][$j]{score}=$diagonal_score;
				$matrix[$i][$j]{pointer}="diagonal";
			}
		} else
		{
			if ($up_score >= $left_score)
			{
				$matrix[$i][$j]{score}=$up_score;
				$matrix[$i][$j]{pointer}="up";
			}
		}

		if ($left_score >= $diagonal_score)
		{
			if ($left_score >= $up_score)
			{
				$matrix[$i][$j]{score}=$left_score;
				$matrix[$i][$j]{pointer}="left";
			}
		}
	}
}
#trace-back

my $align1="";
my $align2="";

my $j=length $seq1;
my $i=length $seq2;

while(1)
{
	last if $matrix[$i][$j]{pointer} eq "none"; #ends at first cell

	if ($matrix[$i][$j]{pointer} eq "diagonal")
	{
		$align1.=substr($seq1,$j-1,1);
		$align2.=substr($seq2,$i-1,1);
		$i--;
		$j--;
	} elsif ($matrix[$i][$j]{pointer} eq 'left')
	{
		$align1.=substr($seq1,$j-1,1);
		$align2.='-';
		$j--;
	} elsif ($matrix[$i][$j]{pointer} eq 'up')
	{
		$align1.='-';
		$align2.=substr($seq2,$i-1,1);
		$i--;
	}
}
$align1=reverse $align1;
$align2=reverse $align2;
print "$align1\n";
print "$align2\n";
