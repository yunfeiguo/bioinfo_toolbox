package YG::Utils;
use Carp;
use strict;
use warnings;

sub parse_sam_record {
    #parse a line of SAM record
    croak "Usage: &parse_sam_record(SAM_RECORD)" unless @_ == 1;
    my @f = split /\t/,$_[0];
    return {	query_name => $f[0],
    	    	flag => $f[1],
		reference_name => $f[2],
		position => $f[3],
		mapping_quality => $f[4],
		cigar => &parse_cigar_str($f[5]),
		next_reference_name => $f[6],
		next_position => $f[7],
		template_length => $f[8],
		sequence => $f[9],
		quality => $f[10],
		tags => &parse_sam_tags(@f[11..$#f]),
	  	};
}
sub parse_sam_tags {
    my $tag = {};
    for my $tag_field(@_) {
	my @f = split /:/,$tag_field or next;
	$tag->{$f[0]} = $f[2];
    }
    return $tag;
}
sub get_idt {
    croak "Usage: &get_idt(SAM_RECORD)" unless @_ == 1;
    #percent identity definition: #matches/total_bp
    my @f = split /\t/, $_[0];
    my $cigar = &parse_cigar_str($f[5]);
    my ($nm) = $_[0] =~ /NM:i:(\d+)/;
    $nm = $nm ? $nm : 0;
    my $ins = defined $cigar->{tabulation}->{I} ? $cigar->{tabulation}->{I} : 0;
    my $del = defined $cigar->{tabulation}->{D} ? $cigar->{tabulation}->{D} : 0;
    my $s = (defined $cigar->{tabulation}->{'M'} ? $cigar->{tabulation}->{'M'} : 0) + 
            (defined $cigar->{tabulation}->{'='} ? $cigar->{tabulation}->{'='} : 0);
    my $clip = (defined $cigar->{tabulation}->{'H'} ? $cigar->{tabulation}->{'H'} : 0) + (defined $cigar->{tabulation}->{'S'} ? $cigar->{tabulation}->{'S'} : 0);
    my $match_len = $s - ($nm - $ins - $del);
    my $aligned_length = $cigar->{len} - $clip;

    #calculate idt only on aligned portion of the read
    return $aligned_length == 0 ? (0,$match_len,$aligned_length) : 
        ($match_len / $aligned_length, $match_len, $aligned_length);
}
sub get_clip {
    my $cigar = shift;
    my $parsed_cigar = &parse_cigar_str($cigar);
    return (defined $parsed_cigar->{tabulation}->{'H'} ? $parsed_cigar->{tabulation}->{'H'} : 0) + (defined $parsed_cigar->{tabulation}->{'S'} ? $parsed_cigar->{tabulation}->{'S'} : 0);
}
sub parse_cigar_str {
    croak "Usage: &parse_cigar_str('3S11M5S')" unless @_ == 1;
    my $cigar_str = shift;
    my @cigars = $cigar_str =~m/(\d+)([MIDNSHP=X])/g;
    my %tabulation;
    my $i = 0;
    my $len = 0;
    while($i <= $#cigars) {
	$tabulation{$cigars[$i+1]} += $cigars[$i];
	$len += $cigars[$i];
	$i += 2;
    }
    my %parsed_cigar;
    $parsed_cigar{array} = \@cigars;
    $parsed_cigar{len} = $len;
    $parsed_cigar{tabulation} = \%tabulation;
    return \%parsed_cigar;
}
sub get_fasta_seq{
croak "Usage: $0 <FASTA> <1:1-1000 2:3-999 ...>\n" unless @ARGV>=2;

my $fa=shift;
my $idx="$fa.fai";
my @regions=@_;

croak "ERROR: index file missing or empty\n" unless -s $idx;

my %contig=&readIdx($idx);


open IN,'<',$fa or croak "Failed to read $fa: $!\n";

for my $i(@regions)
{
    croak "1:1-1000 format expected: $i\n" unless $i=~/(.*):(\d+)-(\d+)/;
    my ($id,$start,$end)= ($1,$2,$3);
    my $len;
    my $offset;
    my $return;
    my $nbreak_start;
    my $nbreak_end;

    #check if contig exists
    croak "$id doesn't exist\n" unless exists $contig{$id};
    #check if the region is out of bound
    croak "Region out of bound or end smaller than start: $i\nContig length: $contig{$id}{length}\n" 
    unless ($end>=1 && 
	$start>=1 && 
	$end >= $start && 
	$end <= $contig{$id}{length});

    my $nchar_ln=$contig{$id}{nchar_ln};

    #atcG\n
    #Ccga\n
    #G and C should have different nbreaks
    $nbreak_start=&getNbreak($start,$nchar_ln);
    $nbreak_end=&getNbreak($end,$nchar_ln);

    #contig start	  position	line breaks
    $offset=$contig{$id}{offset}+($start-1) + $nbreak_start;

    $len= ($end-$start+1) + ($nbreak_end-$nbreak_start);

    seek IN,$offset,0;
    read IN,$return,$len;

    $return=~s/[\r\n]+//g;
    print ">$i\n";
    print "$return\n";
}

close IN;
}

###################SUBROUTINES##########################
sub read_fasta {
    croak "usage: &read_fasta('fasta')" unless @_ == 1;
    #read entire fasta in memory
    my $fa = shift;
    my @seqs;

    my $seq;
    my $id;
    open IN,'<',$fa or die "$!";
    while(<IN>) {
	s/\s+$//;
	if(/^>(.*)/) {
	    if (defined $seq) {
		push @seqs, [$id, $seq] if (defined $id) ;
		$seq = undef;
	    }
	    $id = $1;
	} else {
	    $seq .= $_;
	}
    }
    if (defined $seq) {
	push @seqs, [$id, $seq] if (defined $id) ;
    }
    close IN;
    return \@seqs;
}
sub readIdx
{
    my $fai=shift;
    my %return;
    open IN,"<",$fai or croak "Failed to read $fai: $!\n";
    while (<IN>)
    {
	my @f=split /\t/;
	croak "5 fields expected at line $. of $fai: $_\n" unless @f==5;
	my ($id,$len,$offset,$nchar_ln,$nbyte_ln)=@f;

	$return{$id}={
	    length=>$len, #length of contig
	    offset=>$offset, #offset where first character in that contig appears
	    nchar_ln=>$nchar_ln, #number of characters per line
	    nbyte_ln=>$nbyte_ln, #number of bytes per line
	};
    }
    close IN;
    return %return;
}
sub getNbreak
{
    my $nbreak;
    my ($start,$nchar_ln)=@_;
    if ($start % $nchar_ln == 0)
    {
	$nbreak=$start/$nchar_ln-1;
    } else
    {
	$nbreak=int($start/$nchar_ln);
    }
    return $nbreak;
}
1;
