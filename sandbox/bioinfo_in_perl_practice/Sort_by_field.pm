package Sort_by_field;

use strict;
use warnings;

#take three arguments: filename,delimiter,<numeric or string fileds to be sorted with reverse or normal sort order>. Normal order means from largest to smallest
#for example
#test.fasta,',',(n1n,s2r,s7n)

sub fsort {
	my ($file,$delim,@fields)=@_;
	die "No fields specified\n" unless @fields;
	my @ids;
	my @ops;
	my @orders;
	my @first;
	my @second;
	my $compare_cmd;
	for (@fields) {
		my ($op,$id,$order)=split('',$_);
		if ($op eq 'n') {
			push @ops,'<=>';
		} elsif ($op eq 's') {
			push @ops,'cmp';
		} else {
			die "Illegal column type used, please use 'n' or 's'\n";
		}
		push @ids,$id-1; #array start from 0 while user supplied column number usually starts from 1
		if ($order eq 'n') {
			push @first,'$first';
			push @second,'$second';
			push @orders,'n';
		} elsif ($order eq 'r') {
			push @first,'$second';
			push @second,'$first';
			push @orders,'r';
		} else {
			die "Illegal order used, please use 'n' or 'r'\n";
		}
	}
	my $count=0;
	while (1) {
		$compare_cmd .= "$first[$count]\[$ids[$count]\] $ops[$count] $second[$count]\[$ids[$count]\]";
		$count++;
		last if $count == @fields;
		$compare_cmd .= " or " unless $count == @fields;
	}
	print $compare_cmd,"\n";
	open IN,"<$file" or die "Cannot open file: $!\n";
	print STDERR "Begin reading $file....";
	#group records by first column of interest
	my @input=<IN>;
	print STDERR "Done\n";
	my $match;
	for (0..$ids[0]) {
		if ($_==$ids[0]) {
			$match .= '([^\t]+)';
			$match .= '.*';
			last;
		} else {
			$match .= '[^\t]+';
		}
		$match .= '\t';
	}
	print $match,"\n";
	my %groups;
	for (@input) {
		my ($col1)=/$match/;
		push @{$groups{$col1}},$_;
	}
	print STDERR "Sorting...";
	for (keys %groups) {
		my @group_sorted=sort {&comp($delim,$compare_cmd)} @{$groups{$_}};
		$groups{$_}=\@group_sorted;
	}
	my @sorted;
	if ($orders[0] eq 'n') {
		for (sort {eval "$a $ops[0] $b"} keys %groups) {
			push @sorted,@{$groups{$_}};
		}
	} else {
		for (sort {eval "$b $ops[0] $a"} keys %groups) {
			push @sorted,@{$groups{$_}};
		}
	}
	print STDERR "Done\n";
	return @sorted;
}
sub comp {
	my ($delim,$cmd) =@_;
	my @first = split( $delim, $a );
	my @second = split( $delim, $b );
	return eval $cmd;
}
1;

=head1 Sort_by_field

Sort_by_field: sort numeric or string fields with normal or reverse order

=head1 SYNOPSIS

use Sort_by_field;

print OUTPUT Sort_by_field::fsort("file_name",'delimiter',"n1n","s3r","s3n","n2r");

=head1 DESCRIPTION

How to specify sort method?

'column type'+'column number'+'order'.

For example, we want to sort tab-delimited test.fasta first by string column 1 in normal order, then by numeric column 2 in reverse order, use the following command:

Sort_by_field::fsort("test.fasta",'\t','s1n','n2r');

=cut
