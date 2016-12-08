package Sort_by_field_v2;

use strict;
use warnings;

#take three arguments: filename,delimiter,<numeric or string fileds to be sorted with reverse or normal sort order>. Normal order means from largest to smallest
#for example
#test.fasta,',',(n1n,s2r,s7n)

sub fsort {
	my ($file,$delim,@fields)=@_;
	print STDERR "No fields specified\n" and return 0 unless @fields;
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
			push @first,'${$a}';
			push @second,'${$b}';
			push @orders,'n';
		} elsif ($order eq 'r') {
			push @first,'${$b}';
			push @second,'${$a}';
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
	my @input;
	while (<IN>) {
		my @line=split /\t/;
		push @input,\@line;
	}
	#seek IN,0,0;
	#my @input2=<IN>;
	close IN;
	print STDERR "Done\n";

#	my $match;
#	for (0..$ids[0]) {
#		if ($_==$ids[0]) {
#			$match .= '([^\t]+)';
#			$match .= '.*';
#			last;
#		} else {
#			$match .= '[^\t]+';
#		}
#		$match .= '\t';
#	}
#	print $match,"\n";
#	my %groups;
#	for (1..@input) {
#		my ($col1)= $input2[$_-1]=~/$match/;
#		push @{$groups{$col1}},$input[$_-1];
#	}
#	print STDERR "Sorting...";
#	for (keys %groups) {
#		my @group_sorted=sort { eval $compare_cmd} @{$groups{$_}};
#		$groups{$_}=\@group_sorted;
#	}
#	my @sorted_ref;
#	if ($orders[0] eq 'n') {
#		for (sort {eval "$a $ops[0] $b"} keys %groups) {
#			push @sorted_ref,@{$groups{$_}};
#		}
#	} else {
#		for (sort {eval "$b $ops[0] $a"} keys %groups) {
#			push @sorted_ref,@{$groups{$_}};
#		}
#	}
#	print STDERR "Done\n";
	#convert references to arrays
	my @sorted;
	for (@sorted_ref) {
		my $joined=join '\t',@{$_};
		push @sorted,$joined;
	}
	return @sorted;
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
