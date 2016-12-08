#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Text::Wrap;

die "Usage: $0 <1.subreads_filtered.fasta 2.fasta ...>\n" unless @ARGV>=1;
my %global_fh_count; #store all opened filehandle for output
$Text::Wrap::columns = 71;
#example header
#>m150108_004913_42146_c100718992550000001823154405141501_s1_p0/54499/0_8312
#the part between > and / is the SMRT cell ID, which is used as file name
warn "NOTICE: all files\n";
map {warn "$_\n" } @ARGV;

for my $input(@ARGV) {
    warn "NOTICE: began processing $input\n";
    my $seqio  = Bio::SeqIO->new( '-format' => 'fasta' , -file => $input);
    my %all_fh; #store all opened filehandle for output
    while (my $seqobj = $seqio->next_seq() ) {

	my $fh = &get_fh_byid(&get_smrtid($seqobj->id),\%all_fh,\%global_fh_count);
	#it seems wrapping is necessary for falcon to work
	print $fh ">",$seqobj->id,"\n",wrap('','',$seqobj->seq),"\n";
    }
    #close all opened fh for writing
    close $_ for(values %all_fh);
    warn int(keys %all_fh)," FASTA written for $input\n";
}

warn "Read count for all SMRT cells\n";
my ($key,$value);
print "$key: $value\n" while (($key,$value) = each %global_fh_count);
warn "All done\n";


###############################SUBROUTINES##################################################
sub get_smrtid {
    my $header = shift;
    my ($id) = $header =~ m%^(\w+)/%;
    die "failed to parse $header \n" unless $id;
    return $id;
}
sub get_fh_byid {
    my $id = shift;
    my $fh_hashref = shift;
    my $global_fh_hashref = shift;

    unless(defined $fh_hashref->{$id}) {
	my $file = "$id.fasta";
	if(defined $global_fh_hashref->{$id} and $global_fh_hashref->{$id} >= 1) {
	    warn "WARNING: $id shows up in different input files\n".
	    "WARNING: This may suggest duplicate SMRT IDs.\n".
	    "WARNING: $file will be overwritten.\n";
	}
	open my $fh,'>',$file or die "Failed to open $file: $!\n";
	$fh_hashref->{$id} = $fh;
    }
    if(defined $global_fh_hashref->{$id}) {
	$global_fh_hashref->{$id} ++;
    } else {
	$global_fh_hashref->{$id} = 1;
    }
    return $fh_hashref->{$id};
}
