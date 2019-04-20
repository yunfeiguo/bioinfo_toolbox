#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <fasta>\n" unless @ARGV==1;
warn "Scanning for gaps...\n";
my $fa = shift;
my $count = 0;
open FA,"<",$fa or die "failed to open $fa: $!\n";
my $id = "";
my $seq = "";
while(my $i=<FA>) {
    if ($i =~ /^>/) {
	if ($id) {
	    process($id, $seq);
	}
	($id) = $i =~ />(.*)/;
	$seq = "";
    } else {
	chomp $i;
	$seq .= $i;
    }	
}
if ($id) {
    process($id, $seq);
}
sub process {
	#print "<<@_>>";
	($id, $seq) = @_;
    my $pos = 0; #position of probe
	my $ingap = 0; #whether we are in a gap
	for my $j(split //,$seq) {
	    if ($j eq 'N' or $j eq 'n') {
		print $id,"\t$pos\t" if $ingap == 0;
		$ingap = 1;
	    } else {
		if($ingap) {
		    $ingap = 0;
		    print $pos,"\tGap$count\n";
		    $count++;
		} else {
		    1;
		}
	    }
    $pos++;
}
if($ingap) {
    print $pos,"\tGap$count\n";
    $count++;
}
}
