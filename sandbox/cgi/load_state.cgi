#!/usr/bin/env perl

use strict;
use warnings;
use CGI;
use CGI::Carp qw/fatalsToBrowser/;

my $cgi;
if (open (my $STATE,"state"))
{
	$cgi= new CGI ($STATE);
	close $STATE;
}
print $cgi->header,$cgi->start_html;
print $cgi->p($cgi->param());
print $cgi->end_html;

