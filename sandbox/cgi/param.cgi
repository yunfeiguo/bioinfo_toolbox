#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty;
use English qw/-no_match_vars/;

my $cgi=new CGI;
for ($cgi->param())
{
	$cgi->param($_,ucfirst(lc($cgi->param($_))));
}
$cgi->delete('first');

print $cgi->header();
print $cgi->start_html("Welcome");
print "<h1>Welcome, ",$cgi->param('first')," ",$cgi->param('last'),"</h1>";
print "<p>",$cgi->request_method(),"</p>";
print "<p>",$OSNAME,"</p>";
print $cgi->end_html;
