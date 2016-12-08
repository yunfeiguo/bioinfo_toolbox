#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty;

my $cgi=new CGI;

print $cgi->header;
print $cgi->start_html;
print $cgi->table
(
	{-border=>1,-cellspacing=>3,-cellpadding=>3},
	$cgi->Tr({-align=>'center',-valign=>'top'},[
		$cgi->th(["col1","col2","col3"]),
		]),
	$cgi->Tr({-align=>'center',-valign=>'middle'},[
		$cgi->td(["red","blue",'yellow']),
		$cgi->td(["green","orange","brown"]),
		]),
	$cgi->caption("example")
);
