#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty qw/:standard/;

print header,&generate_form;
for ('feedback','webmaster','press')
{
	param('to',$_);
	print "<p><a href=",self_url,">$_</a>";
}
print p(url(-full=>1),url(-relative=>1),url(-absolute=>1),url(-path=>1),url(-query=>1));
print end_html;

sub generate_form
{
	my $url=url();
	return start_form(-method=>'get',-action=>$url),
	h1("Please enter your name:"),
	p("Last name",textfield('last')),
	p("first name",textfield('first')),
	p(submit),
	end_form;
}
