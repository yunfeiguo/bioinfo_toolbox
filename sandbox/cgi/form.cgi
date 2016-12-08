#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty qw/-debug/;

my $cgi=new CGI::Pretty;
print $cgi->header,$cgi->start_html;
if ($cgi->param('first') && $cgi->param('last'))
{
	my $first=ucfirst(lc($cgi->param('first')));
	my $last=ucfirst(lc($cgi->param('last')));
	print $cgi->h1("hello, $first $last");
}
else
{
	print $cgi->h1("enter your name");
	if ($cgi->param('first') || $cgi->param('last'))
	{
		print $cgi->center($cgi->font({-color=>'red'},"You must enter a",
				($cgi->param('first')?"last":"first"),"name"));
	}
	print &generate_form;
}
print $cgi->end_html;

if (open (my $STATE,">","state"))
{
	$cgi->save($STATE);
	close $STATE;
}

sub generate_form
{
	return 
$cgi->start_form,
$cgi->h1("Please enter your name"),
$cgi->p("Last name",$cgi->textfield("last")),
$cgi->p("First name",$cgi->textfield("first")),
$cgi->p($cgi->submit()),
$cgi->end_form;
}
