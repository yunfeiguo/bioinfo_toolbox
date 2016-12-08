#!/usr/bin/env perl

use strict;
use warnings;
use CGI;
use CGI::Carp qw/fatalsToBrowser/;

my $cgi=new CGI;

my $cookie=$cgi->cookie("myCookie");

if ($cookie)
{
	print $cgi->header();
	#no need to send cookie again
} else
{
	my $value=generate_unique_id();
	$cookie=$cgi->cookie({-name => "myCookie",-value => $value,-expires => "+1d"});
	print $cgi->header({-type => "text/html",-cookie => $cookie});
}

sub generate_unique_id
{
	return sprintf("%08.8x",rand()*0xffffffff);
}
