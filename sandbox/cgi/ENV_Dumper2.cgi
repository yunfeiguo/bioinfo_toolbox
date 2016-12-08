#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty;

my $cgi=new CGI::Pretty;
print $cgi->header('text/html','401 auth failed',-authname=>'Don Ray'),
	$cgi->start_html("env dumper"),
	$cgi->table({-border=>1},
		$cgi->Tr($cgi->th(["name","value"])),
		map { $cgi->Tr($cgi->td([$_,$ENV{$_}])) } sort keys %ENV,
		),
	$cgi->end_html;
