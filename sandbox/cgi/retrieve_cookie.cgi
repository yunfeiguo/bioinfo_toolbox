#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;
use CGI qw/:standard/;

my ($cookies,@cookies,%cookies);
$cookies=$ENV{HTTP_COOKIE} || $ENV{COOKIE};
if ($cookies)
{
	@cookies=split /;\s/,$cookies;
	for (@cookies)
	{
		/([^=]+)=(.*)/ and $cookies{$1}=$2;
	}
}

print header(-type=>"text/html");
for (keys %cookies)
{
	print p("$_ : $cookies{$_}");	
}

#alternative method
my $cgi=new CGI;
print p($cgi->cookie("mycookie"));
