#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;
use CGI;

my $cgi=new CGI;
print $cgi->header;
my $cookie1=$cgi->cookie(-name=>"firstcookie",-value=>"abcde");
print "Cookie 1: $cookie1\n";
