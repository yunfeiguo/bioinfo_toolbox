#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Pretty qw/:standard/;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Cookie;

print header;
my $cookie2=CGI::Cookie->new(-name=>"secondcookie",-value=>"xxxooo");
print "Cookie2: $cookie2\n";
