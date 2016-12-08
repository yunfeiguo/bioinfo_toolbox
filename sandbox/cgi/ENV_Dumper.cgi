#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty;
use lib "/opt/perl/lib/site_perl/5.14.2/";
use URI::Escape;

print "Content-type: text/html\n\n";
print "<html><head><title>Environment Dumper </title></head><body>";
print "<center><table border=1>";
print "<tr><td>$_</td><td>".uri_unescape($ENV{$_})."</td></tr>" for (sort keys %ENV);
print "</table></center></body></html>";
