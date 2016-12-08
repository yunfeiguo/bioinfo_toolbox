#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;

my $date;
my $format="%s";

unless (open DATE,"-|")
{
	exec '/bin/date','-u',"+$format";
}

$date=<DATE>;
close DATE;

print "Content-Type: text/html\n\n";
print "<H1>$$ $date</H1>";

my $result=open (DATE,"-|");
exec '/bin/date','-u',"+$format" unless $result;
$date=<DATE>;
close DATE;
print "<H1>$$ date2: $date</h1>";

open (DATE,"-|") || exec '/bin/date','-u',"+$format";
$date=<DATE>;
close DATE;

print "<h1>$$ date 3: $date</h1>";
