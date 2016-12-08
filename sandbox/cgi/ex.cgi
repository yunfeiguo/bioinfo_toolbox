#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;

print "Content-type: text/plain\n\n";
print "Hello CGI world!\n";
print "You're calling from ",$ENV{SERVER_ADDR},"\n";
