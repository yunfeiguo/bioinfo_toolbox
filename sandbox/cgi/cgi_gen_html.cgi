#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty;

my $cgi=new CGI;

print $cgi->header,$cgi->start_html("Simple html");
print $cgi->center("Centered text");
print $cgi->p("a paragraph");
print $cgi->br;
print CGI::b("bold"),CGI::i("italic");
print $cgi->p("a para",$cgi->sup("a superscript"));

print $cgi->ol({-type=>'9'}),$cgi->li(["ONE","two"]);
print $cgi->end_html;
