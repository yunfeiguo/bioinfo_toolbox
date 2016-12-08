#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty qw/:standard fruit fruits/;

print header("text/html"),
	fruits(
		fruit ({-size=>'small',-color=>'red'},["strawberry","cherry"])
	);
