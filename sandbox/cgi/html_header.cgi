#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty qw/Link myheadertag/;

my $cgi=new CGI;

print $cgi->header;
print $cgi->start_html(
	-title=>'a complex html header',
	-author=>'yunfeiguo@usc.edu',
	-xbase=>'http://example.com',
	-target=>'_map_panel',
	-meta=>	{
		keywords => "cgi header",
		description => "how to make a big header",
		msg => 'hello',
	},
	-style => { src=>'/style/fourthage.css' },
	-head => [
		Link({-rel=>'origin',
			-href=>'http://seqmule.com'}),
		myheadertag({-myattr=>'myvalue'})
		]
	);
	print $cgi->end_html;
