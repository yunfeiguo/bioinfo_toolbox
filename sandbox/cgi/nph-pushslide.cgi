#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Push qw/:standard/;
use CGI::Carp qw/fatalsToBrowser/;

do_push(
	-next_page=>\&show_slide,
	-last_page=>\&go_back,
	-type=>'dynamic',
	-delay=>1,
	);

sub show_slide
{
	my ($cgi,$count)=@_;
	return undef if $count>2;

	my $slide;
	$slide.=h1("This is slide $count");
	return header().start_html("Slide $count").$slide.end_html();
}

sub go_back
{
	my $url;
	#=$ENV{'HTTP_RERFERER'};
	$url='g.cn' unless defined $url; #otherwise default to the homepage

	#generate a 'refresh' header to redirect the client
	return header(-refresh=>"5;URL=$url",-type=>"text/html"),
	start_html("The end"),
	p({-align=>"center"},"Thanks for watching!"),
	end_html();
}
