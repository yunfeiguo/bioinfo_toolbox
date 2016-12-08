#!/usr/bin/env perl

use strict;
use warnings;
use CGI::Push qw/:standard/;
use CGI::Carp qw/fatalsToBrowser/;

my $line="";
my $delay=1;
my $total_delay=1.5;
do_push(
    -next_page=>\&refresh,
	-delay=>$delay,
    -last_page=>\&done,
);

sub refresh
{
	my ($cgi,$count)=@_; #passed in by CGI::Push
	return undef if $count>33;
	my $page=start_html().p("The count is $count");
	if (length($line)>9)
	{
	    $line="";
	} else
	{
	    $line.="*";
	}
	$page.=p($line."\n").end_html();
	push_delay($total_delay-push_delay());
	return $page;
}

sub done
{
    my ($cgi,$count)=@_;
    return start_html."Count stopped on $count".end_html;
}
