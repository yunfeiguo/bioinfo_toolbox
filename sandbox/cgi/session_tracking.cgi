#!/usr/bin/env perl

use strict;
use warnings;
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
use Apache::Session::File;

my $cgi=new CGI;

my $cookie=$cgi->cookie("myCookie");

my (%session,$id);
mkdir '/tmp/sessions';
eval {

tie %session,'Apache::Session::File',$cookie,{Directory=>'/tmp/sessions'};

}; #catch any error

if ($@)
{
    if ($@=~/^Object does not exist in the data store/)
    {
	#create a new one
	tie %session,'Apache::Session::File',undef,{Directory=>'/tmp/sessions'};
	$cookie=undef;
    } else
    {
	print $cgi->header('text/html','503 Service Unavailable');
	die "Error: $@ ($!)";
	exit;
    }
}

unless ($cookie)
{
    #retrieve the new session id from the %session
    $cookie=$cgi->cookie(-name=>"myCookie",-value=>$session{_session_id},-expires=>"+1d");
    print $cgi->header(-type=>"text/html",-cookie=>$cookie);
} 
else 
{
    print $cgi->header;

}

print $cgi->start_html("Session Demo"),
$cgi->h1("Helloe, you are session id ",$session{_session_id}),
$cgi->end_html;

untie %session;
