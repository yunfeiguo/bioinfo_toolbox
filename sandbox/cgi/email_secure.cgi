#!/usr/bin/perl

use warnings;
use strict;
use CGI::Pretty qw/:standard -debug/;
use CGI::Carp qw/fatalsToBrowser/;

$CGI::POST_MAX=1000; #limit size of POST
$ENV{'PATH'}="/bin";
print header,start_html("Mail yourself a greeting");
my $mail_to=param('email');

#check the email address is decent
if (not $mail_to or $mail_to !~ /\@/)
{
    print start_form(-method=>'POST'),
    h2("Please enter an email address"),
    p(textfield('email')),
    p(submit()),
    end_form;
} elsif ($mail_to =~ tr/;<>*|`&$!#[]{}:'"//)
{
    print h2("Invalid address");
} else
{
    if (open MAIL,"|mail $mail_to")
    {
	print MAIL "Hello from email!\n";
	close MAIL;
	print h1("Greeting sent");
    } else
    {
	print h2("Failed");
    }
}
print end_html;
