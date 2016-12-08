#!/usr/bin/env perl

use warnings;
use strict;
use lib "/opt/perl/lib/site_perl/5.14.2/";
use URI::Escape;


print "Content-Type: text/html\n\n";

my $mail_to=$ENV{'QUERY_STRING'};
$mail_to=uri_unescape($mail_to); #decode special char
print "<HTML><HEAD><TITLE>Mail yourself a greeting</TITLE>";
print "</HEAD><BODY><H1>Greeting Sent!</H1>";
print "</BODY></HTML>";

open (my $MAIL,"|mail $mail_to");
print $MAIL "Hello from Email!\n";
close $MAIL;

