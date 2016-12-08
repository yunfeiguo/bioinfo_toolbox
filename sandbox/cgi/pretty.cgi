#!/usr/bin/env perl

use warnings;
use strict;
use CGI::Carp qw/fatalsToBrowser/;
use CGI::Pretty qw/:standard/;

$CGI::Pretty::LINEBREAK="\n\n";
push @CGI::Pretty::AS_IS,qw/LI B A/;
my $cgi=new CGI::Pretty;
print header,
      start_html("Pretty HTML Demo"),
      ol(li(['first','second','third'])),
end_html;
