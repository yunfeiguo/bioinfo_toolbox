#!/usr/bin/perl -T

use strict;
use warnings;
use CGI::Pretty;
use CGI::Carp qw/fatalsToBrowser/;

my $q=new CGI::Pretty;
print $q->header;
print $q->start_html(-title=>"Bob's Auto Parts Order");
print $q->start_form(-action=>"../php/processorder.php",-method=>"post");
print $q->table(
    {-border=>0},
    $q->Tr(
	{-bgcolor=>"blue"},
	$q->td({-width=>150},"Item"),
	$q->td({-width=>15},"Quantity"),
    ),
    $q->Tr( 
	$q->td({-align=>"center"},"Tires"),
	$q->td({-align=>"center",-size=>3,-maxlength=>3},$q->textfield("tireqty")),
	),
    $q->Tr(
	$q->td({-align=>"center"},"Oil"),
	$q->td({-align=>"center",-size=>3,-maxlength=>3},$q->textfield("oilqty")),
	),
    $q->Tr(
	$q->td({-align=>"center"},"Spark Plugs"),
	$q->td({-align=>"center",-size=>3,-maxlength=>3},$q->textfield("sparkqty")),
	),
    $q->Tr(
	$q->td({-colspan=>2 -align=>"center"},$q->submit(-value=>"Submit Order")),
    ),
);
print $q->end_form(),$q->end_html();

