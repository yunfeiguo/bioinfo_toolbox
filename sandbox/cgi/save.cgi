#!/usr/bin/env perl

use warnings;
use strict;

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
use Fcntl qw/:flock/; #for file locking symbols

my $msqfile="/tmp/state.msg";
my $cgi=new CGI;

print $cgi->header(),$cgi->start_html("Stateful CGI Demo");

if (open (LOAD,$msgfile))
{
	flock LOAD,LOCK_SH; #shared lock
		my $oldcgi=new CGI(LOAD);
	flock LOAD,LOCK_UN; #release lock
		close (LOAD);

	if (my $oldmsg=$oldcgi->param('message'))
	{
		print $cgi->p("the prev msg was: $oldmsg");
	}
}

if (my $newmsg=$cgi->param('message'))
{
	print $cgi->p("The current message is: $newmsg");
	if (open (SAVE,"> $msgfile"))
	{
		flock SAVE,LOCK_EX; #exclusive lock
			$cgi->save(SAVE);
		flock SAVE,LOCK_UN;
	}
	else
	{
		print $cgi->font({-color=>'red'},"failed to save: $!");
	}
}
print $cgi->p("Enter a new msg:");
print $cgi->startform(-method=>'GET'),
      $cgi->textfield('message'),
      $cgi->submit({-value=>'Enter'}),
      $cgi->endform();

print $cgi->end_html();
