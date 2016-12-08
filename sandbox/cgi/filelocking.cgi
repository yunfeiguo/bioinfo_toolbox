#!/usr/bin/env perl

use strict;
use warnings;

use CGI::Pretty;
use CGI::Carp qw/fatalsToBrowser/;

use Fcntl qw/:flock :DEFAULT/; #for file locking symbols

my $msgfile="/tmp/state.msg";
my $cgi=new CGI::Pretty;

if ($cgi->param('message') eq 'secret')
{
    print $cgi->redirect("http://biocluster2.med.usc.edu/~yunfeiguo/cgi-bin/form.cgi");
} 
else
{
    print $cgi->header,$cgi->start_html("Stateful CGI Demo");

    if (open (my $LOAD,$msgfile))
    {
	flock ($LOAD,LOCK_SH | LOCK_NB) or die "Cannot read\n"; #shared lock
	my $oldcgi=new CGI::Pretty($LOAD);
	flock $LOAD,LOCK_UN; #release lock
	close $LOAD;

	if (my $oldmsg=$oldcgi->param('message'))
	{
	    print $cgi->p("The last msg was: $oldmsg");
	}

    }


    if (my $newmsg=$cgi->param('message'))
    {
	print $cgi->p("The current message is: $newmsg");
	if (open (my $SAVE,">",$msgfile))
	{
	    flock ($SAVE,LOCK_EX | LOCK_NB) or die "Failed to lock\n"; #exclusive lock
	    $cgi->save($SAVE);
	    flock $SAVE,LOCK_UN; #release lock
	} else
	{
	    print $cgi->font({-color=>'red'},"Failed to savei: $!");
	}
    }

    print $cgi->p("Enter a new message:");
    print $cgi->startform(-method=>'GET'),
    $cgi->textfield('message'), #auto-filled from CGI parameter if sent
    $cgi->submit({-value=>'Enter'}),
    $cgi->endform;

    print $cgi->end_html();
}
