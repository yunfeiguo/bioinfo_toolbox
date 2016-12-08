#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Std;
use File::Find qw/find/;
use File::Spec;

my %options;
getopts("vrnd",\%options);
die "Usage: $0 [-vrnd] <dir> [regex pattern]\n".
"unlinks regular files in a directory and subdiretories if required, optionally with a certain regex pattern".
"-v verbose\n".
"-r recursive, still only removes files but no directories\n".
"-d add a delay to reduce burden on file system\n".
"-n dry-run\n"
unless @ARGV == 1 or @ARGV == 2;
warn "remove regular files only\n";

my $dir = shift @ARGV;
my $pattern = shift @ARGV;
my $filecount = 0;
my $skipcount = 0;

if($options{r}) {
    find(
	{
	    wanted          => sub{&rm($File::Find::name)},
	    no_chdir        => 1,
	    bydepth		=> 1,
	},$dir);
} else {
    opendir my $dh,$dir or die "failed to open dir: $!\n";
    while(readdir $dh)
    {
	my $target = File::Spec->catfile($dir,$_); 
	&rm($target);
    }

    closedir $dh;
}
warn "All done.\n$filecount removed,$skipcount skipped\n";
sub rm {
    my $target = shift;
    sleep 1 if $options{d} and time % 5 == 0;
    if(-f $target and $target =~ /$pattern/) {
	if($options{n}? (print $target,"\n"):(unlink $target)) {
	    warn "$target removed\n" if $options{v};
	    $filecount++;
	} else {
	    die "failed to remove $target $!\n";
	}
    } else {
	warn "$target not a file, skipping ...\n" if $options{v};
	$skipcount++;
    }
}
