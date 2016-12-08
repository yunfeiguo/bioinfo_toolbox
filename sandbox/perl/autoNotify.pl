#!/usr/bin/env perl
#run a command, send output to email

use strict;
use warnings;

my $email="guoyunfei1989\@gmail.com";
my $cmd = "/home/yunfeiguo/bin/show_falcon_progress /home/yunfeiguo/projects/PacBio_reference_genome/falcon_aln/hx1_20150716";
#my $cmd = "date";
my $tmp = "/tmp/$$.stdout";
!system("$cmd > $tmp") or die "$cmd: $!\n";
!system("cat $tmp | /bin/mail -s '$cmd' $email") or die "mail: $!\n";
