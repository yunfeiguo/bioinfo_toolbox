#!/usr/bin/env perl
#show falcon jobs progress in each folder separately

use strict;
use warnings;
use File::Spec;

my $usage = "Usage: $0 [current wd (by default) or falcon root dir containing 0-rawreads,1-preads_ovl,2-asm-falcon]\n";
die $usage unless @ARGV <= 1;
warn $usage if @ARGV == 0;
my $dir;
if(@ARGV == 0) {
    $dir = $ENV{PWD};
} else {
    $dir = shift @ARGV;
}
my @subdir = qw/0-rawreads 1-preads_ovl 2-asm-falcon/;
@subdir = map {File::Spec->catfile($dir,$_)} @subdir;
my $output;

#check dir
for my $i($dir,@subdir) {
    die "$i not found\n" unless -d $i;
}

for my $i(@subdir) {
    my $total = 0;
    my $donetotal = 0;
#count total jobs
    warn "Processing $i\n";
    warn "Counting total number of jobs, this may take a while\n";
    
    opendir (DIR, $i) or die "Error: cannot read from dir $i\n";
    my @alldir = readdir(DIR);
    my @founddir = grep { m/^m_/ || m/^job_/ } @alldir;
    $total += @founddir;
    closedir (DIR);
    
    #{my $subtotal = `find $i -name 'm_*' -type d -o -name 'job_*' -type d| wc -l`;
	#$subtotal =~ s/(\d+).*/$1/s; #use 's' to let . match \n
	#$total += $subtotal;
    #}

    #rdb_build x2, da_done x2, cns x1, p_merge x1, falcon_asm x1;
    #these jobs don't have a separate folder, and will be run every time
    #so we just add them to total number of jobs
    if($i =~ /0-rawreads/) {
	$total += 3;
	#count preads c_xxxx.sh jobs
	my $pread_dir = File::Spec->catfile($dir,"0-rawreads","preads");
	if (-d $pread_dir) {
		opendir (DIR, $pread_dir) or die "error: cannot read from dir $pread_dir\n";
		my @consensusjobs = grep {m/^c_.+\.sh$/} readdir (DIR);
		$total += @consensusjobs;
		#my $consensus_jobs = `find $pread_dir -name 'c_*.sh' -type f| wc -l`;
		#$consensus_jobs =~ s/(\d+).*/$1/s; #use 's' to let . match \n
		#$total += $consensus_jobs;
	}
    } elsif ($i =~ /1-preads/) {
	$total += 3;
    } elsif ($i =~ /2-asm/) {
	$total += 1;
    }

#count finished jobs
    warn "Counting jobs that are done in $i, this may take a while\n";
    
    for my $nextdir (@founddir) {
    	opendir (DIR, "$i/$nextdir") or die "Error: cannot read from dir $nextdir\n";
    	my @foundfile = grep { m/_done$/ } readdir (DIR);
    	$donetotal += @foundfile;
    }
    	
    
    #{my $subtotal = `find $i -name '*_done' -type f | wc -l`;
#	$subtotal =~ s/(\d+).*/$1/s; #use 's' to let . match \n
#	$donetotal += $subtotal;
    #}
    
    
#summary
    my $percentage = sprintf("%.2f",$donetotal/$total*100);
    if($donetotal == 0) {
    	#when there are 0 jobs done, we don't show total jobs as this number will change later.
	$output .= "NOTICE: 0 jobs done in $i\n";
    } else {
        $output .= "NOTICE: $total jobs in $i.\n";
	$output .= "NOTICE: $donetotal ($percentage %) jobs done in $i.\n";
    }
}

print $output;
