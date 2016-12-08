package Hash;

use strict;
use warnings;

sub perlhash
{
    my $hash=0;
    for (split //,shift)
    {
	$hash=$hash*33+ord($_);
    }
    return $hash;
}

sub TIEHASH
{
    my $h= {
	keys=>0,
	buckets => [ [],[],[],[],[],[],[] ],
	current => [ undef,undef]
    };
    return bless $h,shift;
}

#empty an existing hash
sub CLEAR
{
    my ($h)=@_;
    $h->{keys}=0;
    @{$h->{buckets}}=([],[],[],[],[],[],[]);
    @{$h->{current}}=(undef,undef);
}

#look up a specified key in a hash
sub lookup
{
    my ($h,$key)=@_;
    
    my $buckets=$h->{buckets};
    my $bucket=perlhash($key) % @$buckets;

    my $entries=@{ $buckets->[$bucket] };
    my $entry;
    if ($entries > 0)
    {
	#look for correct entry inside the bucket
	$entry=0;
	while ($buckets->[$bucket][$entry][0] ne $key)
	{
	    if (++$entry == $entries)
	    {
		#None of the entries in the bucket matched
		$bucket=$entry=undef;
		last;
	    }
	}
    }
    else
    {
	#the relevant bucket was empty
	$bucket=$entry=undef;
    }
	return ($bucket,$entry);

}

#check whether a key exists in a hash
sub EXISTS
{
    my ($h,$key)=@_;
    my ($bucket,$entry)=&lookup($h,$key);

    #if $bucket is undefined, the key doesn't exist
    return defined $bucket;
}

#delete a key-value pair from a hash
sub DELETE
{
    my ($h,$key)=@_;
    my $buckets=$h->{buckets};
    my ($bucket,$entry)=&lookup($h,$key);

    if (defined $bucket)
    {
	#remove the entry from the bucket, and return its value
	$entry=splice (@{$buckets->[$bucket]},$entry,1);
	return $entry->[1];
    }
    else
    {
	return undef;
    }
}

#store a key-value pair in a hash
sub STORE
{
    my ($h,$key,$val)=@_;

    my $buckets=$h->{buckets};
    my ($bucket,$entry)=&lookup($h,$key);

    if (defined $bucket)
    {
	$buckets->[$bucket][$entry][1]=$val;
    }
    else
    {
	$h->{keys}++;
	$bucket=&perlhash($key) % @$buckets;
	push @{$buckets->[$bucket]},[$key,$val];

	#expand the hash if all the buckets are full
	if ($h->{keys} > @$buckets)
	{
	    #double the number of buckets
	    my $newbuckets=[];
	    push (@$newbuckets,[]) for 1..2*@$buckets;

	    #redistribute keys
	    for my $entry(map {@$_} @$buckets)
	    {
		my $bucket=&perlhash($entry->[0]) % @$newbuckets;
		push @{$newbuckets->[$bucket]},$entry;
	    }
	    $h->{buckets}=$newbuckets;
	}
    }
}

sub FIRSTKEY
{
    my $h=shift;
    @{ $h->{current} } = (undef,undef);
    return $h->NEXTKEY(@_);
}

sub NEXTKEY
{
    my $h=shift;
    my $buckets=$h->{buckets};
    my $current=$h->{current};

    my ($bucket,$entry)=@{ $current};

    if (!defined $bucket || $entry+1 == @{$buckets->[$bucket] } )
    {
	FIND_NEXT_BUCKET:
	do {
	    if (++$current->[0] == @$buckets)
	    {
		@{$current}= (undef,undef);
		return undef;
	    }
	} while (@{ $buckets->[$current->[0]] } == 0);
	$current->[1]=0;
    } else
    {
	$current->[1]++;
    }
    return $buckets->[$current->[0]][$current->[1]][0];
}

1;
