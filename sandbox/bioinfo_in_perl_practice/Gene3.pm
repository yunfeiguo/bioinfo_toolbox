package Gene3;
#
# third version of Gene.pm
#
use Carp;

#Class data methods
{
	my $_count=0;
	sub get_count { $_count}
	sub _incr_count { ++$_count}
	sub _decr_count { --$_count}
}

#Constructor

sub new {
	my ($class, %arg) = @_;
	my $self={};
	$self =bless {	_name => $arg{name} || croak ("no name"),
			_organism => $arg{organism} || croak ("no no organism"),
			_chromosome => $arg{chromosome} || "????",
			_pdbref => $arg{pdbref} || "????",
			_author => $arg{author} || "unknown",
			_date => $arg{date} || "????",
		},$class;
	$class->_incr_count;
	return $self;
	}
#AUTOLOAD subroutine will take place of the previous accessor and mutator methods.
our $AUTOLOAD; #this statement is valid after perl5.6

sub AUTOLOAD {
	my ($self,$new_value) = @_;
	#print "\$AUTOLOAD: $AUTOLOAD\n";
	my ($operation,$attribute) = ($AUTOLOAD =~ /.*::(get|set)(_\w+)$/);
	unless ($operation && $attribute) {
		croak ("$AUTOLOAD method is not in recognized form (get|set)_attribute.");
	}
	unless (exists $self->{$attribute}) {
		croak "Attribute $attribute does not exist.",ref $self;
	}
	#print "\$operation:$operation,\$attribute:$attribute\n";
	no strict 'refs';
	if ($operation eq 'get') {
		*{$AUTOLOAD} = sub {shift->{$attribute}};
	} elsif ($operation eq 'set') {
		*{$AUTOLOAD} = sub {shift->{$attribute}=shift};
		$self->{$attribute}=$new_value;
	}
	use strict 'refs';
	return $self->{$attribute};
}

sub citation {
	my($self,$author,$date)=@_;
	$self->set_author($author) if $author;
	$self->set_date($date) if $date;
	return ($self->{_author},$self->{_date});
}
##Accessor
#
#sub get_name { $_[0]->{_name}}
#sub get_organism {	$_[0]->{_organism}}
#sub get_chromosome { $_[0]->{_chromosome}}
#sub get_pdbref { $_[0]->{_pdbref}}
#
##Mutator
#
#sub set_name {
#	my ($self,$name) =@_;
#	$self->{_name}=$name if $name;
#}
#sub set_organism {
#	my ($self,$organism) = @_;
#	$self->{_organism}=$organism if $organism;
#}
#sub set_chromosome {
#	my ($self,$chromosome) = @_;
#	$self->{_chromosome}=$chromosome if $chromosome;
#}
#sub set_pdbref {
#	my ($self,$pdbref) = @_;
#	$self->{_pdbref}=$pdbref if $pdbref;
#}

sub DESTROY {
	my ($self)=@_;
	$self->_decr_count;
}

1;
