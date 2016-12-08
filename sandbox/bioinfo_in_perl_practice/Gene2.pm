package Gene2;
#
# second version of Gene.pm
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
		},$class;
	$class->_incr_count;
	return $self;
	}

#Accessor

sub get_name { $_[0]->{_name}}
sub get_organism {	$_[0]->{_organism}}
sub get_chromosome { $_[0]->{_chromosome}}
sub get_pdbref { $_[0]->{_pdbref}}

#Mutator

sub set_name {
	my ($self,$name) =@_;
	$self->{_name}=$name if $name;
}
sub set_organism {
	my ($self,$organism) = @_;
	$self->{_organism}=$organism if $organism;
}
sub set_chromosome {
	my ($self,$chromosome) = @_;
	$self->{_chromosome}=$chromosome if $chromosome;
}
sub set_pdbref {
	my ($self,$pdbref) = @_;
	$self->{_pdbref}=$pdbref if $pdbref;
}

1;
