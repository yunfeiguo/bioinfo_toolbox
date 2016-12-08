package Gene;
#
# Fourth and final version of Gene.pm
#
use Carp;

#Class data methods
{
	my $_count=0;
	sub get_count { $_count};
	sub _incr_count { ++$_count};
	sub _decr_count { --$_count};

	my %_attribute_property= (	
		_name => ["????","read.required"],
		_organism => ["????","read.write.required"],
		_chromosome => ["????","read.write"],
		_pdbref => ["????","read.write"],
		_author => ["????","read.write"],
		_date => ["????","read.write"],
	);

	#Class helper methods
	sub _all_attributes {
		return keys %_attribute_property;
	}
	sub _default_attribute {
		my($self,$attribute)=@_;
		return $_attribute_property{$attribute}[0];
	}

	sub _permission {
		my ($self,$attribute,$permission) =@_;
		return $_attribute_property{$attribute}[1]=~/$permission/;
	}
}
#Constructor

sub new {
	my ($class, %arg) = @_;
	my $self=bless {},$class;
	print $self->_all_attributes,"\n";
	for $attribute($self->_all_attributes) {
		# regexp will only return matched values under array context!!!!
	my ($explicit_attribute)=($attribute =~ /_(.*)/);
		if ($arg{$explicit_attribute}) {
			$self->{$attribute}=$arg{$explicit_attribute};
		} elsif ($self->_permission($attribute,'required')) {
			croak "Attribute $explicit_attribute must be specified!",ref($self);
		} else {
			$self->{$attribute}=$self->_default_attribute($attribute);
		}
	}
	$class->_incr_count;
	return $self;
}

sub clone {
	my ($self,%arg) = @_;
	my $class=ref $self;
	my %arg_clone;
	for $attribute($self->_all_attributes) {
		my ($explicit_attribute)=($attribute=~/_(.*)/);
		$arg{$explicit_attribute} ? ($arg_clone{$explicit_attribute}=$arg{$explicit_attribute}):($arg_clone{$explicit_attribute}=$self->{$attribute}); 
	}
	my $self_clone=$class->new(%arg_clone);
	return $self_clone;
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
		croak "Attribute $attribute does not exist.",ref($self);
	}
	#print "\$operation:$operation,\$attribute:$attribute\n";
	no strict 'refs';
	if ($operation eq 'get') {
		if ($self->_permission($attribute,'read') ) {
			*{$AUTOLOAD} = sub {shift->{$attribute}};
			return $self->{$attribute};
		} else {
			croak "Action $operation denied for attribute $attribute",ref $self;
		}
	} elsif ($operation eq 'set') {
		if ($self->_permission($attribute,'write')) {
			*{$AUTOLOAD} = sub {shift->{$attribute}=shift};
			$self->{$attribute}=$new_value;
		} else {
			croak "Action $operation denied for attribute $attribute",ref $self;
		}

	}
	use strict 'refs';
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
	#pay attention to the context!!!
	my ($self)=@_;
	$self->_decr_count;
}

1;

=head1 Gene

Gene: objects for Genes with a minimum set of attributes

=head1 Synopsis

use Gene;

if (1) {
	my $my_gene=Gene->new(	
		name => "MYC",
		organism => "Homo Sapiens",
		pdbref => "PDB0001",
	);
	print $my_gene->get_name,"\n";
	print $my_gene->get_organism,"\n";
	print $my_gene->get_chromosome,"\n";
	print $my_gene->get_pdbref,"\n";
	#test AUTOLOAD failure
	print $my_gene->get_exon,"\n";

	$my_gene->set_name("C-MYC");
	$my_gene->set_organism("Baboon");
	$my_gene->set_chromosome("20q");
	$my_gene->set_pdbref("RS110");

	print $my_gene->get_name,"\n";
	print $my_gene->get_organism,"\n";
	print $my_gene->get_chromosome,"\n";
	print $my_gene->get_pdbref,"\n";
	print "Citation:",$my_gene->citation('Yunfei','2012'),"\n";
	my $your_gene=$my_gene->clone(	
		organism => "Mouse",
		name => "Aging",
		chromosome => "9",
	);
	print $your_gene->get_name,"\n";
	print $your_gene->get_organism,"\n";
	print $your_gene->get_chromosome,"\n";
	print $your_gene->get_pdbref,"\n";

	print "total count:",Gene3->get_count,"\n";
}
print "total count outside scope: ",Gene3->get_count,"\n";

=head1 AUTHOR

Unknown

=head1 COPYRIGHT

Copyright (c) 2012, Owe you, Inc.

=cut
