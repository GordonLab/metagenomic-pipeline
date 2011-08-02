package Read;

=head1 Read

Read: a perl class for creating Read objects for metagenomic annotation.

=head1 SYNOPSIS

use Read;

my$read1 = Read->new(
   header => '>FLTSILN01EJKEM';
   sequence => 'TAATTGATGGTGGGTCAACTACCGTATATCCCAAAGA';
);

=head1 Version

Version 2.0, Last modified June 22, 2011.

=cut

sub version{
  return "2.0";
}

sub date_last_modified{
  return "June 22, 2011";
}

# Begin code

use strict;
use warnings;
use Carp;
use Constants; # a package in Modules, containing constants

=head1 Object Attributes: Accesors(default value), Mutators(expected values)

Required (no mutation allowed):
  header: get_header(????)
  sequence: get_sequence(????)

Optional:
  hostGenomeHit: get_hostGenomeHit(????), set_hostGenomeHit(genomeHitID)
  keggHit: get_keggHit(????), set_keggHit(keggHitID) : This is intended to be the best hit from the database with a KO annotation
  keggBlastHit: get_keggBlastHit(????), set_keggBlastHit(keggHitID) : This is intended to be the best hit in the database, regardless of whether or not it has a KO annotation
  cogHit: get_cogHit(????), set_cogHit(cogHitID)
  cogBlastHit: get_cogBlastHit(????), set_cogBlastHit(cogBlasthitID)
  replicate: get_replicate(0), set_replicate(1)
  gutGenomeHit: get_gutGenomeHit(????), set_gutGenomeHit(gutHitID)

=cut

# Class data and methods
{
  # a list of all attribtues for the object with default values
  my %_attributes = (
      _header          => '????',   
      _sequence	       => '????',
      _hostGenomeHit   => '????',
      _replicate       => '0', # this is 1 if it needs to be removed
      _keggHit         => '????',
      _keggBlastHit    => '????',
      _meropsHit       => '????',
      _cogHit          => '????',
      _cogBlastHit     => '????',
      _gutGenomeHit    => '????',
  );
  
  # Return a list (array) of all attributes
  # @attributes = $read->_all_attributes();
  sub _all_attributes {
    keys %_attributes;
  }

  # Return the value of an attribute [pass attribute as argument]
  # This is really an object internal method and shouldn't be accessed
  #   from calling program.
  # From outside, better to use get_header(), etc.
  sub _attribute_value {
    my($self,$attribute) = @_;
    $_attributes{$attribute}; 
  }
}

# The constructor method
# Called from class, e.g. $read1 <- Read->new();
sub new{
  my($class,%arg) = @_;
  # Create a new empty object (a hash) with bless function
  my$self = bless {},$class;

  # The following line calls the _all_attributes subroutine to get all possible attributes defined above
  foreach my $attribute ($self->_all_attributes()) {
    # E.g., attribute = "_header", argument="header",...
    # this is sort of a tricky substitution (=~s) which doesn't require
    #   you to set argument = attribute in a previous step.
    # return is the contents of the parentheses, (.*)
    my($argument) = ($attribute =~/^_(.*)/); # removes the leading _

    # Reject the call if we don't have a sequence and header defined
    unless($arg{header} and $arg{sequence}){
      croak("Arguments \"header\" and \"sequence\" are required to create a read object.");
    }

    # reject headers that aren't fasta format
    # the extra clause 'and $arg{$argument}' insures that we enter 
    #   this loop iff a header value was passed.
    if ($argument eq "header" and $arg{$argument}) {
      unless ($arg{$argument}=~/^>.*/) {
	croak "$argument attribute \'$arg{$argument}\' is not expected fasta header format";
      }
    }

    # Reject sequences that aren't DNA (with N's)
    if ($argument eq "sequence" and $arg{$argument}) {
      if ($arg{$argument}=~/[^ACTGN]/i) {
	croak "$argument attribute \'$arg{$argument}\' is not DNA sequence (contains elements other than A,C,T,G,N";
      }
    }
  
    # If the attribute had an argument passed, add it to object
    if (exists $arg{$argument}) {
      $self->{$attribute} = $arg{$argument};
    }
    # else set attributes to default
    else {
      $self->{$attribute} = $self->_attribute_value($attribute);
    }
  }

  return $self;
}

# ACCESSORS
# e.g, $header -> get_header();
# two arguments, the object ($_[0]) and the attribute to get
sub get_header        {$_[0] -> {_header}} 
sub get_sequence      {$_[0] -> {_sequence}}
sub get_hostGenomeHit {$_[0] -> {_hostGenomeHit}}
sub get_keggHit       {$_[0] -> {_keggHit}}
sub get_keggBlastHit  {$_[0] -> {_keggBlastHit}}
sub get_meropsHit     {$_[0] -> {_meropsHit}}
sub get_cogHit        {$_[0] -> {_cogHit}}
sub get_cogBlastHit   {$_[0] -> {_cogBlastHit}}
sub get_replicate     {$_[0] -> {_replicate}}
sub get_gutGenomeHit  {$_[0] -> {_gutGenomeHit}}

# MUTATORS
# header and sequence can't be mutated

sub set_hostGenomeHit {
  my($self,$value) = @_;
  $self->{_hostGenomeHit} =$value if $value; # will only set if value is supplied  
}
sub set_keggHit {
  my($self,$value) = @_;
  $self->{_keggHit} =$value if $value; # will only set if value is supplied
}
sub set_keggBlastHit{
  my($self,$value) = @_;
  $self->{_keggBlastHit} =$value if $value; # will only set if value is supplied
}
sub set_meropsHit {
  my($self,$value) = @_;
  $self->{_meropsHit} =$value if $value; # will only set if value is supplied
}
sub set_cogHit {
  my($self,$value) = @_;
  $self->{_cogHit} = $value if $value; # will only set if value is supplied
}
sub set_cogBlastHit{
  my($self,$value) = @_;
  $self->{_cogBlastHit} =$value if $value; # will only set if value is supplied
}
sub set_replicate {
  my($self,$value) = @_;
  $self->{_replicate} =$value if $value; # will only set if value is supplied
}
sub set_gutGenomeHit{
  my($self,$value) = @_;
  $self->{_gutGenomeHit} =$value if $value; # will only set if value is supplied
}

=head1 Additional Methods

=over 5

=item $clonedRead = $read->clone();

Returns an exact copy of the calling object (read). Intended use is to make a copy that can be manipulated, e.g. change attributes in copy without changing the attribute in the original object.

=cut

sub clone{
  my($caller) = @_;

  # Extract the class name
  my$class = ref($caller);
  
  # Create the new object
  my$self = bless{},$class;

  # Make except copy of each attribute
  foreach my$attribute ($self->_all_attributes()){
    $self->{$attribute} = $caller->{$attribute};
  }

  return $self;
}

=item $length = $object->sequenceLength();
  
Returns the length of the sequence stored in the object.

=cut

sub sequenceLength{
  my($self) = @_;  # operate on self object
  length($self->{_sequence});
}

=item $binary = $object->tooManyN();

Returns 1 if the sequence contains more than 2 N's, or 2 adjacent N's.

=cut

sub tooManyN{ 
  my($self)=@_;
  if ( ($self->{_sequence}=~s/N/N/ig >2) or ($self->{_sequence}=~/N{2,}/)){
    return 1;
  }
  else{
    return 0;
  }
}

=item $binary = $object->tooShort();

Returns true if the sequence length is less than constant MINIMUM_SEQUENCE_LENGTH in the Constants.pm module

=cut

sub tooShort{ # set isFinal==0 if pruneShort is true
  my($self)=@_;
  if ($self->sequenceLength()< Constants::MINIMUM_SEQUENCE_LENGTH()){
    return 1;
  }
  else {
    return 0;
  }
}

=item $binary = $object->passFilters();

Returns 0 if the sequence is too short, has too many N's, has a match in the host genome, or is marked as a replicate. Else returns 1

=cut

sub passFilters{
  my($self)=@_;
  if($self->tooShort() or $self->tooManyN() or $self->get_replicate() or ($self->get_hostGenomeHit() ne "????")){
    return 0;
  }
  else {
    return 1;
  }
}

# sub sequenceNucleotideCounts {
#   my($self) = @_;
#   my%Hash = ('A'=>($self->{_sequence}=~ s/A/A/ig),'C'=>($self->{_sequence}=~ s/C/C/ig),'G'=>($self->{_sequence}=~ s/G/G/ig),'T'=>($self->{_sequence}=~ s/T/T/ig),'N'=>($self->{_sequence}=~ s/N/N/ig),'Other'=>($self->{_sequence}=~ s/[^ACTGN]/X/ig)); # note the ^ in Other invocation matches anything except what follows, e.g anything not ACTGN
#   return %Hash;
# }

1;

=back

=head1 AUTHOR

Brian Muegge, motivated by Chapter 3 of Tisdale "Mastering Perl for Bioinformatics"

=head1 COPYRIGHT

Copyright (c) 2011, Washington University School of Medicine

=cut
 
