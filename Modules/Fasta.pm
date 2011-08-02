package Fasta;

=head1 Fasta

Fasta: a perl class for parsing and writing fasta files, creating a collection of read objects

=head1 SYNOPSIS

use Fasta;

my$fasta = Fasta->new(file => 'FastaFile.fa');

=head1 Version

Version 2.0, Last modified June 22, 2011.

=cut

sub version{
  return "2.0";
}

sub date_last_modified{
  return "June 22, 2011";
}

#Begin Code

use strict;
use warnings;
use Carp;
use Read;

=head1 Object Attributes

Required:
  file (a DNA fasta file): get_file()
    
Created:
  reads, a Hash of header--> read object pairs: get_reads()
  readOrder, an array of headers in the order encountered: get_readOrder()

=cut

# Class data and methods
{
  # a list of all attributes with default values
  # Must start with a fasta file.
  # Each header/sequence pair will be used to create a read object
  # The ordered list of created objects will be stored in _reads attribute
  #   as an array. 
  my %_attributes = (
      _file          => '??',
      _reads         => {},
      _readOrder     => [],
   );       

  # Return a list (array) of all attributes
  sub _all_attributes {
    keys %_attributes;
  }
  
  # Return the value of an attribute
  sub _attribute_value {
    my($self,$attribute) = @_;
    $_attributes{$attribute}; 
  }
}

# Constructor
sub new{
  my($class,%arg) = @_;
  
  #Create a new object
  my$self = bless{},$class;

  # Must start with a fasta file
  unless($arg{file}){
    croak("No fasta file passed as argument 'file' in Fasta objection creation.");
  }
  
  # Set the fasta attribute
  $self->{_file} = $arg{file};

  # Check if fasta format. Die if not
  unless($self->is_dna_fasta()){
    croak("The fasta file is not expected fasta format (headers beginning with > and followed immediately by word; sequences consisting only of ACTGN)");
  }

  # Create the read collection
  $self->parse_fasta();

  return $self;
}

# Accessors 
sub get_file    {$_[0] -> {_file}} 
sub get_reads   {$_[0] -> {_reads}} # returns a hash
sub get_readOrder {$_[0] -> {_readOrder}} # returns the headers in order

# No mutators allowed

=head1 Additional Methods

=over 5

=item $number = $obj->read_number();

Returns the number of reads associated with the object.

=cut

sub read_number{
  my($self) = @_;
  return scalar @{$self->{_readOrder}};
}

=item $clonedFasta = $FastaObject->clone();

Creates a clone of every read in the object, and returns a new object pointing to those clones in _reads and _readOrder. Intended use is to operate on this clone (e.g. remove filtered reads) and not manipulate the original object.

Is there a better way to do this? Problem is that attribute _reads points directly to the read objects, so I think each of them has to be cloned as well

=cut

sub clone{
  my($caller) = @_;

  # Extract the class name
  my$class = ref($caller);
  
  # Create the new object
  my$self = bless{},$class;

  $self->{_file} = $caller->{_file};

  # for every read, clone a copy of the read, and add to this clone's readOrder and reads Hash
  # This is tricky: if you just copy _readOrder or _reads to $self,
  #   then modifications of the clone will still affect the original
  foreach my$read (@{$caller->{_readOrder}}){
    # clone the read from the caller hash, and put clone of it in self hash
    my$clonedRead = ${$caller->get_reads()}{$read}->clone();
    ${$self->{_reads}}{$read} = $clonedRead;
    push(@{$self->{_readOrder}},$read);
  }

  return $self;
}

=item $obj->remove_filtered_reads()

Operates directly on the Read objects associated with the fasta object, and removes all reads with Read->passFilters()==0. The intended use of this method is to remove reads from a cloned object, then use that cloned object for calculating average length or printing fasta files.

=cut

sub remove_filtered_reads{
  my($self) = @_;

  # Removing from hash with delete function is trivial
  # delete can work on arrays, but not in a useful way
  #   - delete removes the element but doesn't change the index numbering
  # splice function will delete the element and change the index numbering
  # But callings splice while iterating means indexes are changed before 
  #   iteration can access through all elements
  # E.g., if iterating from 1 to 10, calling splice on element 7 immediately
  #   changes what was index 8-->7, 9-->8, etc. Thus when the iteration
  #   advances to 8, the element that was originally 8 is skipped
  # So need to iterate in reverse order, e.g. 5000 to 1, so that splicing
  #  element 4997 (and making 4998-->4997, 4999-->4998, 5000-->4999) only
  #  changes the indices of elements already processed

  # From Max element to 0
  for (my$i=($self->read_number()-1);$i>=0;$i--) {
    #  readOrder[$i] gives us header, 
    #    and header keys the hash pointing to read objects
    if ($self->{_reads}{$self->{_readOrder}[$i]}->passFilters()==0) {
      delete $self->{_reads}{$self->{_readOrder}[$i]};
      # splice takes three arguments: the array to operate on, the element
      #  to start remove at, and the number of elements to remove
      splice (@{$self->{_readOrder}},$i,1);
    }
  }
  
  return $self;
}

=item $obj->remove_reads(%reads_to_remove)

Given a hash with headers of reads to remove as keys, the method will remove those read objects from the Fasta object. The intended use of this method is to remove reads from a cloned object, then use that cloned object for calculating average length or printing fasta files.

=cut

sub remove_reads{
  my($self,%reads_to_remove) = @_;

  # Removing from hash with delete function is trivial
  # delete can work on arrays, but not in a useful way
  #   - delete removes the element but doesn't change the index numbering
  # splice function will delete the element and change the index numbering
  # But callings splice while iterating means indexes are changed before 
  #   iteration can access through all elements
  # E.g., if iterating from 1 to 10, calling splice on element 7 immediately
  #   changes what was index 8-->7, 9-->8, etc. Thus when the iteration
  #   advances to 8, the element that was originally 8 is skipped
  # So need to iterate in reverse order, e.g. 5000 to 1, so that splicing
  #  element 4997 (and making 4998-->4997, 4999-->4998, 5000-->4999) only
  #  changes the indices of elements already processed  

  # From Max element to 0
  for (my$i=($self->read_number()-1);$i>=0;$i--) {
    #  readOrder[$i] gives us header, 
    #    and header keys the hash pointing to read objects
    if (defined $reads_to_remove{$self->{_readOrder}[$i]}) {
      delete $self->{_reads}{$self->{_readOrder}[$i]}; # delete from the hash
      # splice takes three arguments: the array to operate on, the element
      #  to start remove at, and the number of elements to remove
      splice (@{$self->{_readOrder}},$i,1);
    }
  }
  return $self;
}

=item $self->parse_fasta();

Method is called by the constructor. Parses the input fasta file and creates a read object from each encountered header/sequence pair. Stores the head in attribute _readOrder, and stores the header->object map in attribute _reads.

=cut

sub parse_fasta{
  my($self) = @_;

  # Initialize variables
  my ($sequence, $header) = "";

  #Read in the fasta file
  open(FASTA,$self->{_file}) or croak("Cannot open the fasta file $self->{_file}");

  # Grab each header/sequence and make a read. We will add to the self{_reads} array as we go
  while (my$line = <FASTA>) {
    chomp($line);

    # If line begins with letters, it is sequence
    if ($line =~ /^\w/) {
      $sequence .= uc($line);
    } elsif ($line =~ /^>/) { # else, lines beginning with > are headers

      # We see the first header before we see sequence
      # So store headers until all the sequence is seen, then create the read object
      # The following loop won't be triggered when we see the first header,
      #    since at the point $header = undef

      if (defined($header)) {
	push @{$self->{_readOrder}},$header;
	$self->{_reads}{$header} = Read->new(header=>$header,sequence=>$sequence);	 
      }
      
      # Keep only the first "word" of the header as key
      my@temp= split /\s+/,$line;
      $header=$temp[0];
      
      # Reset the sequence to empty for the next round
      $sequence = "";
    }
  }
  
  # Catches the last header and sequence that weren't caught in elsif loop
  push @{$self->{_readOrder}},$header;
  $self->{_reads}{$header} = Read->new(header=>$header,sequence=>$sequence);

  close FASTA;
}

=item $logical = $object->is_dna_fasta();

Returns true iff the multifasta meets some basic requirements for multifasta file formats with DNA sequence. Called from the constructor.

=cut

sub is_dna_fasta {
  my($self) = @_;

  open(IN,$self->{_file}) or die "$!";
  while (my$line=<IN>) {
    chomp($line);
    # note that I allow header lines to have multiple words,
    #   even though we won't keep all words in the read object
    # Should I add an attribute to read object to keep the entire
    #   header line? This may be helpful for parsing?
    unless(($line=~/^>\w+/)or($line=~/^[ACTGN]+$/i)){
      return 0;
    }
  }
  return 1;
}

=item $string = $object->put_fasta();

Write a (multi)fasta file to the return string. The string will include every sequence in the fasta object's read collection, in the order that the original Fasta file was parsed. Formatted with not more than 80 bases per line.

=cut

sub put_fasta {
  my($self) = @_;

  my$out = ""; # will hold the string of fasta formatted text
  my $lineLength = 80;

  foreach(@{$self->{_readOrder}}){
    $out.="$_\n";
    for (my$pos = 0;$pos< (${$self->{_reads}}{$_})->sequenceLength(); $pos += $lineLength){
      $out.= substr(${$self->{_reads}}{$_}->get_sequence(),$pos,$lineLength);
      $out.= "\n";
    }
  }

  return $out;
}

=item my($meanLength,$standardDeviation) = $obj->lengthMeanStdev();

Returns the mean and standard deviation of the sequences in the fasta object's read collection.

=cut

sub lengthMeanStdev{
  my($self) = @_;

  my($mean,$lengthSum,$deviation,$squareDeviation,$sumDeviation,$variance,$stdev) = "";

  foreach my$read (@{$self->get_readOrder()}){
    my$object = ${$self->get_reads()}{$read};
    $lengthSum += $object->sequenceLength(),
  }

  $mean = $lengthSum / scalar(@{$self->{_readOrder}});

  # calculate standard deviation
  foreach (@{$self->{_readOrder}}) {
    $deviation = (${$self->{_reads}}{$_}->sequenceLength()) - $mean;
    $squareDeviation = $deviation*$deviation;
    $sumDeviation += $squareDeviation;
  }
  $variance = $sumDeviation / (scalar(@{$self->{_readOrder}}) -1);
  $stdev = sqrt($variance);

  return($mean,$stdev);
}

=item  $text = $object->lengthHistogram();

This subroutine takes all of the sequences in the fasta object, bins them into a histogram by length, and returns a string containing the raw histogram data in ascending order. The string will be two columns: the first contains the sequence length and the second contains the number of reads at that length.

=cut 
    
sub make_fasta_histogram{
  my($self) = @_;

  # Initialize the Histogram hash and Return string
  my%Histogram = ();
  my$return = "SequenceLength\tCount\n";

  #  get length of each sequence and add it to the histogram
  foreach(@{$self->{_readOrder}}){
    my$sequenceLength = ${$self->{_reads}}{$_}->sequenceLength();
    if (defined $Histogram{$sequenceLength}) {
      $Histogram{$sequenceLength}++;
    } else {
      $Histogram{$sequenceLength} = 1;
    }
  }

  # sort numerically and append to return string
  foreach (sort {$a<=>$b} keys %Histogram) {
    $return.="$_\t$Histogram{$_}\n";
  }

  return $return;
}

=item $fastaObject->make_fasta_flat_file($jobName);

Subroutine will print to $jobName\_FastaFlatFile.txt a matrix with attributes for each read in the fasta object.

=cut

sub make_fasta_flat_file{
  my($self,$jobName) = @_;
  
  open(OUT,">$jobName\_FastaFlatFile.txt") or die "Could not open the file $jobName\_FastaFlatFile.txt: $!";
  
  print OUT "Header\tSequenceLength\tTooShort\tTooManyN\tReplicate\tHostGenomeHit\tPassFilter\tKeggBlastHit\tKeggAnnotatedHit\tMeropsHit\tCogBlastHit\tCogAnnotatedHit\tGutGenomeHit\n";

  foreach my$read (@{$self->get_readOrder()}){
    my$object = ${$self->get_reads()}{$read};
    print OUT $object->get_header(),"\t",$object->sequenceLength(),"\t",$object->tooShort(),"\t",$object->tooManyN(),"\t",$object->get_replicate(),"\t",$object->get_hostGenomeHit(),"\t",$object->passFilters(),"\t",$object->get_keggBlastHit(),"\t",$object->get_keggHit(),"\t",$object->get_meropsHit(),"\t",$object->get_cogBlastHit(),"\t",$object->get_cogHit(),"\t",$object->get_gutGenomeHit(),"\n";
  }
  close OUT;
  
  return "$jobName\_FastaFlatFile.txt";
}

=back

=head1 AUTHOR

   Brian Muegge

=head1 COPYRIGHT

Copyright (c) 2011, Washington University School of Medicine

=cut
