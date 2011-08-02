package SFF;

=head1 SFF

SFF: a perl module for processing SFF files in the Metagenomic Pipeline. Currently very simple

=head1 Version

Version 2.0, Last Modified June 22, 2011

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
use Constants;
use MessagesFileHandling;

=head1 Methods

=over 5

=item $fnaFileName = extract_MIDfna_from_SFF($sffstring,$midNumber,$jobName);

Given a mid number (from the 454 MID kit) and 1 or more sfffiles, this
subroutine will extract the sequences matching that MID and convert to
a .fna fasta file.

=cut

sub extract_MIDfna_from_SFF{

  my($self,$sffString,$midNumber,$jobName,$outfile) = @_;

  # Make an sff file with just your mid
  my$sfffileCommand = "sfffile -o $jobName -s $midNumber -mcf ".Constants::RL_MID_CONFIG_PARSE()." $sffString";
  MessagesFileHandling->append_to_file("Sff Splitting:\t$sfffileCommand\n",$outfile);
  system ("$sfffileCommand");
  
  # RL and Standard have different mid naming conventions
  if($midNumber =~ /RL(\d+)/){
    system "sffinfo -seq $jobName\.RLMIDSRL$1\.sff >$jobName\.fna";
    system "gzip --best $jobName\.RLMIDSRL$1\.sff";
    MessagesFileHandling->append_to_file("Compressed sample sff:\t$jobName\.RLMIDSRL$1\.sff.gz\n",$outfile);
  }
  elsif($midNumber =~/MID(\d+)/){
    system "sffinfo -seq $jobName\.GSMIDSMID$1\.sff >$jobName\.fna";
    system "gzip --best $jobName\.GSMIDSMID$1\.sff";
    MessagesFileHandling->append_to_file("Compressed sample sff:\t$jobName\.GSMIDSMID$1\.sff.gz\n",$outfile);
  }
  else{
    die "Unexpected mid number $midNumber.  Expected form RL# or MID#.  Show this message to Brian.  $!";
  }

  MessagesFileHandling->append_to_file("Fasta for sff $sffString in file:\t$jobName\.fna\n",$outfile);
  return "$jobName\.fna";
}

=item my$sffFile = extract_from_sff(headerList=>$textFileWithFastaHeaders,sffString=>$sffString,letter=>'i'|'e',output_filename=>$string(optional))

(subroutine written by Sathish in June 2011)

This subroutine generates a new sfffile using a list of headers to include or exclude. The immediate application is to generate an sff file free of human sequences for deposition to public databases.

The following arguments are required in the invocation hash:
 headerList => text file of accession numbers. Each line of the text file should contain one fasta header, with or without a leading >
 sffString => SFF file(s) from which the subset should be withdrawn
 letter => 'i' or 'e'. Argument i means that the header list passed should be included in the final sff (e.g., these and only these headers will be returned). Argument e means that the header list passed should be excluded from the final sff (e.g., the returned sff will contain all of the reads in the starting sffString except for the reads in this header list)
 
The following arguments are optional. 
 output_filename => name of output SFF. If not explicitly specified, a name will be created based on input arguments

The return is the name of the sff file created by the subroutine.

=cut

sub extract_from_sff{  
  
  my($self,%arg) = @_;  

  unless($arg{headerList} and $arg{sffString} and $arg{letter}){
  croak("Must provide arguments \'headerList\', \'sffString\', \'letter\' to subroutine extract_from_sff in SFF package. Please read the documentation and try again.\n");
  }
  
  # croak if any passed options are not one of the expected possible values
  foreach (keys %arg){
    unless($_ =~ /^headerList$|^sffString$|^letter$|^output_filename$/){
    croak ("Unexpected argument $_ passed to extract_from_sff subroutine. Please read documentation and try again. Arguments are case sensitive.");
    }
  }
  
  unless(($arg{letter} eq "e") or ($arg{letter} eq "i")){
    croak("Argument 'letter' must have value 'e' or 'i' in subroutine extract_from_sff in SFF package. Please try again. Value is case sensitive.");
  }

  # Specifies default naming format for output, which will be overridden by the optional output_filename argument if supplied
  my$output_filename = "$arg{headerList}\_";
  $output_filename .= "excluded_from_"       if ($arg{letter} eq "e");
  $output_filename .= "included_in_"         if ($arg{letter} eq "i");
  $output_filename .= "\_$arg{sffString}";
 
  ($output_filename=$arg{output_filename}) if defined $arg{output_filename};

  # Make an sfffile with headers you want to keep or exclude, specified by i/e, from an input text file consisting of accession ids
 
  my$sfffile_command = "sfffile -$arg{letter} $arg{headerList} -o $output_filename $arg{sffString} ";
  system ("$sfffile_command");
  return "$output_filename";
}

1;

=back

=head1 Author

Brian Muegge and Sathish Subramanian

=head1 COPYRIGHT

Copyright (c) 2011, Washington University School of Medicine

=cut
