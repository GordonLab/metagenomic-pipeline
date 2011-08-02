package MessagesFileHandling;

=head1 MessagesFileHandling

MessagesFileHandling: a Perl module for managing messages and waiting for files

=head1 Version

Version 2.0, Last modified June 22,2011.

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
use Constants; # a package in Modules, containing constants

=head1 Methods

=over 5

=item usage()

Prints a usage statement for the script, and exits the code. Need to add example invocations (after all coding is done, maybe Sathish could do this as part of a Tutorial?) 

=cut

sub usage{
  print "\n\tUsage: perl Metagenomic_Pipeline_2.0.pl\n".
    "\n\t   Required arguments\n".
    "\t   **************************************\n".
    "\t   -fna <Fasta File> OR -sff <SFFfile(s) surrounded by single quotes> -mid <mid number to extract>\n".
    "\t   -name <prefix for output files>\n".
    "\t   Your mid must be of form \"MID#\" for Titanium Standard or \"RL#\" for Rapid Library\n".
    "\n\t   Optional arguments\n".
    "\t   **************************************\n".
    "\t   -host 1  : prune out host sequences. Hostdb can be specied in Constants.pm\n".
    "\t   -kegg 1 : Submit sequences for kegg analysis pipeline\n".
    "\t   -cog 1 : Submit sequences for cog analysis pipeline. COG annotation is only ok as of June 2011\n".
    "\t   -gut 1 : Blast sequences against human gut genomes database. Database specified in constants\n".
    "\t   -merops 1: Blast sequences against MEROPS protease subunit database.\n".
    "\t   -dereplicate 1: Invoked the dereplicator script to get rid of pyrosequencer replicates.\n".
    "\t   -directory 'name' : Specify a directory to store output. Default current working directory.\n".
    "\t   -histogram 1: will generate a sequence length histogram file of quality reads\n".
    "\t      Also prints to runSummary the mean length of the sequences in the fasta\n".
    "\t   -sequence 1: Prints a fasta file of only those sequences that passed filters\n".
    "\n\t   Comment on memory requirements:\n".
    "\t   **************************************\n".
    "\t   The COG and KEGG annotation methods require large modules to be loaded. Each module\n".
    "\t   requires approximately 750 MB of memory. If you need to use both annotations, be sure\n".
    "\t   to add an extra 1GB to the memory request at job submission.\n".
    "\n\t   Examples:\n".
    "\t   **************************************\n".
    "\t    TODO TODO TODO \n".
    "\n";
  exit(1);
}

=item $newName = wait_for_file($fileName);

This sets up a loop that periodically checks for the existence of a file (often the final concatenated blast output from split_run_check_combine.pl).  Once the file is found, another loop to triggered to insure that the file size isn't changing. Both the WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS are parameters set in the Constants.pm package.

Once the file is found to have been created, and is no longer changing file size, the subroutine returns the name of the discovered file. This may be different than the passed $fileName argument, as the $fileName argument can include regular expression arguments, and the returned name will be the actual file name.

=cut

sub wait_for_file{
  my($fileName) = @_;

  # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

  # This prints to STDOUT
  print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

  # Check for the finished file every WAIT_TIME minutes
  until (-e "$fileName") {
    print "$fileName does not exist yet at ".localtime()."\n";
    sleep Constants::WAIT_TIME();
  }

  print "File $fileName found at ".localtime()."\n";

  # make sure the file is no longer being modified (changing size)
  my$time = 0;
  my$size = -s $fileName;

  until ($time == Constants::NO_CHANGE_IN_FILE_FOR_X_SECONDS()) {
    if ($size == -s $fileName) {
      $time++;
      sleep 1;
      if ($time%5 == 0) {
	print "No change in file size for $time seconds\n";
      }
    } else {
      $time = 0;
      $size = -s $fileName;
      print "The file size changed, sleeping for 1 minute\t";
      sleep 60;
      print"Waking up, try again\n";
    }
  }
  print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
  return $fileName;
}

# Maybe TODO: add a subroutine that looks for recheck jobs from split_run_check_combine
#  If they exist, resubmit it with doubled memory and long queue?
# Not quite sure how to do this - talk to J? Would require capturing the c#### code from his script

=item $Date_and_Time = get_Date_and_Time()

When invoked, this subroutine will return a string giving the date and time of the system in a "readable" format. The current format is Hour:Minute on Month Date, Year, e.g 12:05 on May 25, 2011.

=cut 

sub get_Date_and_Time{
  # Get current time and format for unique directory name                       
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime time;

  # sample return of call 'localtime time'                                      
  # 10 31 9 23 8 109 3 265 1                                                    

  # Convert the month number to a word                                          
  my@Month = ("January","February","March","April","May","June","July","August","September","October","November","December");
  return "$hour:".sprintf("%02d",$min)." on $Month[$mon] $mday, ".($year+1900);
}

=item ($jobName,$outfile) = start_outfile($outdirectory,$prefix,$commandLineArguments)

Will create or append to a file $prefix\_runSummary.txt in the $outdirectory. This "outfile" will be the record of arguments passed to the script, as well as the commands called and files created by methods in the various packages. I will probably also use this outfile to store information such as number of reads with kegg annotation, etc.

=cut 

sub start_outfile{

  my($self,$directory,$name,$commandline) = @_;

  # Create a file for the run summary. Remove trailing / from directory, if passed.  Join the directory to name for subsequent invocations with jobName
  $directory =~ s/\/$//;
  my$jobName = "$directory"."/"."$name";  # return this newly created name
  my$out = "$jobName\_runSummary.txt";

  # See if the directory already exists.  if not, create it.
  if(!(-d "$directory")){
    mkdir "$directory", 0755 or die "Can't create directory $directory : $!";
  }

  # modify the command line so it can be printed to tab delimited format more easily
  $commandline=~s/^\-//; # get rid of hyphen before first option passed
  $commandline=~s/\-/\t/g; # replace remaining hyphens (separating options) with tabs

  # Print a summary of invocation
  open (OUT, ">>$out") or die "Could not open the file $out: $!";
  print OUT "Metagenomic_Pipeline_2.0.pl\n";
  print OUT "Version:\t".main::version()."\n";
  print OUT "Script invoked at:\t".get_Date_and_Time()."\n";
  print OUT "Command Line Arguments:\t$commandline\n";
  #print each command on it's own line for later parsing
  my@options = split/\t/,$commandline;
  foreach my$option (sort @options){
    my@temp = split/\s/,$option;
    my$name = shift @temp;
    my$string = join(' ',@temp);
    print OUT "Option $name:\t$string\n";
  }

  close OUT;

  return ($jobName,$out);
}

=item append_to_file($message,$file)

Opens the filename passed (by appending to preexisting file), and prints the message passed to that file.  Then the script closes the file. Nothing is returned

=cut

sub append_to_file{

  my($self,$message,$file) = @_;

  open (OUT, ">>$file") or die "Could not open the file $file: $!";
  print OUT "$message";
  close OUT;

  return;
}

=item $concatenatedFile = catfiles($out_file_name,\@fileNamesToCat)

A simple subroutine with two arguments. First is the desired filename for the output. The second argument is a dereferenced array, where each element of the array is the name of a file that should be concatenated. The return is the name of the newly created file, which stores all of the lines in the passed array files, in the order passed. This subroutine does not remove the array files, though it is adviced that user delete those files after this call to avoid duplication of space.

=cut

sub catfiles{
  my($self,$filename,$files) = @_;
  # dereference
  my@files = @{$files};

  # Open the output to store the concatenation results
  open(OUT,">$filename") or die "$!";
  # for every file in the passed array, write each line to output
  foreach my$text (@files){
    open(IN,"<$text") or die "$!";
    while(my$line = <IN>){
      print OUT $line; 
    }
    close IN;
  }
  close OUT;
  
  return $filename;
}

1;

=back

=head1 AUTHOR

Brian Muegge

=head1 COPYRIGHT

Copyright (c) 2011, Washington University School of Medicine

=cut
