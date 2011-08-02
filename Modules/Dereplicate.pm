package Dereplicate;

=head1 Dereplicate

Dereplicate: a subroutine to call the Dereplicator program, parse it's output, and assign replicate status to read objects as attributes.

=head1 Version

Version 2.0; Last Modified June 22, 2011

=cut;

sub version{
  return "2.0";
}

sub last_modified{
  return "June 22, 2011";
}

# Begin code

use strict;
use warnings;
use Carp;
use MessagesFileHandling;
use Constants;

=head1 Methods

=over 5

=item my$replicate_Summary_file = Dereplicate->dereplicate(fastaObj=>$FastaObject,jobName=>$jobName,outfile=>$outfile}

The subroutine initiates a call to the dereplicator script, and waits for the script to finish. Once the output is written, the subroutine parses the output and identifies which of the passed reads have been flagged as replicates. It changes the _replicate attribute of those reads to "1". The return is a reference to the text file that contains the summary from the deplicator code. All other output from the code is deleted. The parameters for running the code can be altered in the Constants file.

=cut

sub dereplicate{

  my($self,%arg) = @_;

  croak("Subroutine dereplicate requires hash arguments fastaObj, jobName, and outfile\n") unless (defined $arg{fastaObj} and defined $arg{jobName} and defined $arg{outfile});

  my$command = "python ".Constants::EXTRACT_REPLICATES()." ".$arg{fastaObj}->get_file()." ".Constants::REPLICATE_PERCENT()." ".Constants::REPLICATE_LENGTH()." ".Constants::REPLICATE_START()." $arg{jobName}\_Dereplicate";
  
  MessagesFileHandling->append_to_file("Dereplicator call:\t$command\n",$arg{outfile});

  # write this as a job file and qsub the job
  system("echo \"$command\" >$arg{jobName}\_DereplicateJob.txt");

  system("qsub -l h_vmem=2G $arg{jobName}\_DereplicateJob.txt");

  my$derepSummary = MessagesFileHandling::wait_for_file("$arg{jobName}\_Dereplicate/extracted_clusters.cluster_summary");

  # GoodHeaders will have definitions for reads that are NOT replicates, and also for the single sequence (longest) that is used to represent replicate clusters
  my%GoodHeaders = process_cluster_summary($derepSummary);

  foreach(@{$arg{fastaObj}->get_readOrder()}){
    unless(defined $GoodHeaders{$_}){
      # A value of 1 in attribute _replicate marks the read object as a replicate which should be removed
      ${$arg{fastaObj}->get_reads()}{$_}->set_replicate(1);
    }
  }
  
  # Move the summary file to main directory, then remove all other files
  
  system("mv $derepSummary $arg{jobName}\_dereplicate_summary.txt");

  # remove all remaining files in tmp and main
  unlink glob "$arg{jobName}\_Dereplicate/tmp/*";
  rmdir "$arg{jobName}\_Dereplicate/tmp";
  unlink glob "$arg{jobName}\_Dereplicate/*";
  rmdir "$arg{jobName}\_Dereplicate";
 
  return "$arg{jobName}\_dereplicate_summary.txt";
}

=item %ReplicateHash = process_cluster_summary($summary_file);

Intended as an internal method only. Given the summary file from the dereplicator script, this parser will return a hash where headers for every non-replicate, or the single header representing the representative sequence for a cluster, are returned as the key. Thus headers that are replicated within the cluster are not found in the hash. Used by the main script to determine which reads get marked as replicate.

=cut

sub process_cluster_summary{
  my($summary) = @_;

  my%return = ();
  open(IN,"<$summary") or die "Could not open the file $summary:$!";
  while(my$line=<IN>){
    chomp($line);
    # EXPECTED FORMAT:
    #Cluster	Ref sequence	Num of seq
    #1	FLTSILN01EJKEM	9
    #2	FLTSILN01BIMFS	1
    #3	FLTSILN02MOW0U	8
    if($line=~/^(\d+)\t(.*)\t(\d+)$/){
      # Add carrot so it matches the Fasta object headers
      $return{">".$2} = $3;
    }
  }
  close IN;
  return %return;
}

1;

=back

=head1 AUTHOR

Brian Muegge

=head1 COPYRIGHT

Copyright (c) 2011, Washington University School of Medicine

=cut
