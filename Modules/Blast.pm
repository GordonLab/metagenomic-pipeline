package Blast;


=head1 Blast

Blast: a Perl module for creating and parsing tab delimited blast files.

=head 1 Version

Version 2.0, Last modified July 18,2011.

=cut

sub version{
  return "2.0";
}

sub date_last_modified{
  return "July 18, 2011";
}

# Begin code

use strict;
use warnings;
use Carp;
use Constants;
use MessagesFileHandling;

=head1 Methods

=over 12

=item $binary = is_tab_blast(blastout.txt)

Simple check to insure that file is in tab delimited format (blast option -m 8 or -m 9). Returns 1 if it is expected tab delimited format, 0 otherwise.

This script runs fairly slowly with large blast returns. May be ways to optimize with more sophisicated algorithm than regex matching.

=cut

sub is_tab_blast{
  my($file) = @_;

  open(IN,"<$file") or carp ("Could not open file $file");
  
  while(my$line=<IN>){
    chomp($line);
    
    # -m 9 format has 4 header lines separating each query's return
    # BLASTX 2.2.21 [Jun-14-2009]
    # Query: FLTSILN01EJKEM
    # Database: /srv/cgs/data/bmuegge/TempKEGGdir/KEGG_v52_KOSeq
    # Fields: Query id, Subject id, % identity, alignment length, mismatches,
    #   gap openings, q. start, q. end, s. start, s. end, e-value, bit score
    
    # -m 8 and -m 9 have 12 returns for each blast line as specified in 
    #    Fields above

    # For some weird reason, the space separating e-value and bit score
    #    is a tab plus one space. So \d+\t doesn't match that combo.
    # And since I chomped the line, the bit-score doesn't have a trailing tab
    # I found this method is fairly slow. Nick Semenkovich helped me to realize that bounding regex terms in parentheses stores them in memory and may contribute to slow speed. So I took out all parentheses, which did speed things up about 2 fold (though still "slow"). Also add to account for the decimal point in the third column (percent ID) and 12th column (bit score), as well as the e in the evalue (11th column)

 # BIG speed-up pointed out by Nick Semenkovich. Perl regex .* matches any character except newline. Tab, \t, is included in match. Because of the * quantifier, this implements a "greedy algorithm" that matches the entire line, then steps backwards iteratively until the rest of the expression can evaluate true. So .*\t in a tab delimited file will basically first get all columns, then step back and do all columns minus one, then all columns minus two, etc. To improve speed, use the special quantifier-? combination (e.g. +?), which means match the preceeding one or more times until the following character is found. Here, .+?\t means match anything until a tab is found. This takes run time from 40 seconds to 5 seconds for a modestly large blast return.

    unless (($line=~/^# (BLASTX|Query:|Database:|Fields:)/)or($line=~/^.+?\t.+?\t[0-9\.]+\t(\d+\t){7}[0-9e\.\-]+\t\s?[0-9\.]+/)) {
      return 0;
    }  
  }
  close IN;

  return 1;
}

=item $blastOutFile = make_blast(file=>$fna,suffix=>'string',db=><path to blast db>,program=>'blastx,outfile=>$outfile)

Given a fasta file and blast parameters, will invoke split_run_check_combine.pl to generate the desired tab delimited blast job. The following arguments are required in the invocation hash:

file => fasta file to be blasted
suffix => string to append to the blast output
db => path to blast DB (ideally a defined constant)
program => blast program to use, e.g. blastx or blastn
outfile => a file to append messages to
 
The following arguments are optional. Default values come from the Constants package. Note below lists expected value/type at command line 

b=>'#" : blast option -b limiting number of returns.
e=>'#" : blast option -e with minimum e-value.
m=>'9' : blast option -m. Default value is 8.
n=>'#' : number of sequences for each small fasta.
v=>'#' : number of one line returns allowed, blast options -v

Lastly, there are several options that may be invoked.

Q=>'#' : specifies which translation table to us in blastx. Recommend table 11 for bacterial table
F=> 'string' : Changes default filtering parameters. Recommended to use a 'soft filtering' approach to preven query matches from getting broken into multiple small HSP's (high-scoreing seqment pairs) [e.g. http://ab.inf.uni-tuebingen.de/software/megan/how-to-use-blast/ ]. For blastn searches, use argument F=>'m D', and for blastx searches use F=>'m S'. Script will add the necessary double quotes when passing the argument to blast.
z=> '#' : changes effective size of database. Intended to be used when database is searched volume by volume and e-values need to be modified to reflect total size.

L=>'1' : submit job with long queue. Default is submission to short queue (less than one hour)
M=>'#" : extra memory requirement for jobs. Default is 1 GB memory.

=cut

sub make_blast{
  my($self,%arg) = @_;

  unless($arg{file} and $arg{suffix} and $arg{db} and $arg{program} and $arg{outfile}){
    croak("Must provide arguments \'file\', \'suffix\', \'db\', \'outfile\'and \'program\' to subroutine make_blast in Blast package. Please read the documentation and try again.");
  }

  # croak if any passed options are not one of the expected possible values
  foreach (keys %arg){
    unless($_ =~ /^file$|^suffix$|^db$|^program$|^outfile$|^b$|^e$|^m$|^n$|^F$|^Q$|^L$|^M$|^z$|^v$/){
      croak ("Unexpected argument $_ passed to make_blast subroutine. Please read documentation and try again. Arguments are case sensitive.");
    }
  }

  # Use passed options to override defaults, if provided
  my($b,$e,$m,$n,$v) = (Constants::BLAST_B(),Constants::BLAST_E(),Constants::BLAST_M(),Constants::NUM_SEQ(),Constants::BLAST_V());

  ($b=$arg{b}) if defined $arg{b};
  ($e=$arg{e}) if defined $arg{e};
  ($n=$arg{n}) if defined $arg{n};
  ($v=$arg{v}) if defined $arg{v};
  if(defined $arg{m}){
    if ($arg{m}==9){
      $m=$arg{m};
    }
    else {
      croak("Argument m has value $arg{m}. Expected value is 9. Default value is ".Constants::BLAST_M().". No other blast output formats are supported at this time.");
    }
  }

  # Check the format of F argument,if provided
  if(defined $arg{F}){
    unless(($arg{F}eq'm D') or ($arg{F} eq 'm S')){
      croak("Argument F must be either \'m D\' or \'m S\'. Please try again.");
    }
  }

  my$blastCommand = "perl ".Constants::SPLIT_RUN_CHECK()." -i $arg{file} -s $arg{suffix} -n $n -p 'blastall -p $arg{program} -d $arg{db} -i INCLUDE_INFILE -b $b -m $m -e $e -v $v";

  my$extra_args = "";
  $extra_args .= " -Q $arg{Q}" if defined ($arg{Q});
  $extra_args .= " -z $arg{z}" if defined ($arg{z});
  $extra_args .= " -F \"$arg{F}\"" if defined ($arg{F});

  # The following will evaluate to a single quote if $extra_args is empty
  $blastCommand .= "$extra_args'";

  $blastCommand .= " -L 1" if (defined $arg{'L'});
  $blastCommand .= " -M $arg{M}" if (defined $arg{'M'});

  MessagesFileHandling->append_to_file("Submitting blast job:\t$blastCommand\n",$arg{outfile});

  system("$blastCommand");

  my$blastOut = MessagesFileHandling::wait_for_file("$arg{file}\.$arg{suffix}");

  MessagesFileHandling->append_to_file("Blast results written to:\t$blastOut\n",$arg{outfile});
  return $blastOut;
}

=item %Blast_Results = parse_tab_blast_file(file=>$blastFile,fastaObj=>$fastaObject,outfile=$outfile)

Input is tab delimited blast output, fasta object for these reads, and an outfile. Options can be passed to override the blast parser parameters from the Constants module:

e=>'#"     : Minimum e-value. Default Constants::BLAST_E (1e-5)
id =>'#'   : Minimum blast percent identity. Default Constants::MIN_ID (50)
align=>'#' : Minimum percent of query length aligned. Default CONSTANTS::MIN_ALIGN (50)
score=>'#" : Minimum bit score. Default Constants::MIN_SCORE (50)

Return is a Hash, where each keys (queryID) points to an Array-of-Arrays containing all passed blast returns in array format. 

Intended use is that this hash of returns will be passed to the annotation parsers.

=cut

sub parse_tab_blast_file{
  my($self,%arg) = @_;

  # Reject if not tab delimited
  croak("The passed file $arg{file} is not expected tab delimited format") unless (is_tab_blast($arg{file}));

  # must pass the fasta object and the blast file
  croak("Required arguments are file,fastaObj, and outfile to method parse_tab_blast_file") unless ($arg{file} and $arg{fastaObj} and $arg{outfile});

  # Use passed options to override defaults, if provided
  my($minE,$minID,$minAlign,$minScore) = (Constants::BLAST_E(),Constants::BLAST_MINID(),Constants::BLAST_MINALIGN(),Constants::BLAST_MINSCORE());

  ($minE=$arg{e})         if $arg{e};
  ($minID=$arg{id})       if $arg{id};
  ($minAlign=$arg{align}) if $arg{align};
  ($minScore=$arg{score}) if $arg{score};

  # print the parsing parameters to record
  MessagesFileHandling->append_to_file("Blast file $arg{file} being parsed to include returns with e score <= $minE, percent identity >= $minID, percent align >= $minAlign, and score >= $minScore\n",$arg{outfile});

  # Initialize hashes for tracking the returns
  # This hash is to make sure we don't store hsp's split into two lines
  my%query_subject=();
  # Query is Hash of Array-of-Arrays: first key is queryID, second key is an AoA.
  #   The Each element of the AoA stores an array with the blast returns 
  my(%query) = ();

  open(IN,$arg{file}) or croak("Could not open the blast file $arg{file}");
  while (my$line=<IN>) {
    chomp($line);
    unless($line=~/^#/){
      my@results = split/\t/,$line;
      
      # Calculate the percent of aligned length usings the read object   
      # The syntax is a little complicated, details here:
      #  * Brackets bounding arg{fastaObj} show the scope of the variable named with $
      #      ${arg{fastaObj}}
      #  This is necessary to seperate the object name from the method arrow ->
      #  * The hash call is bounded with ${} to force hash interpretation
      #      ${$... ...reads()}
      #  * The hash key in our object includes the carrot, while the blast return
      #  does not include the carrot. Thus we add carrot back with '>'.$results[0]
      #  Blast returns 7 and 6 are query stop and query start, respectively

      my$percentAlign = 100*(abs($results[7]-$results[6])+1)/
           ${${arg{fastaObj}}->get_reads()}{'>'.$results[0]}->sequenceLength();

      # Keep the read if meets thresholds
      if ( $results[10]  <= $minE and
	   $results[2]   >= $minID and
	   $percentAlign >= $minAlign and
	   $results[11]  >= $minScore ) {
	
	# This will ignore any hsp's split into several lines, 
	#   where each has same query and subject ID 
	unless (exists $query_subject{$results[0]}{$results[1]}){
	  $query_subject{$results[0]}{$results[1]} = 1;
	   
	  # Store the results
	  push @{$query{$results[0]}}, [@results];  
	}
      }
    }
  }
  close IN;

  return %query;
}

=item $bestHit = best_blast_hit(@blastReturn)

Input should be a 2D array returned by parse_tab_blast_file, corresponding to the blast matches for a single query sequence. This is the value stored in the parse_tab_blast_file hash for a key (queryID). This 2D array represents all of the passed filter matches for a given read/query. Each 'row' (primary index) of the 2D array is a different return (note that multiple HSP's from a single subject should have been removed). The 12 elements in the row (secondary index) are the 12 items returned by blast tab delimited formats (-m 8 and -m 9).

The output is the subjectID for the highest scoring alignment. In the event that more than one subjectID had the exact same e-value and bit-score ties for highest scoreing status, the return is the list of all subjectID's that had the equivalent alignment score, separated by spaces.

=cut

sub best_blast_hit{
  # no need for $self argument here, because passed argument is a simple array and not a hash
  my(@blastReturn) = @_;

  # Set the scores to the zero element 
  my($bestE,$bestScore) = ($blastReturn[0][10],$blastReturn[0][11]);
  my($bestHit) = "$blastReturn[0][1]";
 
  # Check other elements. If better, replace previous score. If tie, add to present
  for (my$i=1;$i<scalar(@blastReturn);$i++) {
    if ($bestE == $blastReturn[$i][10] and $bestScore == $blastReturn[$i][11]) {
      # the first return is the best blast hit: store subjectID as value of header
      $bestHit .= " $blastReturn[$i][1]";
    }
    elsif ($bestE > $blastReturn[$i][10]){
      ($bestE,$bestScore) = ($blastReturn[$i][10],$blastReturn[$i][11]);
      ($bestHit) = "$blastReturn[$i][1]";
    }
  }
  return $bestHit;
}

=item $bestHit = best_annotated_blast_hit(\@blastReturn,\%LookupHash)

There are two inputs to this subroutine, each of which needs to be passed dereferenced.

1) The first input should be a 2D array returned by parse_tab_blast_file, corresponding to the blast matches for a single query sequence. This is the value stored in the parse_tab_blast_file hash for a key (queryID).This 2D array represents all of the passed filter matches for a given read/query. Each 'row' (primary index) of the 2D array is a different return (note that multiple HSP's from a single subject should have been removed). The 12 elements in the row (secondary index) are the 12 items returned by blast tab delimited formats (-m 8 and -m 9).

2) The second input is a lookup hash, where keys are subjectID's from the blast database that are "annotated" (for instance, have a KEGG or COG annotation). The value stored in this hash is not used by the subroutine, but should probably contain the appropriate annotation.

The output of the subroutine is the subjectID for the highest scoring alignment for which the subjectID has a definition in the lookup hash. If the highest scoring alignment in the array does not have a defined value in the lookup hash, it is ignored and the subroutine moves to the next alignment. This proceeds until an alignment is found with a definition in the lookup hash. In the event that more than one subjectID have the exact same e-value and bit-score ties as the highest scoreing status, and the subjectID's have annotations in the lookup hash, the return is the list of all subjectID's that had the equivalent alignment score, separated by spaces.

=cut

sub best_annotated_blast_hit{
  my($blastReturn,$Lookup) = @_;
  # Dereference
  my@blastReturn = @$blastReturn;
  my%Lookup = %$Lookup;

  # Set the bestE and bestScore to zero
  my($bestE,$bestScore) = (0,0);
  my($bestHit) = "";
  
  # Iterate through every return: if it has a defined lookup, process
  for (my$i=0;$i<scalar(@blastReturn);$i++) {  
    if(defined $Lookup{$blastReturn[0][1]}){
      # Set the E, score, and hit if this is the first result with annotation
      if($bestE == 0 and $bestScore == 0){
	($bestE,$bestScore) = ($blastReturn[$i][10],$blastReturn[$i][11]);
	$bestHit = "$blastReturn[$i][1]";
      }
      # E and Score are already defined. Add this return if equal
      elsif ($bestE == $blastReturn[$i][10] and $bestScore == $blastReturn[$i][11]) {
	# the first return is the best blast hit: store subjectID as value of header
	$bestHit .= " $blastReturn[$i][1]";
      }
      elsif ($bestE > $blastReturn[$i][10]){
	($bestE,$bestScore) = ($blastReturn[$i][10],$blastReturn[$i][11]);
	($bestHit) = "$blastReturn[$i][1]";
      }
    }
  }

  return "$bestHit";
}

=item $bestHit = best_annotated_hit_pipeInHeader(@blastReturn)

This is a bit derivative of the other best_blast_hit methods, may be a way to make this more general

I've created a modified KEGG, COG, and MEROPS database where each header ends with |.... if the entry is annotated with a KO, COG, etc. There could be multiple |... patterns at the end of the text

I checked the databases before I created the header modifications, and to the best of my knowledge, no subject ID has a |

Rather than load the very large gene=>annotations lookup hashes with each return, I'll just keep blast hits that end with the specified pattern

There is only one input to this subroutine, the 2D array returned by parse_tab_blast_file, corresponding to the blast matches for a single query sequence. This is the value stored in the parse_tab_blast_file hash for a key (queryID).This 2D array represents all of the passed filter matches for a given read/query. Each 'row' (primary index) of the 2D array is a different return (note that multiple HSP's from a single subject should have been removed). The 12 elements in the row (secondary index) are the 12 items returned by blast tab delimited formats (-m 8 and -m 9).

The output of the subroutine is the subjectID for the highest scoring alignment for which the subjectID has an annotation subroutine. If the highest scoring alignment in the array does not have a defined annotation value, it is ignored and the subroutine moves to the next alignment. This proceeds until an alignment is found with a definition. In the event that more than one subjectID have the exact same e-value and bit-score ties as the highest scoreing status, and the subjectID's have annotations, the return is the list of all subjectID's that had the equivalent alignment score, separated by spaces.

=cut

sub best_annotated_hit_pipeInHeader{
  # no need for self as this is a passed array
  my(@blastReturn) = @_;
 
  # Set the bestE and bestScore to zero
  my($bestE,$bestScore) = (0,0);
  my($bestHit) = "";
  
  # Iterate through every return: if it has a defined lookup, process
  for (my$i=0;$i<scalar(@blastReturn);$i++) {  
    if($blastReturn[$i][1] =~/\|/){
      # Set the E, score, and hit if this is the first result with annotation
      if($bestE == 0 and $bestScore == 0){
	($bestE,$bestScore) = ($blastReturn[$i][10],$blastReturn[$i][11]);
	$bestHit = "$blastReturn[$i][1]";
      }
      # E and Score are already defined. Add this return if equal
      elsif ($bestE == $blastReturn[$i][10] and $bestScore == $blastReturn[$i][11]) {
	# the first return is the best blast hit: store subjectID as value of header
	$bestHit .= " $blastReturn[$i][1]";
      }
      elsif ($bestE > $blastReturn[$i][10]){
	($bestE,$bestScore) = ($blastReturn[$i][10],$blastReturn[$i][11]);
	($bestHit) = "$blastReturn[$i][1]";
      }
    }
  }

  return "$bestHit";
}

=item $aggregateBlastResults = make_multiple_blast(db=>$DB,suffix=>$suffix,output=>$string,other arguments for make_blast...)

This subroutine is intended to generate separate blast jobs for every volume of a multi-volume blast database. The benefit of this approach is that searches against single volumes (properly formatted) can be run on nodes with only 1 GB of memory and thus take advantage of more nodes on the cluster.

The required arguments to the subroutine are db, suffix, and output. Other arguments can be passed that will be passed unmodified to the make_blast method.

At the end of this method, the output files made by blast against each small volume will be aggregated into a single file, and this file name will be the return value of the subroutine. The small files will be deleted.

Highly recommended to pass option -z #### to correct the effective database length to reflect the size of the entire database. For more information, see the Gordon Lab wiki page "Managing big blast databases efficiently."

=cut

sub make_multiple_blast{
  my($self,%arg) = @_;
  
  # For this subroutine, I need db, suffix, and output.
  # Other arguments are passed to make_blast and will cause errors there if incorrect
  unless($arg{db} and $arg{suffix} and $arg{output}){
    croak("Must provide a value for 'db', 'output', and 'suffix' to make_multiple_blast subroutine argument hash in Blast package. Please read documentation and try again.\n");
  }
  
  # Capture the starting db and suffix, which will get modified later
  my$starting_db = $arg{db};
  my$starting_suffix = $arg{suffix};

  # Capture the output name, and remove it from the hash as it should not be passed to make_blast
  my$output = $arg{output};
  delete $arg{output};

  # initialize the list that will hold the temporary blast output files
  my@returnedBlastFiles = ();

  # Find the collection of volumes that will be searched.
  my@volumes = glob("$starting_db*\.psq $starting_db*\.nsq");
  
  # Exit if there are no volumes to search
  if(scalar(@volumes) == 0){
    croak("No BLAST volumes were found in subroutine make_multiple_blast using database $arg{db} (package Blast.pm). Please try again.\n");
  }
  # Else search each volume
  else {
    my$count = 0;
    foreach(@volumes){      
      # remove the trailing protein or nucleotide distinction to get down to the prefix.## stem
      $_ =~s/\.psq$//g;
      $_ =~s/\.nsq$//g;

      # Update the argument list for db and suffix, and pass this to make_blast routine
      $arg{db} = $_;
      $arg{suffix} = "$starting_suffix\_$count";
      my$smallFile = Blast->make_blast(%arg);
      push(@returnedBlastFiles,$smallFile);
      $count++;
    }
  }

  # We now have a number of files with portions of the blast return.
  # aggregate them to a single file with catfiles method, and remove the small fiels
  my$BigFile = MessagesFileHandling->catfiles($output,\@returnedBlastFiles);
  unlink(@returnedBlastFiles);
  return $BigFile;
}

=item $smallBlastResults = best_return_per_query(blastFile=>$hostBlast,jobName=>$jobName,suffix=$string);

This method was written to deal with the very large blast results being returned by host blasting. The problem is that blasting against these genomes can lead to a very large number of returned lines when a read has homology to a repetitive element in the host genome.

This script will read in the blast results returns, find only the best blast result for each query (based on e-value), and print out to a new file where for each queryID there is one and only one return (=best return). In the event that there are multiple returns with equivalent e-values, the first one encountered will be used.

The method will delete the passed file to save space.

=cut 

sub best_return_per_query{
  my($self,%arg) = @_;

  unless($arg{blastFile} and $arg{jobName} and $arg{suffix}){
    die "Need arguments 'blastFile' and 'jobName' to subroutine best_return_per_query. $!";
  }
  
  my%Query_BestE = ();
  my%Query_Line = ();

  open(IN, "<$arg{blastFile}") or die "Could not open the file $arg{blastFile} : $!";
  while(my$line=<IN>){
    chomp($line);
    unless($line=~/^#/){
      my@temp = split/\t/,$line;
      if(not defined $Query_Line{$temp[0]}){
	$Query_Line{$temp[0]} = $line;
	$Query_BestE{$temp[0]} = $temp[10];
      }
      else { # already defined
	if($temp[10] < $Query_BestE{$temp[0]}){
	  $Query_Line{$temp[0]} = $line;
	  $Query_BestE{$temp[0]} = $temp[10];
	}
      }
    }
  }
  close IN;

  open(OUT,">$arg{jobName}\_$arg{suffix}") or die "Could not open the file $arg{jobName}\_$arg{suffix} : $!";
  foreach my$query (sort keys %Query_Line){
    print OUT "$Query_Line{$query}\n";
  }
  close OUT;

  # remove the large starting blast file
  unlink($arg{blastFile});

  return("$arg{jobName}\_$arg{suffix}");
}

1;

=back

=head1 AUTHOR

Brian Muegge, Nick Semenkovich

=head1 COPYRIGHT

Copyright (c) 2011, Washington University School of Medicine

=cut
