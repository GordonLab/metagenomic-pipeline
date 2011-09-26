#!/usr/bin/env perl

use warnings;
use strict;

# Metagenomic_Pipeline_2.0.pl
# Written by Brian Muegge, 2011
# Last Modified: July 2011
# Purpose: Wrapper for 454 Metagenomic annotation, linking with various modules

# TODO : change all "dies" in the subroutines of packages to "croak" - die tells which line in the package called the error, while croak
#  tells which line in the calling script triggered the error, which is much more useful for debugging.

sub version{ return "2.0";}

use Getopt::Long;   # To capture options from command line
use Cwd;            # Gets current working directory, which can be optionally overridden with new directory for output

# Invoke the various modules and packages that will be used
use lib './Modules'; # Packages and Object classes are stored in this directory
use Annotation; # contains counting routines for KEGG and COG summarization
use Blast; # Generate and parse blast jobs, and assign annotation
use Constants; # contains constant parameters
use Dereplicate;
use Fasta; # Fasta object class. Used to create and manage read objects frome very metagenomic read.
use MessagesFileHandling; # Package needed for listener scripts, outfile management, etc
use SFF; # Package needed to split sff files by mid sequence

system('./.git-status-tester.sh');

# Initialize variables for input and options
my($sffString,$fna,$name,$midNumber) = ""; # will contain strings
my$directory = Cwd->cwd(); # If an argument is passed to option -dir, the output will be written to that directory instead of cwd
my($kegg,$cog,$gutGenome,$merops,$host,$dereplicate,$histogram,$sequence,$noFilter) = "";  # flags with default value FALSE

# Capture the options passed so we can print them for our records
my$commandline = join(" ",@ARGV);

# Process the options passed
&get_options();

# Print usage message if the passed options don't generate an adequate job command
MessagesFileHandling->usage() unless ((($sffString && $midNumber && $name) or ($fna && $name))and not($sffString && $fna));

# Create an outfile and print the invocation record. 
my($jobName,$outfile) = MessagesFileHandling->start_outfile($directory,$name,$commandline);

# 1) Split sff, if neccessary (if sff and mid passed). Create a fasta file and compress the starting sff
($fna = SFF->extract_MIDfna_from_SFF($sffString,$midNumber,$jobName,$outfile)) if $sffString;

# 2) Create a fasta object from the sequences
# A Read object will be created for every sequence in fasta, and fasta object has methods to access all read objects.
my$FastaObject = Fasta->new(file=>$fna);
MessagesFileHandling->append_to_file("Number of reads in starting fasta file $fna:\t".$FastaObject->read_number()."\n",$outfile);

# 3) Pruning: hostGenome and Dereplication
if($host){
  # In the future, add an option for people who want to blast against more than one host genome?
  # note that this generates a HUGE blast output file, because a single read can hit multiple repetitive elements in the genome
  # This was causing the program to crash in parsing because too much memory was being used
  # To limit this, I'm passing b=1 and v=1 (because any high quality hit is enough, I don't need multiple)
  # I've also created a new method that will read in the resulting blast file, find only the best hit per queryID,
  #   (based on evalue) and print that file for parsing (and also delete the first large blast file). The parse 
  #   thresholds will still be enacted in the parse_tab_blast_file method
  my$hostBlast = Blast->make_blast(file=>$fna,suffix=>'host',db=>Constants::HOSTDB(),program=>'blastn',outfile=>$outfile,M=>'2',F=>'m D',b=>'1',v=>'1');
  
  # Limit the file to only the best return per query ID
  my$smallHostBlast = Blast->best_return_per_query(blastFile=>$hostBlast,jobName=>$jobName,suffix=>'smallHostBlast.txt');
  MessagesFileHandling->append_to_file("The best return from Host BLAST for each queryID was extracted and written to $smallHostBlast.\n",$outfile);
  
  #parse. For human contaminant removal, require id to be 75% instead of normal 50%
  my%hostBlastHash = Blast->parse_tab_blast_file(file=>$smallHostBlast,fastaObj=>$FastaObject,outfile=>$outfile,id=>'75');

  system ("gzip --best $smallHostBlast");
  # update object
  foreach my$readID (keys %hostBlastHash){
    ${$FastaObject->get_reads()}{">".$readID}-> set_hostGenomeHit(Blast::best_blast_hit(@{$hostBlastHash{$readID}}));
  }
  # remove the blast hash to free up memory
  undef %hostBlastHash;
}

if($dereplicate){
  # Dereplicate: gets rid of known pyrosequencing artifacts from artifical replication in emulsion PCR
  # replicate reads will be noted by setting _replicate attribute == 1 in the read object
  my$deRepSummaryFile = Dereplicate->dereplicate(fastaObj=>$FastaObject,jobName=>$jobName,outfile=>$outfile);
  MessagesFileHandling->append_to_file("Dereplicator summary file:\t$deRepSummaryFile\n",$outfile);
}

# 4) Functional Annotation: cogg, kegg, merops, gutgenomes, other?
if($kegg){

# Use the following line if you want to search by volume, e.g. against v58_modHeader
#  my$keggBlast = Blast->make_multiple_blast(file=>$fna,suffix=>'kegg',db=>Constants::KEGGDB(),output=>"$jobName\_keggBlastOut.txt",program=>'blastx',outfile=>$outfile,F=>'m S',z=>Constants::KEGGDBLENGTH(),Q=>'11',b=>'10',v=>'10');

# Use the following to search against the annotation only database, which can load on a single 2 GB node
  my$keggBlast = Blast->make_blast(file=>$fna,suffix=>'kegg',db=>Constants::KEGGDB(),program=>'blastx',outfile=>$outfile,F=>'m S',Q=>'11',z=>Constants::KEGGDBLENGTH(),M=>'2');

  my%keggBlastHash = Blast->parse_tab_blast_file(file=>$keggBlast,fastaObj=>$FastaObject,outfile=>$outfile);
  system("gzip --best $keggBlast");

  # Update the read objects
  foreach my$readID (keys %keggBlastHash){
# Use the following line i.f.f you need to store best hits that aren't annotated
#    ${$FastaObject->get_reads()}{">".$readID}->set_keggBlastHit(Blast::best_blast_hit(@{$keggBlastHash{$readID}}));
    ${$FastaObject->get_reads()}{">".$readID}->set_keggHit(Blast::best_annotated_hit_pipeInHeader(@{$keggBlastHash{$readID}}));   
  }
  undef %keggBlastHash;
}

if($cog){

# Use the following to blast by volume, e.g against STRING/9.0_modHeader
#  my$cogBlast = Blast->make_multiple_blast(file=>$fna,suffix=>'cog',db=>Constants::COGDB(),output=>"$jobName\_cogBlastOut.txt",program=>'blastx',outfile=>$outfile,F=>'m S',z=>Constants::COGDBLENGTH(),Q=>'11',b=>'10',v=>'10');

# Use the following to search against the annotation only database, which can load on a single 2 GB node
  my$cogBlast = Blast->make_blast(file=>$fna,suffix=>'cog',db=>Constants::COGDB(),program=>'blastx',outfile=>$outfile,F=>'m S',Q=>'11',z=>Constants::COGDBLENGTH(),M=>'2');
  my%cogBlastHash = Blast->parse_tab_blast_file(file=>$cogBlast,fastaObj=>$FastaObject,outfile=>$outfile);
  system("gzip --best $cogBlast");

  # Update the read objects
  foreach my$readID (keys %cogBlastHash){
#    ${$FastaObject->get_reads()}{">".$readID}->set_cogBlastHit(Blast::best_blast_hit(@{$cogBlastHash{$readID}}));
    ${$FastaObject->get_reads()}{">".$readID}->set_cogHit(Blast::best_annotated_hit_pipeInHeader(@{$cogBlastHash{$readID}}));   
  }
  undef %cogBlastHash;
}

if($merops){
  my$meropsBlast = Blast->make_blast(file=>$fna,suffix=>'merops',db=>Constants::MEROPSDB(),program=>'blastx',outfile=>$outfile,F=>'m S',Q=>'11');
  my%meropsBlastHash = Blast->parse_tab_blast_file(file=>$meropsBlast,fastaObj=>$FastaObject,outfile=>$outfile);
  system("gzip --best $meropsBlast");
  foreach my$readID (keys %meropsBlastHash){
    # All reads in merops database are annotated, so no distinction between annotated reads and not annotated reads
    ${$FastaObject->get_reads()}{">".$readID}->set_meropsHit(Blast::best_blast_hit(@{$meropsBlastHash{$readID}}));
  }
  undef %meropsBlastHash;
}

# 5) Taxonomy. 
# TODO - add PhymmBl?
if($gutGenome){
  # using the blast parameters from Arumugam Nature 2011
  my$gutGenomeBlast = Blast->make_blast(file=>$fna,suffix=>'gutGenomeBlastOut.txt',db=>Constants::GUTDB(),program=>'blastn',outfile=>$outfile,e=>'1e-20',F=>'m D',b=>'10',n=>'2000');
  my%gutBlastHash = Blast->parse_tab_blast_file(file=>$gutGenomeBlast,fastaObj=>$FastaObject,outfile=>$outfile,e=>'1e-20',align=>'80');
  system("gzip --best $gutGenomeBlast");
  foreach my$readID (keys %gutBlastHash){
    ${$FastaObject->get_reads()}{">".$readID}->set_gutGenomeHit(Blast::best_blast_hit(@{$gutBlastHash{$readID}}));
  }
  undef %gutBlastHash;
}

# 6) Summary and outfiles (including histograms and stats)

# note the key thing here is to operate either on all reads or just "high quality" reads (in a cloned object) as needed

# First, make a flat file summarizing the entire project
# Be sure that if major modules are added (e.g MEROPS, PhymmBL) that the corresponding hard coded line in make_fasta_flat_file gets updated
my$FlatFile = $FastaObject->make_fasta_flat_file($jobName);
system("gzip --best $FlatFile");
MessagesFileHandling->append_to_file("Flat file with all per read records in:\t$FlatFile\.gz\n",$outfile);

# Someone may want to pass a fasta file directly with no filtering, e.g. no removal of short reads
# If so, they should pass $noFilter == 1 as option
unless($noFilter){
  $FastaObject->remove_filtered_reads();
  MessagesFileHandling->append_to_file("Number of high quality reads after applying filters:\t".$FastaObject->read_number()."\n",$outfile);
}

my($meanLength,$stdev) = $FastaObject->lengthMeanStdev();
MessagesFileHandling->append_to_file("The mean length of high quality sequences is:\t$meanLength\n",$outfile);
MessagesFileHandling->append_to_file("The standard deviation of high quality sequences is:\t".sprintf("%.2f",$stdev)."\n",$outfile);

if($histogram){
  # Make a histogram of high quality sequence length
  open (OUT,">$jobName\_Histogram.txt") or die "$!";
  print OUT $FastaObject->make_fasta_histogram();
  close OUT;
  MessagesFileHandling->append_to_file("Histogram in file:\t$jobName\_Histogram.txt\n",$outfile);
}

# Write a fasta of quality sequences
if($sequence){
  open(OUT,">$jobName\_qualitySeqs.fa") or die "$!";
  print OUT $FastaObject->put_fasta();
  close OUT;

  # Compress
  system("gzip --best $jobName\_qualitySeqs.fa");

  MessagesFileHandling->append_to_file("Passed filter, high quality sequences in file:\t$jobName\_qualitySeqs.fa.gz\n",$outfile);
}

# Compress the fna file
system("gzip --best $fna");

# Summarize the kegg annotation, if invoked
if($kegg){
  my($KOcount,$ECcount,$PathwayCount) = Annotation->get_kegg_summary(fastaObj=>$FastaObject,jobName=>$jobName,name=>$name);
  MessagesFileHandling->append_to_file("KO count file:\t$KOcount\nEC count file:\t$ECcount\nPathway count file:\t$PathwayCount\n",$outfile);
}
# summarize the gut genome hits, if invoked
if($gutGenome){
  my($GenomeCount,$GenusCount) = Annotation->get_gut_genome_summary(fastaObj=>$FastaObject,jobName=>$jobName,name=>$name);
  MessagesFileHandling->append_to_file("Bacterial genome count file:\t$GenomeCount\nBacterial species count file:\t$GenusCount\n",$outfile);
}
if($merops){
  my$MEROPSsummary = Annotation->get_merops_summary(fastaObj=>$FastaObject,jobName=>$jobName,name=>$name);
  MessagesFileHandling->append_to_file("Merops summary file:\t$MEROPSsummary\n",$outfile);
}
if($cog){
  my($COGcount,$CategoryCount) = Annotation->get_cog_summary(fastaObj=>$FastaObject,jobName=>$jobName,name=>$name);
  MessagesFileHandling->append_to_file("COG count summary:\t$COGcount\nCOG Category count summary:\t$CategoryCount\n",$outfile);
}

###############
# Subroutines
###############

sub get_options{
  GetOptions ( "sff=s" => \$sffString,
	       "fna=s" => \$fna,
	       "mid=s" => \$midNumber,
	       "name=s" => \$name,
	       "directory=s"=> \$directory,
	       "kegg=i" => \$kegg,
	       "cog=i" => \$cog,
	       "gut=i" => \$gutGenome,
	       "host=i" => \$host,
	       "merops=i" => \$merops,
	       "histogram=i" => \$histogram,
	       "dereplicate=i" => \$dereplicate,
	       "sequence=i" => \$sequence,
	       "nofilter=i" =>\$noFilter
	       );
  return 1;
}
