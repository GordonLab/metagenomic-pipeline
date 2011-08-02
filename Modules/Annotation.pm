package Annotation;


=head1 Annotation

Annotation: a perl module for database/schema specific methods

=head 1 Version

Version 2.0, Last modified May 11,2011.

=cut

sub version{
  return "2.0";
}

sub date_last_modified{
  return "May 11, 2011";
}

# Begin code

use strict;
use warnings;
use Carp;
use Constants;
use Fasta;

=head1 Methods

=over 12

=item ($COGcount,$CategoryCount) = get_cog_summary(fastaObj=>$FastaObject,$jobName=>$jobName,name=>$name);

Returns the summary of COG hits at the level of COGID's and COG categories

=cut

sub get_cog_summary{
  my($self,%arg) = @_;
  unless($arg{fastaObj} and $arg{jobName} and $arg{name}){
    die "Need fastaObj, name,  and jobName as hash arguments to subroutine get_cog_summary : $!";
  }

  # Initialize the counting hashes
  my(%COG_Count,%Category_Count) = ();

  # Initialize the COG specific hashes
  require COGMODULE;
  my%COG_Lookup = COGMODULE->COG_Lookup();
  my%COG_Category = COGMODULE->COG_Category();

  # Populate count hashes with all possible COG's and Categories
  foreach my$cogID (sort keys %COG_Lookup){
    $COG_Count{$cogID} = 0;
  }
  foreach my$category (sort keys %Category_Count){
    $Category_Count{$category} = 0;
  }

  # Process the blast results
  my($notAssignedCount) = 0;
  
  # Count the hits from the reads onto the genomes and genuses
  foreach my$readID (@{$arg{fastaObj}->get_readOrder()}){
    my$annotation = ${$arg{fastaObj}->get_reads()}{$readID}->get_cogHit();
    if($annotation eq '????'){
      $notAssignedCount+=1;
    }
    else{
      # There could be more than one equivalent returns
      my@genes = split/\s+/,$annotation;
      my$n = scalar(@genes); # if there are n equivalently ranked genes, each one only adds 1/n to ko, ec counts

      foreach my$gene (@genes){
	# Get all COG's associated with gene
	my@COGs = split/\|C/,$gene;

	for(my$i=1;$i<=$#COGs;$i++){
	  my$cog = "C".$COGs[$i];

	  # increment the COG's count
	  $COG_Count{$cog}+= 1/$n;

	  my$category = $COG_Lookup{$cog};
	  my@categories = split(//,$category);
	  foreach(@categories){
	    $Category_Count{$_}+= (1/$n);
	  }
	}
      }
    }
  }

  # Print the summary, and generate the GenusSummary hash

  # note on floating point precision.
  # Some of the values are stored as 47.9999999999, or 50.00000001
  # This is what happens when multiple fractions are added, e.g. 1/3 three times
  # To make this more readable, I will sprintf everything to convert the .99999 and .000001 to .00000,
  #   then add 0 to turn into integer (not sure why this works, but it does)
  
  open (OUT,">$arg{jobName}\_COGcount.txt") or die "Could not open the file $arg{jobName}\_COGcount.txt: $!";
  print OUT "COG\t$arg{name}\n";
  foreach (sort {$a cmp $b} keys %COG_Count){
    my$count = sprintf("%.3f",$COG_Count{$_});
    $count+=0;
    print OUT "$_\t$count\n";
  }
  close OUT;

  open (OUT,">$arg{jobName}\_COGcategoryCount.txt") or die "Could not open the file $arg{jobName}\_COGcategoryCount.txt: $!";
  print OUT "COG Category\tDefinition\t$arg{name}\n";
  foreach (sort {$a cmp $b} keys %Category_Count){
    my$count = sprintf("%.3f",$Category_Count{$_});
    $count+=0;
    if(! defined $COG_Category{$_}{'LevelTwo'}){
      warn "No definition for $_ in cog_category\n";
    }
    else{
      print OUT "$_\t$COG_Category{$_}{'LevelTwo'}\t$count\n"; 
    }
  }
  close OUT;

  return("$arg{jobName}\_COGcount.txt","$arg{jobName}\_COGcategoryCount.txt");
}

=item $MeropsCount = get_merops_summary(fastaObj=>$FastaObject,jobName=>$jobName,name=>$name);

Every MEROPS read is associated with a six digit family ID. This routine returns a matrix of counts by familyID

=cut

sub get_merops_summary{
  my($self,%arg) = @_;
  
  unless($arg{fastaObj} and $arg{jobName} and $arg{name}){
    die "Need fastaObj, name, and jobName as hash arguments to subroutine get_merops_summary : $!";
  }

  # Initialize the hashes
  my(%Family_Count,%Family_Function) = ();

  # Get list of all possible families
  my$number_to_function = "/home/comp/jglab/bmuegge/MEROPS/MEROPS_9.5_numberToFunction.txt";
  open(IN,"<$number_to_function") or die "Could not open the file $number_to_function : $!";
  while(my$line=<IN>){
    unless($line=~/^#/){
      chomp($line);
      my($number,$function) = split/\t/,$line;
      $Family_Count{$number} = 0;
      $Family_Function{$number} = $function;
    }
  }
  close IN;

  my($notAssignedCount) = 0;
  
  # Count the hits from the reads onto the genomes and genuses
  foreach my$readID (@{$arg{fastaObj}->get_readOrder()}){
    my$annotation = ${$arg{fastaObj}->get_reads()}{$readID}->get_meropsHit();
    if($annotation eq '????'){
      $notAssignedCount+=1;
    }
    else{
      # There could be more than one equivalent returns
      my@genes = split/\s+/,$annotation;
      my$n = scalar(@genes); # if there are n equivalently ranked genes, each one only adds 1/n to ko, ec counts

      foreach my$gene (@genes){
	# In merops database, each entry has one and only one family assigned
	my($text,$family) = split/\|/,$gene;
	$Family_Count{$family}+=(1/$n);
      }
    }
  }

  # Print the summary, and generate the GenusSummary hash
  open (OUT,">$arg{jobName}\_meropsCount.txt") or die "Could not open the file $arg{jobName}\_meropsCount.txt: $!";
  print OUT "#MEROPS Family\t$arg{name}\n";
  foreach (sort keys %Family_Count){
    my$count = sprintf("%.3f",$Family_Count{$_});
    $count+=0;
    print OUT "$_\t$count\n"; 
  }
  close OUT;
  
  return "$arg{jobName}\_meropsCount.txt";
}


=item (gutGenomeCountFile,GenusCountFile) = get_gut_genome_summary(fastaObj=>$FastaObject,jobName=>$jobName,name=>$name);

This subroutine will generate a simple text file with the relative abundance of each of the gut genomes in the fasta file, normalized for genome size. The methods acts on an annotated Fasta Object, and generally should be called on a Fasta Object that has already had poor quality sequences removed.

Two text files are required, and are defined in the Constants package. GeneToGenome.txt should contain one line for every read in the GutGenome database. The first entry on the line is the fasta header for the read, and the second entry (tab delimited) is the full name of the genome to which that read belongs. The second file is GenomeSize and is a simple tab delimited file where each line has the full Genome name in column one, and the genome size in column 2.

For each genome, a count will be made of how many times a read was annotated as hitting that genome. In the event that a single read had n equivalent best blast hits in the gut database, the count will be incremented by 1/n instead of 1.

The final table reports the genome in column 1, and the count in column 2. Column 3 is the normalized relative abundance of the genome. This is calculated as (Reads assigned to Genome / Total Number of Reads)*(1/Genome Size). For reads not assigned to a genome (="unknown"), the same formula is used with Genome Size= average of all other genome sizes. The fraction is then normalized to sum to 100.

Two returns. One is file name for counts per genome. the other is file name for counts per genus (in event that multiple genomes are in same genus/species)

=cut

sub get_gut_genome_summary{
  my($self,%arg) = @_;
  
  unless($arg{fastaObj} and $arg{jobName} and $arg{name}){
    die "Need fastaObj, name, and jobName as hash arguments to subroutine get_gut_genome_summary : $!";
  }

  # Initialize the hashes
  my(%ID_to_Genome,%Name_to_Genus,%ID_Count,%GenomeSize) = ();
  
  my$totalSize = 0; # this is needed to determine average genome size
  my$notAssignedCount = 0; # this will be used for reads with no match in the gut genomes
  
  # GUTDATA is the output of script bmuegge/Microbialomics_Database_Stuff/all_genome_metadata.pl
  open(IN,Constants::GUTDATA) or die "$!";
  while(my$line=<IN>){
    chomp($line);
    unless($line=~/^#/){
      my($genomeID,$genomeName,$genusSpecies,$genomeLength) = split/\t/,$line;
      $ID_to_Genome{$genomeID} = $genomeName;
      $Name_to_Genus{$genomeName} = $genusSpecies;
      $ID_Count{$genomeID} = 0;
      $GenomeSize{$genomeID} = $genomeLength;
      $totalSize += $genomeLength;
    }
  }
  close IN;

  # Calculate the average genome size - this will be used for hits with no assignment
  my$averageSize = $totalSize / scalar(keys %GenomeSize);
  
  # Count the hits from the reads onto the genomes and genuses
  foreach my$readID (@{$arg{fastaObj}->get_readOrder()}){
    my$annotation = ${$arg{fastaObj}->get_reads()}{$readID}->get_gutGenomeHit();
    if($annotation eq '????'){
      $notAssignedCount+=1;
    }
    else{
      # Return may have more than 1 equivalent hits, separated by spaces
      my@temp = split/\s+/,$annotation;
      my$n = scalar(@temp); # If the read has more than 1 annotation, each is incremented by 1/n so
                            #  that the total contribution of the reads sums to 1
      foreach(@temp){
	$ID_Count{$_}+= (1/$n);
      }
    }
  }

  # my$totalHits = 0;
  # foreach (keys %ID_Count){
  #   $totalHits+=$ID_Count{$_};
  # }
  # $totalHits+=$notAssignedCount;

  # print "The total number of hits assigned is $totalHits. The total number of reads is ".$arg{fastaObj}->read_number()."\n\n";
  # exit;


  # Calculate the relative abundance of each genome, normalized for genome size
  # Formula for relative size = (count / total reads)*(1/genomeSize)
  # Then normalize so column sums to 100% (with no assigned counts included)
  # First, we will sum the relatizeSizes. Then go back and normalize to 100%

  # print array will have one "row" for each genome. Column 1 is name, columne 2 is count, column 3 is size, column 4 is relative abund
  my@printArray = ();
  my($genomeNumber,$totalRelativeSize) = (0,0);
  foreach my$id (sort keys %ID_Count){
    my$relativeSize = ($ID_Count{$id} / $arg{fastaObj}->read_number())*(1/$GenomeSize{$id});
    $printArray[$genomeNumber][0] = $ID_to_Genome{$id};
    $printArray[$genomeNumber][1] = $ID_Count{$id};
    $printArray[$genomeNumber][2] = $GenomeSize{$id};
    $printArray[$genomeNumber][3] = $relativeSize;

    $genomeNumber++;
    $totalRelativeSize += $relativeSize;
  }

  # Add the 'not assigned' counts to the totalRelativeSize
  my$relativeSize = ($notAssignedCount/$arg{fastaObj}->read_number()) * (1/$averageSize);
  $totalRelativeSize += $relativeSize;
  $printArray[$genomeNumber][0] = "Not Assigned";
  $printArray[$genomeNumber][1] = $notAssignedCount;
  $printArray[$genomeNumber][2] = $averageSize;
  $printArray[$genomeNumber][3] = $relativeSize;

  # Print the summary, and generate the GenusSummary hash
  my%GenusCount = ();
  open (OUT,">$arg{jobName}\_gutGenomeCount.txt") or die "Could not open the file $arg{jobName}\_gutGenomeCount.txt: $!";
  print OUT "Geneom\t$arg{name}\t$arg{name}\t$arg{name}\n";
  print OUT "Genome\tCount\tGenomeLength\tNormalized Relative Abundance\n";
  for(my$i=0; $i<=$#printArray;$i++){
    print OUT "$printArray[$i][0]\t";
    my$formattedcount = sprintf("%.4f",$printArray[$i][1]);
    $formattedcount+=0;
    print OUT "$formattedcount\t"; 
    print OUT "$printArray[$i][2]\t"; 
    print OUT sprintf("%.4f",(100*$printArray[$i][3])/$totalRelativeSize);
    print OUT "\n";

    if($printArray[$i][0] eq "Not Assigned"){
      $GenusCount{"Not Assigned"} = 100*$printArray[$i][3]/$totalRelativeSize;
    }
    elsif(defined $GenusCount{$Name_to_Genus{$printArray[$i][0]}}){
      $GenusCount{$Name_to_Genus{$printArray[$i][0]}}+=(100*$printArray[$i][3]/$totalRelativeSize);
    }
    else{
      $GenusCount{$Name_to_Genus{$printArray[$i][0]}} =(100*$printArray[$i][3]/$totalRelativeSize);
    }
  }
  close OUT;

  open (OUT,">$arg{jobName}\_GenusCount.txt") or die "Could not open the file $arg{jobName}\_GenusCount.txt: $!";
  print OUT "Genus\t$arg{name}\n";
  print OUT "Genus\tNormalized Relative Abundance\n";

  foreach(sort keys %GenusCount){
    unless($_ eq "Not Assigned"){
      print OUT "$_\t".sprintf("%.4f",$GenusCount{$_})."\n";
    }
  }
  print OUT "Not Assigned\t$GenusCount{'Not Assigned'}\n";
  close OUT;

  return ("$arg{jobName}\_gutGenomeCount.txt","$arg{jobName}\_GenusCount.txt");
}

=item (KOcountsFile.txt,ECcountsFile.txt,PathwayCountsFile.txt) = get_kegg_summary();

For all annotated reads, sumamrizes KOs, ECs, and Pathways. Return is three filenames of the newly created files.

=cut

sub get_kegg_summary{
  my($self,%arg) = @_;

  require KEGGMODULE;
  
  unless($arg{fastaObj} and $arg{jobName} and $arg{name}){
    die "get_kegg_summary needs fastaObj, name,  and jobName : $!";
  }

  # initialize the hashes
  my%KO_HOH = KEGGMODULE->KO_HOH(); # TwoD Hash. GeneID is first key. SecondKey is EC, LEvelOne,LevelTwo,PathwayNum
  my%PathwayNum_to_Name = KEGGMODULE->PathwayNum_to_Name(); # PathwayNum is first Key, Pathway Name is second
    
  # Initialize the counter hashes
  my(%KOcount,%ECcount,%PathwayCount) = ();

  # initialize the counts in each counter hash to zero
  foreach my$ko (keys %KO_HOH){
    $KOcount{$ko} = 0;
    
    # Initialize every possible EC to zero
    # Get the list of possible KO's from the second level hash keyed by 'EC'
    if ($KO_HOH{$ko}{'EC'} =~ /:/){
      my@temp = split/:/,$KO_HOH{$ko}{'EC'};
      foreach (@temp){
	$ECcount{$_} = 0;
      }
    }
    else { # one and only one EC
      $ECcount{$KO_HOH{$ko}{'EC'}} = 0;
    }
  }

  foreach my$pathwayNum (keys %PathwayNum_to_Name){
    $PathwayCount{$pathwayNum} = 0;
  }

  my($totalHits,$noAnnotation) = (0,0);
  # count the hits in the FastaObject
  foreach my$readID (@{$arg{fastaObj}->get_readOrder()}){

    my$annotation = ${$arg{fastaObj}->get_reads()}{$readID}->get_keggHit();
    if($annotation eq '????'){
      $noAnnotation +=1;
    }
    else{

      # There could be more than one equivalent returns
      my@genes = split/\s+/,$annotation;
      my$n = scalar(@genes); # if there are n equivalently ranked genes, each one only adds 1/n to ko, ec counts

      foreach my$gene (@genes){

	# Get all KO's associated with gene
	my@KOs = split/\|K/,$gene;

	for(my$i=1;$i<=$#KOs;$i++){
	  my$ko = "K".$KOs[$i];

	  # increment the KO count
	  $KOcount{$ko}+= 1/$n;

	  # Look up the EC(s) associated with this KO 
	  # I know from my construction of the KO HOH that each EC/pathway is stored once and only once per KO
	  my$ec = $KO_HOH{$ko}{'EC'};
	  my@ecs = split/:/,$ec;
	  foreach my$item (@ecs){
	    $ECcount{$item} += 1/$n;
	  }

	  # lookup and increment the pathways associated with this KO
	  my$pathwayNum = $KO_HOH{$ko}{'PathwayNum'};
	  my@pathways = split/:/,$pathwayNum;
	  foreach my$item (@pathways){
	    $PathwayCount{$item} += 1/$n;
	  }
	}
      }
    }
  }

  # Print the summary files
  # note on floating point precision.
  # Some of the values are stored as 47.9999999999, or 50.00000001
  # This is what happens when multiple fractions are added, e.g. 1/3 three times
  # To make this more readable, I will sprintf everything to convert the .99999 and .000001 to .00000,
  #   then add 0 to turn into integer (not sure why this works, but it does)
  
  open (OUT,">$arg{jobName}\_KOcount.txt") or die "Could not open the file $arg{jobName}\_KOcount.txt: $!";
  print OUT "KO\t$arg{name}\n";
  foreach (sort {$a cmp $b} keys %KOcount){
    my$count = sprintf("%.3f",$KOcount{$_});
    $count+=0;
    print OUT "$_\t$count\n";
  }
  close OUT;

  open (OUT,">$arg{jobName}\_ECcount.txt") or die "Could not open the file $arg{jobName}\_ECcount.txt: $!";
  print OUT "EC\t$arg{name}\n";
  foreach (sort {$a cmp $b} keys %ECcount){
    my$count = sprintf("%.3f",$ECcount{$_});
    $count+=0;
    print OUT "$_\t$count\n"; 
  }
  close OUT;

  open (OUT,">$arg{jobName}\_PathwayCounts.txt") or die "Could not open the file $arg{jobName}\_PathwayCounts.txt: $!";
  print OUT "PathwayNumber\tPathwayName\t$arg{name}\n";
  foreach (sort {$a cmp $b} keys %PathwayCount){
    my$count = sprintf("%.3f",$PathwayCount{$_});
    $count+=0;

    if(defined $PathwayNum_to_Name{$_}){
      print OUT "$_\t$PathwayNum_to_Name{$_}\t$PathwayCount{$_}\n";
    }
    else {
      print OUT "$_\t\t$PathwayCount{$_}\n";
      unless($_ eq "Not Defined"){
	warn "No Pathway Definition for pathwayNum $_\n";
      }
    }
  }
  close OUT;

  # undefine the big hashes made to free up memory
  undef %KO_HOH;
  
  return ("$arg{jobName}\_KOcounts.txt","$arg{jobName}\_ECcounts.txt","$arg{jobName}\_PathwayCounts.txt");
}

1;
