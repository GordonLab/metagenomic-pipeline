package Constants;

=head1 Constants

Constants: a text file to manage constant parameters used by the Metagenomic Pipeline.

=head1 Version

Version 2.0; Last Modified July 22, 2011.

=cut

sub last_modified {
  return "July 22, 2011";
}

sub version {
  return "2.0";
}

# begin code

use strict;
use warnings;

# Blast Databases
use constant BLASTPROG => 'blastall';
#use constant COGDB => '/srv/cgs/data/STRING/9.0_modHeader/COG';
use constant COGDB => '/srv/cgs/data/STRING/9.0_onlyAnnotated/COG';
use constant COGDBLENGTH => '1903618658';
#use constant KEGGDB => '/srv/cgs/data/KEGG/version58/KEGG58';
#use constant KEGGDB => '/srv/cgs/data/KEGG/v58_smallVolumes/KEGG58';
#use constant KEGGDB => '/srv/cgs/data/KEGG/v58_ModHeader/KEGG58';
use constant KEGGDB => '/srv/cgs/data/KEGG/v58_onlyAnnotated/KEGG';
use constant KEGGDBLENGTH => '2214788408';
use constant MEROPSDB => '/srv/cgs/data/MEROPS/version9.5/MEROPS';
#use constant NRDB => '/srv/cgs/data/nr/nr';
use constant GUTDB => '/srv/cgs/data/GutGenomes/128_Genomes/128_Genomes';
use constant GUTDATA => '/srv/cgs/data/GutGenomes/128_Genomes/128_Genomes_Metadata.txt';
#use constant HOSTDB => '/srv/cgs/data/Hsa/Hsa_ref.fasta';
use constant HOSTDB => '/srv/cgs/data/human_genomic_05272010/human_genomic';

#use constant HUMAN => '/srv/cgs/data/Hsa/Hsa_ref.fasta';
#use constant MOUSE => '/srv/cgs/data/bmuegge/Mouse/Mouse';

# Gut Genome resources
use constant GeneToGenome => '/srv/cgs/data/GutGenomes/128genomes_May2011/GeneToGenome.txt';
use constant GenomeSize => '/srv/cgs/data/GutGenomes/128genomes_May2011/GenomeSize.txt';

# Blast parameters
use constant BLAST_E => 1e-5; # minimum e-score
use constant BLAST_B => 10; # number of returns
use constant BLAST_V => 10; # number of one line returns
use constant BLAST_M => 8; # tab delimited output
use constant BLAST_MINID => 50; # minimum percent identity required between query and subject
use constant BLAST_MINALIGN => 0; # minimum alignment length required between query and subject, as percentage of total query legnth
use constant BLAST_MINSCORE => 50; # bit score from blast

# Dereplicator program and paramaters
use constant EXTRACT_REPLICATES => '/home/comp/jglab/bmuegge/bin/replicates/scripts/extract_replicates.py';
use constant REPLICATE_LENGTH => 0;
use constant REPLICATE_PERCENT => 0.97;   
use constant REPLICATE_START => 20; # these three numbers are Tanya's standard values

# Fasta de-multiplexing and quality filtering
use constant MINIMUM_SEQUENCE_LENGTH => 60;
use constant RL_MID_CONFIG_PARSE => '/home/comp/jglab/bmuegge/Metagenomic_Pipeline/MIDConfig_bdmMod_RL0.parse'; 

# Listener parameters
use constant WAIT_TIME => 120; # 300 seconds = 5 minutes
use constant NO_CHANGE_IN_FILE_FOR_X_SECONDS => 60;

# split_run_check and parameters
use constant SPLIT_RUN_CHECK => '/home/comp/jglab/bmuegge/annotation/split_run_check_combine.pl';
use constant NUM_SEQ => 150;  # I'm setting a low number so searches against large KEGG and COG databases can definitely get done in < 1 hour despite 2Gb RAM load

1;

=head1 AUTHOR

Brian Muegge

=head1 COPYRIGHT

Copyright (c) 2011, Washington University School of Medicine

=cut
