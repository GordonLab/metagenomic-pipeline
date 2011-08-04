#!/usr/bin/perl -w


#    Script: split_run_check_combine.pl
#    ___________________________________________________________________________
#
#    Version 0.1
#
#    Copyright (C) 2008-2009 Brian Muegge and Jeremiah Faith
#
#    http://hamlet.wustl.edu/microbialomics_dev/docs_db/
#
#    About: License
#
#       Licensed under the GNU General Public License
#
#        This program is free software; you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation; either version 2 of the License, or
#        (at your option) any later version.
#
#        This program is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with this program; if not, visit http://www.gnu.org/licenses/gpl.txt
#        or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
#        Boston, MA  02111-1307  USA.
#
# Topic: Getting started
# 
#
# (start code) 
#  perl split_run_check_combine.pl -i <input (in fasta or scarf format)> -p <program and parameters (put in single quotes; include INCLUDE_INFILE)> -n <num_lines_per_job>
#
#  # Example (will blastp proteins.faa against the string database using 40 sequences per batch job)
#  perl split_run_check_combine.pl -i proteins.faa -p '/srv/cgs/local/i386/ncbi-blast/latest/bin/blastall -p blastp -d /srv/cgs/data/STRING_v8.0/protein.sequences.v8.0.fa -i INCLUDE_INFILE -m 8 -e 10e-10' -n 40
# (end)
#
# This script aims to provide a robust/generalizable framework for starting and checking jobs on the 
# cluster. It splits a FASTA or scarf file into many smaller files, runs a user-defined program on each small file, checks to make sure that each of the programs
# run on the smaller files completed successfully, and if they did complete successfully, it combines all of the outputs into a single file.
# 
# You specify where to place the infile name for each subfile in your program using the INCLUDE_INFILE.  You should try to set -n so that jobs on the cluster take around 15 minutes.
#
# Topic: How it works
# This script splits a scarf file or a fasta file into many smaller files.  It run a particular program 
# creating a successful-completion file *.complete when finished and a *.start file before running (in case the job doesn't finish)
# run a job that waits for the first job to finish and checks to make sure everything is ok (otherwise it should re-run the failed jobs, but this isn't implemented yet)
# if all of the jobs are ok, it should combine all of the outputs and remove the temp files from the split
#
# IMPORTANT: in your program program parameters put the word INCLUDE_INFILE whereever you want the infile to be placed
# this script assumes your program prints the output to STDOUT
#

use strict;
#use Getopt::Long;
use Getopt::Long qw(:config no_ignore_case);
use FindBin qw($Bin);

my $infile;
my $num_lines_per_job;
my $program;
my $mode = 'split';
my $code = "c" . int(rand() * 100000); # used to track things
my $file_suffix;
my $infile_suffix;
my $this_program = "$Bin/split_run_check_combine.pl";
my $outfile;
my $memory;
my $use_long_queue;
my $node_arch;
my $JOB_ID;

GetOptions ( #"s" => \$start_job,
            "infile=s"   => \$infile,
            "num_lines_per_job=i"   => \$num_lines_per_job,      # integer
            "code=s"   => \$code,      # integer
	    "params=s"     => \$program,
	    "suffix=s" => \$file_suffix, # used after combining the files to 
	    "Suffix_infile=s" => \$infile_suffix, 
	    "mode=s" => \$mode,
	    "Memory=i" => \$memory,
	    "LongQueue" => \$use_long_queue,
	    "arch_node=s" => \$node_arch,						#Mod to define the node arch
	    "outfile=s" => \$outfile,
	    "jobid=i" => \$JOB_ID,
            );  # flag


#print "have memory $memory long queue $use_long_queue\n";

#$code = $code .'_'. $file_suffix if $file_suffix;



if ($mode eq 'split') { # in split mode, split the file, start the job and start a job to check the jobs 
	usage() unless $infile && $program && $num_lines_per_job;
#	die "usage: perl $this_program -i <input (in fasta or scarf format)> " .
#	"-p <program and parameters (put in single quotes; include INCLUDE_INFILE)> -n <num_lines_per_job> -s [file suffix for joined results]\n" unless $infile && $program && $num_lines_per_job;

	print STDERR "code is $code\n";
	mkdir "SGE";
	
	# split the big sequence file into a bunch of little ones
	my ($infile_suffix, @split_files) = split_jobs($infile, $num_lines_per_job, $code);
	my $job_id = run_jobs(\@split_files, $program, $code, $memory, $use_long_queue, $node_arch);

	my $num_jobs = scalar @split_files;
	append_job_check($code, $num_jobs, $job_id, $infile, $file_suffix, $infile_suffix);
}
elsif ($mode eq 'recheck') { # in split mode, split the file, start the job and start a job to check the jobs 
	# blah
	my ($num_jobs, $infile, $file_suffix, $infile_suffix, $lines_jobs, $jobs_to_rerun) = get_recheck_jobs($code, 1);

	# remove all of the complete and start jobs so they can be properly recreated in "run" mode
	system("rm $code.*.complete $code.*.start");

#	print "have num jobs $num_jobs\nhave rerun jobs:\n$jobs_to_rerun\n";

	# rerun the failed jobs and check for completion
	my $job_id = run_jobs("", "", $code, $memory, $use_long_queue, $node_arch, $jobs_to_rerun);
	append_job_check($code, $num_jobs, $job_id, $infile, $file_suffix, $infile_suffix);
}
elsif ($mode eq 'run') { # run mode is called by the jobs in split mode; it simply runs each job of the user's program; creating a *.start file and a *.complete file which are used by the check mode 
	my $hostname = my $host = `hostname`;
#	print "in run on host $hostname\n";

	# create a start file with the stuff to be executed in case something fails and the job needs to be restarted

	my $sub_id = parse_sub_id($program, $code);

	# recommended by Brian (thanks to the heads up from Alejandro)
	# this allows you to use just the program name (and NOT the full path) on the different nodes
	# and it imports the correct environment for the node (i.e. for 32 bit it looks in the 32-bit path)
	# The line is commented since apparently is no longer needed on the CGS cluster. Each node should be able to identify the adequate path
	# eval('cgsenv -f -b');

	system("echo '$program' > $code.$sub_id.start");
	# run the program with the parameters passed by the split mode
	my $status = system("$program");

	# check return status
	if ( $? == -1 ) {  # program never executed
		my $error_msg  = sprintf "Error: command ($program) failed to execute: $!\n";
		system("echo '$error_msg' > $code.$sub_id.failed");
	}
	elsif ($status != 0) { # program executed but failed for some reason
		my $exit_value = $? >> 8;
		my $error_msg = sprintf "command failed (did not exit with status zero); exited with value %d", $exit_value;
		system("echo '$error_msg' > $code.$sub_id.failed");
	}
	else { # if the program completed successfully write a complete file
		my $exit_value = $? >> 8;
		printf "command exited with value %d", $exit_value;
		system("touch $code.$sub_id.complete");
	}
}
elsif ($mode eq 'check') { # in check mode, we just make sure everyfile has a *.start and *.complete and concatenate the outputs
	my $num_files = $num_lines_per_job;
	my ($failed, $files_to_delete) = check_files($code, $num_files, $JOB_ID);
	print "to delete @$files_to_delete\n";
	if (scalar @$failed) {
		create_recheck_files($failed, $code, $files_to_delete, $num_files, $infile, $file_suffix, $infile_suffix);
		unlink(@$files_to_delete);
		die "Error missing some files: @$failed\n";
	}
	else {
		combine_files($code, $num_files, $infile, $outfile, $file_suffix, $files_to_delete, $infile_suffix)
	}
}

sub get_recheck_jobs { 
	my $code = shift;
	my $get_rerun_jobs = shift;

	my $recheck_file = "$code.recheck";

	open (IN, $recheck_file) or die "can't open $recheck_file: $!\n";
	my $num_lines = <IN>;
	chomp $num_lines;
	my $infile = <IN>;
	chomp $infile;
	my $file_suffix = <IN>;
	chomp $file_suffix;
	my $infile_suffix = <IN>;
	chomp $infile_suffix;

	my $text = "";	
	my @lines;
	while (<IN>) {
		chomp;
		s/\s+//g;
		$text .= $_;
		@lines = split /,/, $text;
	}
	close IN;

	my @lines_jobs;
	for my $l (@lines) {
		my ($line, $subjob) = split /:/, $l;
		push @lines_jobs, [$line,$subjob];
	}

	# make things easy to check for using a hash
	my %need_to_recheck;
	for my $l (@lines_jobs) {
		$need_to_recheck{$l->[0]} = 1;
	}


	my $recheck_jobs = "";
	my $old_submission = "$code.sh";


	if ($get_rerun_jobs) {
		# go into the old submission file and grab out the jobs that failed
		open (IN, $old_submission) or die "can't open $old_submission: $!\n";
		my $count = 0;
		while (<IN>) {
			if ($need_to_recheck{$count}) {
				$recheck_jobs .= $_;
			}
			$count++;
		}
		close IN;
	}
#	print "lines is @lines from $text from file $recheck_file jobs are:\n$recheck_jobs\n";

	return $num_lines, $infile, $file_suffix, $infile_suffix, \@lines_jobs, $recheck_jobs;
}

sub create_recheck_files { 
	my ($failed, $code, $files_to_delete, $num_files, $infile, $file_suffix, $infile_suffix) = @_;


	my $recheck_file = "$code.recheck";
	open (OUT, ">$recheck_file") or die "can't open $recheck_file: $!\n";
	my $to_check = join ",", @$failed;
	
	# first line is the number of files in the original query
	print OUT "$num_files\n";
	# second line is the name of the original infile
	print OUT "$infile\n";
	# third line is the file_suffix to be used by the outfile
	print OUT "$file_suffix\n";
	# fourth line is the file_suffix to be used by the outfile
	print OUT "$infile_suffix\n";
	# second line is id of each subfile for the resubmission
	print OUT "$to_check\n";
}

# concatenate all of the results in order
#sub combine_files { 
#	my ($code, $num_files, $infile, $outfile, $files_to_delete) = @_;
#
#	my $outfiles = "";
#	for my $i (0 .. $num_files-1) {
#		push @$files_to_delete, "$infile.$code.$i";
##		system("rm $infile.$code.$i"); # remove the sequence subset file
#		$outfiles .= " $infile.$code.$i.out.gz";
#		my $gz_outfile = "$infile.$code.$i.out.gz";
#		system("zcat $gz_outfile >> $outfile"); # concatenate all of the files into outfile
#		push @$files_to_delete, $gz_outfile;
##		system("rm $gz_outfile"); # remove the sequence subset file
#	}
#
	#my $recheck_file = "$code.recheck";
	# second line is id of each subfile for the resubmission
	#print OUT "$to_check\n";
#}#

# concatenate all of the results in order
sub combine_files { 
	my ($code, $num_files, $infile, $outfile, $file_suffix, $files_to_delete, $infile_suffix) = @_;
	#	combine_files($code, $num_files, $infile, $outfile, $file_suffix, $files_to_delete)

	my $outfiles = "";
	for my $i (0 .. $num_files-1) {
#		system("rm $infile.$code.$i"); # remove the sequence subset file
		my $gz_outfile;
		if ($infile_suffix) {
			push @$files_to_delete, "$infile.$code.$i.$infile_suffix";
			$outfiles .= " $infile.$code.$i.$infile_suffix.out.gz";
			$gz_outfile = "$infile.$code.$i.$infile_suffix.out.gz";
		}
		else {
			push @$files_to_delete, "$infile.$code.$i";
			$outfiles .= " $infile.$code.$i.out.gz";
			$gz_outfile = "$infile.$code.$i.out.gz";
		}

		if ($i == 0) {
			system("zcat $gz_outfile > $outfile"); # concatenate all of the files into outfile (first file creates a new file)
		}
		else {
			system("zcat $gz_outfile >> $outfile"); # concatenate all of the files into outfile
		}
		push @$files_to_delete, $gz_outfile;
#		system("rm $gz_outfile"); # remove the sequence subset file
	}

	my $recheck_file = "$code.recheck";
	if (-e $recheck_file) { # don't need the recheck file anymore
		push @$files_to_delete, $recheck_file;
	}

	push @$files_to_delete, "$code.check.sh";
	push @$files_to_delete, "$code.sh";
	push @$files_to_delete, "$code.check_status";

	print "deleting @$files_to_delete\n";

	# now that everything is joined, delete all of the files
	unlink(@$files_to_delete);
	#foreach my $f (@$files_to_delete){ unlink($f) or die "Can't delete $f : $!";} 
}

# make sure that all of the files are there otherwise return the missing ones so they can be rerun if desired
sub check_files {
	my ($code, $num, $JOB_ID) = @_;

	my $out = "$code.check_status";
	open (OUT, ">$out") or die "can't open $out: $!\n";
	my @failed_jobs;
	my $num_failed = 0;
	my @files_to_delete;
	my $sge_err;
	my $infile;


	# CHECK FROM HERE, NEED TO MAKE SURE THIS WORKS

	my $num_to_check = $num-1;

	# not going through the normal way, need to only scan the recheck files
	my ($num_jobs, $lines_jobs);
	if (-e "$code.recheck") {
	        ($num_jobs, $infile, $file_suffix, $infile_suffix, $lines_jobs) = get_recheck_jobs($code);

		print "using recheck\n";
		$num_to_check = $#$lines_jobs;	
	}

#	for my $i (0 .. $num-1) {
	for my $i (0 .. $num_to_check) {
		my $line = $i;
		my $subjob = $i;

		# if working from a recheck, need to use a different coordinant system
		if ($lines_jobs) {
#			$line  = $lines_jobs->[$i][0];
			$subjob = $lines_jobs->[$i][1];
		}

		my ($start, $complete, $sge_err_OK) = (0,0,0);
		my $start_file = "$code.$subjob.start";
		if (-e $start_file) { 
			#print OUT "found\n";  # printing success cluttered everything up
			$start=1;
		}
		else { 
			print OUT "$start_file not found\n"; 
		}

		my $complete_file = "$code.$subjob.complete";
		if (-e $complete_file) { 
	#		print OUT "found\n"; # printing success cluttered everything up
			$complete=1;
		}
		else { 
			print OUT "$complete_file not found\n"; 
		}

		# CHECK FOR COMPLETION OF TIME STAMPS IN THE SGE ERROR FILE
		
		my $job_element = $line+1;   # the SGE number scheme is 1 based NOT 0 based
		my $sge_err_file = "./SGE/$code.sh.e$JOB_ID.$job_element";
		my $sge_out_file = "./SGE/$code.sh.o$JOB_ID.$job_element";

		#my $sge_err_file = "~/.sge/$code.sh.e$JOB_ID.$job_element";
		#my $sge_out_file = "~/.sge/$code.sh.o$JOB_ID.$job_element";
		($sge_err_OK, $sge_err) = check_sge_error_file($sge_err_file);
		if (!$sge_err_OK) {
			print OUT "sge error for $sge_err_file (job $code.$subjob)\n$sge_err\n";
		}

		if ($start && $complete && $sge_err_OK) {  # remove the check files if we have all three
			push @files_to_delete, $start_file;
			push @files_to_delete, $complete_file;
			push @files_to_delete, $sge_err_file;
			push @files_to_delete, $sge_out_file;
			print "adding $start_file $complete\n";
		}
		else {
			push @failed_jobs, "$line:$subjob";  # save the jobs that failed
		}
	}

	if (scalar @failed_jobs) {
		print OUT "\n\nSummary: failed jobs @failed_jobs\n";
	}
	else {
		print OUT "\n\nSummary: no failed jobs\n";
	}
	close OUT;

	return \@failed_jobs, \@files_to_delete;
}

# not perfect, but should work
sub check_sge_error_file {
	my $in = shift;
	my $OK = 1;
	my $NOT_OK = 0;

	my ($real, $user, $sys, $finished);
	my $text;

	# check if the file exists
	# NOTE I do all of these file readings and checkings using linux commands
	# because I'm using the ~ (tilde) and I want the BASH shell to sort out the home directories
	unless (`if test -e $in; then printf "1"; fi`) {
#		print "HERE no file\n";
		return $NOT_OK, "SGE error file $in doesn't exist";	
	}

	my @lines = `cat $in`;
	for (@lines) {
#		print "reading:\t$_";
		$text .= $_;
		if (/^real\b/) { $real=1; }
		elsif (/^user\b/) { $user=1; }
		elsif (/^sys\b/) { $sys=1; }
		elsif (/^Finished:/) { $finished=1; }
	}
	close IN;


#	print "have R $real U $user S $sys F $finished for $in\n";
	if ($real && $user && $sys && $finished) {
		return $OK, "";
	}
	else {
		return $NOT_OK, $text;
	}
}

# start a job to check that all of the previous jobs completed; this job starts after the first job-array completes
sub append_job_check {
	my ($code, $num_jobs, $job_id, $infile, $file_suffix, $infile_suffix) = @_;

	my $outfile = $infile . ".$code.combined";
	if ($file_suffix) {
		$outfile = $infile . ".$file_suffix";
	}
	
	# have the output checked
	my $program = "$this_program -c $code -n $num_jobs -m check -i $infile -o $outfile -j $job_id";
	$program .= " -s $file_suffix" if $file_suffix;
	$program .= " -S $infile_suffix" if $infile_suffix;

	system("echo '$program' > $code.check.sh");
	#my $qsub = 'qsub -e $HOME/.sge -o $HOME/.sge -V -cwd -r y -hold_jid' .  " $job_id";
#	my $qsub = 'qsub -e $PWD/SGE -o $PWD/SGE -cwd -l arch=lx24-amd64 -r y -hold_jid' .  " $job_id";
	my $qsub = 'qsub -e $PWD/SGE -o $PWD/SGE -cwd -r y -hold_jid' .  " $job_id";
	system("$qsub $code.check.sh");
}

sub parse_sub_id {
	my $text = shift;
	my $code = shift;

	if ($text =~ /\.$code\.(\d+)/) {
		return $1;
	}

	return;
}


#$outfile = "$scarf_file.$suffix";
#print SH_OUT "$program -q $outfile -D $db $params | gzip > $outfile.out.gz\n";

#close SH_OUT;

# start the job on the cluster
#if ($start_job) {
#	system("nq $sh_out | qsub");
#}

# set up a batch file to run all of the jobs with nq;  rather than run the user's job directly, 
# we're going to all this script to run the job, so we can create a *.start and *.complete file
# to report on the progress/failure of each job
sub run_jobs {
	my ($files, $program, $code, $memory, $use_long_queue, $node_arch, $recheck_jobs) = @_;

	# make a batch file
	my $batch_file = $code . ".sh";
	print "opening $batch_file: $!\n";
	open (OUT, ">$batch_file") or die "can't open: $batch_file: $!\n";

	if ($recheck_jobs) {
		print OUT "$recheck_jobs";
	}
	else {
		for my $f (@$files) {
			my $command = $program;	
			$command =~ s/INCLUDE_INFILE/$f/g;
			if ($command =~ /INCLUDE_OUTFILE/){							#Modified to add the possibility of an output file on the command line of the program
				$command =~ s/INCLUDE_OUTFILE/$f.out/g;
				print OUT "$this_program -c $code -m run -p '$command ; gzip $f.out'\n";
			}else{	
				print OUT "$this_program -c $code -m run -p '$command | gzip > $f.out.gz'\n";
			}
		}
	}

	close OUT;

	my $params = '-e $PWD/SGE -o $PWD/SGE ';
	$params .= "-l h_vmem=$memory" . "G" if $memory;
	$params .= " -P long" if $use_long_queue;

#	$node_arch = "arch=lx24-amd64";
	#$node_arch = "";
	$params .= " -l $node_arch" if $node_arch;
	

#	print "using params $params\n";
	system("nq $batch_file | qsub $params");

	my $id;
	while (!($id = get_job_id($code))) {
		sleep(10);
	}

	return $id;
}

sub get_job_id {
	my $code = shift;

	print STDERR "getting job id\n";
#	sleep(10);

	my @jobs = `qstat`;

#	my $job_id;
#	my @stuff;
	my $jobid;
	for my $j (@jobs) {
		$j =~ s/^\s+//;
		$j =~ s/\s+$//;
		if ($j =~ /$code/) {
			my ($jid, @stuff) = split /\s+/, $j;	
			$jobid = $jid;
#			print "have jobid $jobid from @jobs\n";
			print STDERR "found job id $jobid\n";
			return $jobid;
		}
	}
	print STDERR "found job id $jobid\n";

	return;	
}

# split the file into smaller files with $num_lines_per_job in each file
sub split_jobs {
	my ($infile, $num_lines_per_job, $code) = @_;

	my $file_type = get_file_type($infile);
	my $numlines = get_num_seqs($file_type, $infile);
	my $num_jobs = int($numlines / $num_lines_per_job);
	printf "splitting $numlines sequence $infile into %d jobs\n", $num_jobs+1;
	my @files; # of format infile: outfile
	my $file_suffix;


	if ($file_type eq 'scarf') {
		@files = split_scarf($infile, $num_lines_per_job, $code);
	}
	elsif ($file_type eq 'fasta') {
		($file_suffix, @files) = split_fasta($infile, $num_lines_per_job, $code);
	}

	return $file_suffix, @files;
}


# split a scarf file into smaller files with $num_lines_per_job in each file
sub split_scarf {
	my ($scarf_file, $num_lines_per_job, $code) = @_;
	my @files;
	my $suffix = 0;
	my $line_num=1;

	open (IN, $scarf_file) or die "can't open $scarf_file: $!\n";

	my $prefix = "$scarf_file.$code";
	my $outfile = "$prefix.$suffix";
	# open the first file
	open (OUT, ">$outfile") or die "can't open $scarf_file: $!\n";
	push @files, $outfile;
	while(<IN>) {
		if (0 == ($line_num++ % $num_lines_per_job)) {  # reached the desired number for a file
			close OUT;	
			$suffix++;
			$outfile = "$prefix.$suffix";
			open (OUT, ">$outfile") or die "can't open $outfile: $!\n"; # open a new file
			push @files, $outfile; # store the file name
		}
		print OUT $_;
	}
	close OUT;

	return @files;
}

sub get_file_suffix {
	my $file = shift;

	my @pieces = split /\./, $file;
	if (scalar @pieces > 1) {
		return $pieces[$#pieces];
	}

	return;
}

# split a fasta file into smaller files with $num_lines_per_job in each file
sub split_fasta {
	my ($fasta_file, $num_lines_per_job, $code) = @_;
	my @files;
	my $suffix = 0;
	my $line_num=1;
	my $file_suffix = get_file_suffix($fasta_file);

	open (IN, $fasta_file) or die "can't open $fasta_file: $!\n";

	my $prefix = "$fasta_file.$code";
	my $outfile = "$prefix.$suffix";
	$outfile .= ".$file_suffix" if ($file_suffix); # retain the user's file suffix
	open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
	push @files, $outfile; # store the file name
	while(<IN>) {
		if (/^>/ && 0 == ($line_num++ % $num_lines_per_job)) {
			close OUT;
			$suffix++;
			
			$outfile = "$prefix.$suffix";					
			$outfile .= ".$file_suffix" if ($file_suffix); # retain the user's file suffix
			open (OUT, ">$outfile") or die "can't open $outfile: $!\n";
			push @files, $outfile; # store the file name
		}
        	print OUT $_;
        }
	close OUT;

	return $file_suffix, @files;
}



# determine if file is fasta or scarf (i.e. solexa)
sub get_file_type {
	my $in=shift;	
	open (IN, $in) or die "can't open $in: $!\n";

	while (<IN>) {
		next if /^\s*$/; # skip blank lines

		if (/^>/) {
			close IN;
			return 'fasta';
		}
		elsif (/:\d+:\d+:\d+:/) {
			close IN;
			return 'scarf';
		}
		else {
			close IN;
			die "Could not recognize file format for: $in\n";
		}
	}

}


sub get_num_seqs {
	my ($type, $infile) = @_;

	if ($type eq 'fasta') {
		return get_num_seqs_fasta($infile);
	}
	elsif ($type eq 'scarf') {
		return get_num_seqs_scarf($infile);
	}
}

sub get_num_seqs_scarf {
	my $in=shift;

	open (IN, $in) or die "can't open $in: $!\n";
	
	my $line_count=0;
	while (<IN>) { $line_count++; }
	close IN;

	return $line_count;
}

sub get_num_seqs_fasta {
	my $in=shift;

	open (IN, $in) or die "can't open $in: $!\n";
	
	my $line_count=0;
	while (<IN>) { $line_count++ if /^>/; }
	close IN;

	return $line_count;
}


sub usage {
	print "\n\tUsage: perl split_run_check_combine.pl\n".
	"\t\t-i <input (in fasta or scarf format)>\n" .
	"\t\t-p <program and parameters (put in single quotes; include INCLUDE_INFILE)>\n".
	"\t\t-n <num_lines_per_job>\n".
	"\t\t-s [file suffix for joined results]\n".
	"\t\t-M [minimum memory requirements in GB]\n".
	"\t\t-L [use long queue]\n".
	"\t\t-a [architecture of the node needed. i.e.: -a arch=lx24-amd64]\n".
	"\t\t-m [recheck (to recheck, only provide -c <code> and optionally -M and -L)]\n\n".
	"\t\t if the program requires an output file just write INCLUDE_OUTFILE where the outfile should go. The joined output will be still Input.suffix\n\n".
	"\tExamples:\n\n".
	"\t# split test.faa info 5 sequences per file and blast against proteinDB, joining the results into a file called test.faa.blastp\n".
	"\tperl split_run_check_combine.pl -i test.faa -p 'blastall -p blastp -d proteinDB -i INCLUDE_INFILE -m 8 -e 10e-10' -n 5 -s blastp\n\n".
	"\t# try to recover failed jobs from c4789, using the long queue and nodes with at least 6GB of RAM\n".
	"\tperl split_run_check_combine.pl -c c4789 -m recheck -M 6 -L\n".
	"\n";
	exit(1);
}
