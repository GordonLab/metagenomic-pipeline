#!/usr/bin/python

import batch_replicates_config
import fasta
import os
import subprocess
import random
import tempfile
import sys
import time

#from optparse import OptionParser

'''
This script takes a 454 fasta file as input, determines which reads
are artifactual replicates and outputs a filtered set that only
contains legitimate reads and summary files about the sequences that
were filtered out.

Artifactual reads are those that are almost identical and start at the
same position.

Clusters are first identifed using CD-HIT, then clusters are accessed
to determine if they are truly # artifactual by determining if the
initial base pairs of clustered sequences match.

This script requires:
  an installation of cd-hit
  batch_replicates_config.py
  fasta.py
  cd_hit_parse.py
  extract-clusters-html.py


Changing the input to use option flags rather than sys.argv position.  This is not yet implemented.

# This uses the -i flag to indicate the input file name, since it might contain spaces
#parser = OptionParser()
#parser.add_option("-i", "--input", dest="filename", help="FASTA input file name")
#parser.add_option("-o", "--output", dest="dirname", help="Output directory")
#parser.add_option("-c", "--cutoff", dest="cutoff_input", help="The <sequence identity cutoff> is cd-hits global sequence identity calculated as the number of identical amino acids in the alignment divided by the full length of the shorter sequence.  This value should be between 0.85 and 1.0.", action="store", type="float", default="0.9")
#parser.add_option("-l", "--length", dest="length_input", help="The <length difference requirement> is the percent of the sequence that is required to match. 0 is the default and means that there is no length restrictions.  This allows for reads of variable lengths.  1 will require that all the sequences in a cluster are the same length.  This value should be between 0 and 1.0.")

#parser.add_option("-b", "--basepair", dest="bp_input", help="The <initial base pair requirement> is the number of base pairs required to match at the beginning of each sequence.")


#(options, args) = parser.parse_args()
##########
'''

if (len(sys.argv) != 6):

   print """

-------------------

Usage: extract_replicates.py <input filename> <sequence identity cutoff> <length difference requirement> <initial base pair requirement> <output directory>

The input file name should not contain spaces.
   
The <sequence identity cutoff> is cd-hits global sequence identity calculated as the number of identical amino acids in the alignment divided by the full length of the shorter sequence.  A good value to start with is 0.9.
   
The <length difference requirement> is the percent of the sequence that is required to match. 0 is the default and means that there is no length restrictions.  This allows for reads of variable lengths.  1 will require that all the sequences in a cluster are the same length.

The <initial base pair requirement> is the number of base pairs required to match at the beginning of each sequence.  A good value to start with is 3.

-------------------

"""
   sys.exit(1)


filename = sys.argv[1]
cutoff_input = sys.argv[2]
length_input = sys.argv[3]
bp_input = sys.argv[4]
dirname = sys.argv[5]

filename_check = filename.split(" ")
if len(filename_check) > 1:
   print '\nThe input filename should not contain spaces.  Please rename it, and try again.\n'
   sys.exit(2)



#
# Function to make a directory for the output files if it doesn't exit already.
#

def _mkdir(root, newdir):
   try:
      os.mkdir(root)
   except OSError:
      print "\nCannot make the directory", root, ".  Does it already exist?\n"
      sys.exit(2)
   os.mkdir(newdir)

#def _mkdir(root, newdir):
#   if os.path.isdir(newdir):
#      print "\nThe directory '%s' cannot be created.  Does it already exist?" % (root)
#      sys.exit(2)
#   elif os.path.isfile(newdir):
#      print "A file with the name %s already exists.  Please use a different output directory." % (newdir)
#      sys.exit(2)
#   else:
#      head, tail = os.path.split(newdir)
#      if head and not os.path.isdir(head):
#         _mkdir(root, head)
#      if tail:
#         os.mkdir(newdir)


#
# Make the output directory 
#

_mkdir(dirname, dirname+'/tmp')


# Check to make sure the sequence identity cutoff and length difference requirement values
# are within acceptable ranges

try:
    cutoff_test = float(cutoff_input)
    if (cutoff_test > 1.0 or cutoff_test < 0.85):
        print "Please input a cutoff value between 0.85 and 1.0"
        sys.exit(2)
    else:
        cutoff = cutoff_test

except ValueError:
    print "Please input a cutoff value between 0.85 and 1.0"
    sys.exit(2)
 

try:
    length_test = float(length_input)
    if (length_test > 1.0):
        print "Please input a length requirement value between 0 and 1.0"
        sys.exit(2)
    else:
        length = length_test
except ValueError:
    print "Please input a length requirement value between 0 and 1.0"
    sys.exit(2)


# Check to make sure the input file exists and can be opened

try:
   fasta_file = open(filename)
except:
   print "This file could not be opened"
   sys.exit(2)


# Check to make sure the input file is in FASTA format

try:
   fasta_data = fasta.load(fasta_file)
except:
   print 'This file does not seem to be a fasta file.  Please try again with a fasta file'
   sys.exit(0)


# Write out the FASTA input file to a new file, because CD-HIT doesn't handle
# all input file types correctly


new_fasta_file = dirname+'/tmp/input_fasta_file.fa'

try: 
   input_fasta_file = open(new_fasta_file, 'wt')
except:
   print 'Cannot open', new_fasta_file, 'for writing a temporary fasta file'
   sys.exit(2)

for i in fasta_data:
   input_fasta_file.write('>%s\n%s\n' % (i, fasta_data[i]))

input_fasta_file.close()


# The output file from cd-hit is written to a tmp directory for this session

cdhit_command = '%s/cd-hit-est -i %s/tmp/input_fasta_file.fa -o %s/tmp/cdhit_output_temp -c %s -n 8 -s %s -d 0 -M 1000' % (batch_replicates_config.cdhit_dir, dirname, dirname, cutoff, length)


# Run CD-HIT
prog = subprocess.Popen(cdhit_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

(stdout, stderr) = prog.communicate()

if (prog.returncode != 0):
   print 'cd-hit failed'
   print stderr
   sys.exit(2)


# Output cd-hit output and errors to a tmp directory

fp = open(dirname+'/tmp/cd-hit.out', 'w')
fp.write(stdout)
fp.close()

fp = open(dirname+'/tmp/cd-hit.err', 'w')
fp.write(stderr)
fp.close()



# Evaluate CD-HIT files
prog2 = subprocess.Popen('/home/comp/jglab/bmuegge/bin/replicates/scripts/extract-clusters-html.py %s/tmp/cdhit_output_temp.clstr  %s %s/extracted_clusters %s text -i "%s"' % (dirname, filename, dirname, bp_input, filename), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

(stdout, stderr) = prog2.communicate()

if (prog2.returncode != 0):
   print 'extract-clusters-html.py failed'
   print stderr
   sys.exit(2)

#print './extract-clusters-html.py %s/tmp/cdhit_output_temp.clstr  %s %s/extracted_clus\
#ters %s text -i "%s"' % (dirname, filename, dirname, bp_input, filename)


print stdout

fp = open(dirname+'/tmp/extract.err', 'w')
fp.write(stderr)
fp.close()

if stderr:
   print "There was a problem with the analysis.  Check %s/tmp/extract.err for details\n." % (dirname)

else:
   print "Your results are in the directory: %s\n" % (dirname)
   print "See the README file for more information on the output files.\n"




