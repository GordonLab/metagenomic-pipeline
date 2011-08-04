#! /usr/bin/env python

import string
import sys
import re
import fasta
from operator import itemgetter
# from heapq import nlargest

# a library for parsing the cd-hit output
import cdhit_parse

from optparse import OptionParser


# This script takes CD-HIT output and creates four output files
# *.fasta_clusters is a file with all the clusters in fasta format, sorted from clusters with the 
# most sequences to those with the least
# *_unique.fa is a fasta file of all the unique sequences, taking the representative sequence 
# from each cluster
# *.cluster_summary is a summary of the sequences and the number of clusters in each file
# *.cluser_sizes is a list of the number of clusters of each size

"""
Usage: extract-clusters-html.py <filename.clstr> <filename.fa> <output_file>
<initial base pair requirement> <desired output format (text/html)> <input filename>
"""

# This uses the -i flag to indicate the input file name, since it might contain spaces
parser = OptionParser()
parser.add_option("-i", "--input", dest="filename")

(options, args) = parser.parse_args()

# The number of base pairs to use to check the beginning of the sequence
bp_match = int(sys.argv[4])



# The desired output file type - text or html
# By default this is 'text' for the command line scripts and 'html' for the cgi script
output_type = sys.argv[5]

# The input file name
# This is required because the input file name in not conserved through the cgi scripts
# file_input = sys.argv[6]




try:
    # Open the CD-HIT clustered file *.clstr
    cluster_file = open(sys.argv[1], 'r')

    # Open the fasta file used as input for CD-HIT
    fasta_file = open(sys.argv[2], 'r')

    # Output files
    outfile = sys.argv[3]

except:
    print """
Your input files cannot be found.
"""
    sys.exit(2)

# Read in the cd-hit clustered file
try:
    cluster_list = cluster_file.read()
except:
    print 'Cannot open', cluster_file, '\n'
    sys.exit(2)

# Read in the FASTA file and check to make sure it's in FASTA format
try:
    fasta_dict_raw = fasta.load(fasta_file)
except:
    print '\n', sys.argv[2], 'does not appear to be a fasta file\n'
    sys.exit(2)


fasta_dict = {}

# Just take everything before the first space in the first line of the FASTA file as 
# the key.  This is also how cd-hit takes the name, so the keys will match.

for fasta_key in fasta_dict_raw:
    new_fasta = fasta_key.split(' ')
    new_key = new_fasta[0]
    fasta_dict[new_key] = fasta_dict_raw[fasta_key]


# Output file for the list of all the sequences in each cluster
n_output = outfile + '.fasta_clusters'

try:
    output = open(n_output, 'wt')
except:
    print 'Cannot open', n_output, 'for writing'
    sys.exit(2)


# Output file for the summary
n_summary = outfile + '.cluster_summary'
# n_plot = outfile + '.plot'

try:
    output_summary = open(n_summary, 'wt')
#    plot_summary = open(n_plot, 'wt')
except:
    print 'Cannot open', n_summary, 'for writing'
    sys.exit(2)


# Output file for plotting the number of clusters of each size
n_clstr_size = outfile + '.cluster_sizes'

try: 
    output_clstr_size = open(n_clstr_size, 'wt')
except:
    print 'Cannot open', n_clstr_size, 'for writing'
    sys.exit(2)


# Output file for a fasta file of unique sequences
unique = outfile + '_unique.fa'
try:
    output_unique = open(unique, 'w')
except:
    print 'Cannot open', unique, 'for writing'
    sys.exit(2)


# Parse the cd-hit *.clstr file, to extract the information about what sequences
# are in each cluster

cluster_db, cluster_db_count, unique_list = cdhit_parse.read_clusters(cluster_list)

cluster_size_db = {}

cluster_keys = sorted(cluster_db_count.iteritems(), key=itemgetter(1), reverse=True)


# **************  Analyze clusters to see if the first 3 bp match ********

# A function for comparing the beginning base pairs and creating appropriate 
# clusters.  Us bp_match as the number of initial base pairs to check.

def compare_bp(cluster_sequences, fasta_dictionary):

    # The template is the first three bp of the first sequence in the cluster
    template_seq = fasta_dictionary[cluster_sequences[0]]
    template = template_seq[0:bp_match]

    # Creating new clusters after checking the first base pairs
    new_cluster_list = []

    # If a sequence doesn't match in the first bps, then create a new cluster
    revised_cluster_list = []

    bp3_dict = {}

    for seq in cluster_sequences:
        bp3 = fasta_dictionary[seq][0:bp_match]
        bp3_dict[seq] = bp3
        
    new_rkey = 0

    for bp3_key in bp3_dict:
        if bp3_dict[bp3_key] == template:
            new_cluster_list.append(bp3_key)
        elif bp3_dict[bp3_key] != template:
            new_rkey = new_rkey + 1
            revised_cluster_list.append(bp3_key)
    
    return new_cluster_list, revised_cluster_list


# For each cluster go through and check to make sure the first bp_match base pairs
# are the same.  If they're not, make new clusters.  If bp_match is set to 0, then 
# clusters will not be affected by comparing the initial base pairs.

new_key = 0
cluster_set = {}

for cluster in cluster_keys:

    new_key = new_key + 1

    cluster_seqs = cluster_db[cluster[0]]

    (new_cluster_list, revised_cluster_list) = compare_bp(cluster_seqs, fasta_dict)
    
    cluster_set[new_key] = new_cluster_list

    # If some clusters have sequences that don't match, make them into their own clusters

    while revised_cluster_list:
        (new_rev_cl, rev_rev_cl) = compare_bp(revised_cluster_list, fasta_dict)
        
        new_key = new_key + 1

        cluster_set[new_key] = new_rev_cl
    
        if rev_rev_cl:
            revised_cluster_list = rev_rev_cl
            continue
        else:
            break


# print 'new clusters', cluster_set

# Create a dictionary with the cluster number as the key and the number of sequences
# in that cluster as the info
# Create a dictionary with the cluster number as the key and the reference sequence 
# as the info

cluster_num_seq = {}
cluster_ref_seq = {}

for j in cluster_set:

    ref_seq = cluster_set[j][0]
    ref_seq_fasta = fasta_dict[cluster_set[j][0]]
    for s in cluster_set[j]:
        if (len(fasta_dict[s]) > len(ref_seq_fasta)):
            ref_seq = s

    cluster_ref_seq[j] = ref_seq
    cluster_num_seq[j] = len(cluster_set[j])

        

# Use this information for the output

# get a list of the clusters in order of most clusters to least clusters
sorted_clusters = sorted(cluster_num_seq.iteritems(), key=itemgetter(1), reverse=True)

# largest = nlargest(10, cluster_num_seq.iteritems(), itemgetter(1))

largest = sorted_clusters[0:10]

num_seq = float(len(fasta_dict))
num_unique = float(len(cluster_set.keys()))
percent_float = (num_seq - num_unique)/num_seq*100
percent = round(percent_float, 2)



# Output to terminal

# If this is coming from the cgi script the output will be HTML.  If it's from 
# the command line script, it will be text.

if (output_type == 'html'):

    print '<dl><dd>Number of reads:', int(num_seq), '<dd>Number of unique reads:', int(num_unique), '<dd>Percent of reads that are replicates:', percent, '%</dl>'

    print '<p>10 largest clusters <dl>'
    for l in largest:
        lkey = l[0]
        print '<dd>Cluster', lkey, 'Number of sequences:', cluster_num_seq[lkey]
    print '</dl>'

elif (output_type == 'text'):
    print 'Number of reads:', int(num_seq), '\nNumber of unique reads:', int(num_unique), '\nPercent of reads that are replicates:', percent, '%\n'

    print '10 largest clusters\n'
    for l in largest:
        lkey = l[0]
        print 'Cluster', lkey, 'Number of sequences:', cluster_num_seq[lkey]


# ***************

# Output to files

# Output summary
output_summary.write('File analyzed: %s\n454 Replicate Filter version 0.3\nNumber of sequences: %s  Number of unique reads: %s  Percent of repeats %s\n' % (options.filename, num_seq, num_unique, percent,))

output_summary.write('Cluster\tRef sequence\tNum of seq\n')


# ++++++++++

# Writing each sequence out, so that you have a file with all the sequences in
# each cluster

# for incrementing for plot output
n = 0

# Right now the output is in order from most sequences in a cluster to least, except
# where clusters are split after the initial base pair check.
# If you want the output in order of most sequences in a cluster to least for all clusters, uncomment
# this and replace - for cluster_id in cluster_set and uncomment - cluster_id = cluster_id_raw[0]

# for cluster_id_raw in sorted_clusters:

output.write('File analyzed: %s' % (options.filename))

for cluster_id in cluster_set:

    n = n + 1

#    cluster_id = cluster_id_raw[0]
    output.write('\n----------------------------------------\nCluster %s   Reference sequence: %s Number of sequences: %s\n' % (cluster_id, cluster_ref_seq[cluster_id], cluster_num_seq[cluster_id],))
    for item in cluster_set[cluster_id]:
        output.write('>%s\n%s\n' % (item, fasta_dict[item],))

# ++++++++++

# Detail for output summary file
    output_summary.write('%s\t%s\t%s\n' % (cluster_id, cluster_ref_seq[cluster_id], cluster_num_seq[cluster_id],))
    
# ~~~~~~  Create a file that can be used for plotting the distribution of cluster numbers
#         for clusters that have more than one sequence

#    key2_num = int(cluster_num_seq[cluster_id])   

#    if key2_num > 1:
#        plot_summary.write('%s\t%s\n' % (n, sorted_clusters[cluster_id][1],))
    
# ~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~  Create a file that has the number of sequences in a cluster versus
#         the number of clusters there are of that size

    # Determine the number of clusters with a given number of reads in it
    if cluster_size_db.has_key(cluster_num_seq[cluster_id]):
        num_clusters = cluster_size_db[cluster_num_seq[cluster_id]]
        num_clusters = num_clusters + 1
        cluster_size_db[cluster_num_seq[cluster_id]] = num_clusters
    else:
        cluster_size_db[cluster_num_seq[cluster_id]] = 1


cluster_size_keys = sorted(cluster_size_db.iteritems(), reverse=False)

output_clstr_size.write('File analyzed:\n%s\nCluster size\tNumber of clusters\n' % (sys.argv[2],))

for clstr_num_temp in cluster_size_keys:
    clstr_num = clstr_num_temp[0]
    output_clstr_size.write('%s\t%s\n' % (clstr_num, cluster_size_db[clstr_num],))

# ~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~ Create a file with all the unique sequences

for q in cluster_ref_seq:
    output_unique.write('>%s\n%s\n' % (cluster_ref_seq[q], fasta_dict[cluster_ref_seq[q]]))

# ~~~~~~~~~~~~~

output.close()
output_summary.close()
# plot_summary.close()
output_clstr_size.close()
output_unique.close()
