import string, re, sys, cStringIO
from array import array

__printdebug__ = 0

legal_dna = "ACGTNX-"
legal_protein = "ABCDEFGHIKLMNPQRSTUVWXYZ-"

complementTranslation = string.maketrans('ACTG', 'TGAC')

#
# load
#

def load(f, strict=0):
    """
    Loads sequences in FASTA format from the given file object.
    Returns dict.

    The 'strict' option forces the file to be in a stricter FASTA format,
    i.e. it must contain ONLY sequences, and they must be ONLY
    DNA or Protein sequences.
    """
    d= {}

    # read name
    l = f.readline()

    # check for Mac EOL characters.  Note that this assumes that the first
    # line we're getting is representative of the file, which will usually
    # be the case, because 'readline' seems to pay attention to platform-
    # specific EOLs.
    #
    # If we do find a Mac EOL (^M, character # 13) then reconstruct the
    # file into a cStringIO object.  Unless people are dealing with > 20 MB
    # of sequence, that should work fine ;).
    #
    macEOL = '%c' % (13,)
    if macEOL in l:
        l = l + f.read()                # get all of the file
        
        l = l.replace(macEOL, "\n")     # replace EOL characters with \n
        
        f = cStringIO.StringIO(l)       # replace file with StringIO

    #
    # Read through the file and pick out the good lines.  Ignore
    # non-FASTA gobbledygook at the top.
    #

    if not strict:
        while l and l[0] != '>':
            l = f.readline()

    if not l:
        raise "Error! no beginning '>'" 

    try:
        assert(l[0] == '>')
    except:
        errorMsg = """

ERROR!  This doesn't look like a fasta file...

The first several lines are:
--------------------
%s%s%s%s%s
...
--------------------
</pre>

""" % (l, f.readline(), f.readline(), f.readline(), f.readline())
        raise errorMsg

    nextName = l
    sequence = ""

    while nextName:
        name = string.rstrip(nextName[1:])
        while 1:
            l = f.readline()
            if not l:                   # end of file
                nextName = None
                break
            elif l[0] == '>':           # end of sequence
                nextName = l
                break

            sequence = sequence + string.rstrip(l)
        # end while 1

        sequence = string.upper(sequence)
        sequence = string.join(string.split(sequence), '')

        if len(name) and len(sequence):
            if strict and (is_protein(sequence) or is_dna(sequence)):
                d[name] = sequence
            elif not strict:
                d[name] = sequence
            else:
                raise "Error, sequence contains illegal characters."
            
            sequence = ""
            
    # end while name

    return d
# end load

#
# is_protein
#

def is_protein(seq):
    """
    Returns 1 if it contains only legal values for a protein sequence.

    c.f.  http://www.ncbi.nlm.nih.gov/BLAST/fasta.html
    """

    if is_dna(seq):
        return 0
    
    for ch in seq:
        if ch not in legal_protein:
            return 0

    return 1
# is_protein

#
# is_dna
#

def is_dna(seq):
    """
    Returns 1 if it contains only legal values for a DNA sequence.

    c.f.  http://www.ncbi.nlm.nih.gov/BLAST/fasta.html
    """
    for ch in seq:
        if ch not in legal_dna:
            return 0

    return 1
# is_dna

#
# force_dna
#

def force_dna(seq):
    """
    Removes all whitespace & then replaces all non-ACGT characters with
    N.
    """
    seq = string.join(string.split(seq))

    seq2 = []
    for ch in seq:
        if ch not in legal_dna:
            ch = 'N'

        seq2.append(ch)

    return string.join(seq2, '')
# force_dna

#
# reverse_complement
#

def reverse_complement(s):
    s = string.upper(s)
    assert is_dna(s), "Your sequence must be DNA!"

    # now do the reverse translation magic, using this odd 'array' stuff:
    reverse = array('c', s)
    reverse.reverse()
    reverse = string.join(reverse, '')
    
    reverse = string.translate(reverse, complementTranslation)

    return reverse

# 
# gc-content
# This calculates the proportion of G's and C's in a DNA sequence
# Output is a floating point number
# 

def gc_content(sequence):
    sequence = string.upper(sequence)
    assert is_dna(sequence), "Your sequence must be DNA"

    gc_count = 0
    seq_length = float(len(sequence))

    for x in range(0,len(sequence)):
        if sequence[x] in ('G', 'C'):
            gc_count = gc_count + 1

#    print "gc_count", gc_count, "*******\t"
    gc_content = gc_count/seq_length
    return gc_content



