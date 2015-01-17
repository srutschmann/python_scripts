#!/usr/bin/python

import sys
import glob

"""
fastaConcatAlignments.py

Prints the combined output on stdout and the additional information on stderr.

That is,

 $ python fastaConcatAlignments.py >output.fasta  2>additional_output.txt

will send the output on stdout to the "output.fasta" file (overwriting it) and
the additional output on stderr to the "additional_output.txt" file (also
overwriting it).

To invoke, either specify the names of the input files on the commandline, or
run without any arguments - in which case it will pick up all the *fasta files
(in sorted order) in the current directory.

For example:

  $ python fastaConcatAlignments.py cox3_1_something.fasta atp6_1_something.fasta >op.fasta

will process the two fasta files specified in the specified order. The output
will go to op.fasta. The additional information, since it is not redirected
anywhere, will be printed on the terminal.

If the current directory only has those two fasta files, running the program
without any arguments will have the same effect as above:

  $ python fastaConcatAlignments.py >op.fasta

HOWEVER, the files will be processed in the sorted order by their names
- which is the opposite order from that in the previous example.

Since it processes all the *.fasta files in the current directory, one should
be careful when redirecting the output to a .fasta file in the same directory
as it might be processed as part of the next program run as an input file.

SUGGESTION: name the output file somethig else. eg - output.fasta.out
"""

def die(msg):
    print >>sys.stderr, msg
    sys.exit(1)

def parsefile(fileobj):
    """
    a simple state machine parser which yields (tag, sequence) pairs for the
    given filename
    """
    tag = fileobj.readline().strip()
    sequence = []
    assert tag.startswith('>'), "was expecting the first line in the file %s to be a tag" % (fileobj.name,)
    for line in map(str.strip, fileobj.readlines()):
        if line.startswith('>'):
            yield (tag, "".join(sequence))
            tag = line
            sequence = []
        else:
            sequence.append(line)
    # final tag/sequence
    yield (tag, "".join(sequence))


def sorted(lst):
    lst.sort()
    return lst


def main():
    filenames = sys.argv[1:] or sorted(glob.glob("*.fasta"))
    if not filenames:
        die("no files to process")

    # convert all files to dicts with tag => sequence pairs
    filedicts = map(lambda f: dict(parsefile(file(f))), filenames)

    # get a sorted list of all tags across all files
    all_tags = set()
    map(lambda fd: all_tags.update(fd.keys()), filedicts)
    all_tags = sorted(list(all_tags))

    # for each file, find the maximum sequence length
    def max_seq_length(fd):
        return max(map(str.__len__, fd.values()))

    seq_lengths = map(max_seq_length, filedicts)

    # now start writing
    def pad(s, length):
        assert len(s) <= length, "string '%s' is larger than the length after padding" % (s,)
        return s + ('-'* (length - len(s)))

    def tag_string(tag):
        return "".join(map(lambda (fd, seqlen): pad(fd.get(tag, ''), seqlen), zip(filedicts, seq_lengths)))

    for tag in all_tags:
        print tag
        print tag_string(tag)

    # now print the additional output on stderr
    # filename: <length>; <position_range>;
    # Hmm - it seems that all sequences in a particular file are of the same
    # length.
    # If so, the code above which calculates the maximum sequence length for
    # each tag could be much simpler as it just needs to examine a single
    # sequence from each file.
    # Moreover, the maximum sequence length would be a single number (and not
    # per-tag) as it is now.
    # In any case, the current solution is more general than I think what was
    # needed, so I'll let it be (and it works).
    pos = 1
    for filename, fd, seqlen in zip(filenames, filedicts, seq_lengths):
        print >>sys.stderr, "%s: %d; %d-%d" % (filename, seqlen, pos, pos+seqlen-1)
        pos = pos+seqlen

if __name__ == "__main__":
    main()


