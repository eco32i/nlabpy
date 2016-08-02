#
# Simple .fasta parser lifted from earlier versions of scikit-bio
#

def is_empty(line):
    """Returns True empty lines and lines consisting only of whitespace."""
    return (not line) or line.isspace()

def is_fasta_label(line):
    return line.strip().startswith('>')
    

def LabeledRecordFinder(is_label_line, ignore=is_empty):
    """Returns function that returns successive labeled records from file.
    Includes label line in return value. Returns list of relevant lines.
    Skips over any lines for which ignore(line) evaluates True (default is
    to skip empty lines).
    """
    def parser(lines):
        with open(lines, 'r') as lines:
            curr = []
            for l in lines:
                line = l.strip()
                if ignore(line):
                    continue
                # if we find the label, return the previous record
                if is_label_line(line):
                    if curr:
                        yield curr
                        curr = []
                curr.append(line)
            # don't forget to return the last record in the file
            if curr:
                yield curr
    return parser

    
FastaFinder = LabeledRecordFinder(is_fasta_label)


def parse_fasta(infile, finder=FastaFinder):
    """Generator of labels and sequences from a fasta file.
    """

    for rec in finder(infile):
        # first line must be a label line
        if not rec[0].startswith('>'):
            raise ValueError("Found Fasta record without label line: %s" % rec)
        # record must have at least one sequence
        if len(rec) < 2:
            raise ValueError("Found label line without sequences: %s" % rec)

        # remove the label character from the beginning of the label
        label = rec[0][1:].strip()
        seq = ''.join(rec[1:])

        yield label, seq
