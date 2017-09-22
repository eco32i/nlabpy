import gzip
import numpy as np
from itertools import zip_longest

# The following was shamelessly lifted from scikit-bio circa v0.2.0

def _ascii_to_phred(s, offset):
    """Convert ascii to Phred quality score with specified ASCII offset."""
    return np.fromstring(s, dtype='|S1').view(np.int8) - offset


def ascii_to_phred33(s):
    """Convert ascii string to Phred quality score with ASCII offset of 33.
    Standard "Sanger" ASCII offset of 33. This is used by Illumina in CASAVA
    versions after 1.8.0, and most other places. Note that internal Illumina
    files still use offset of 64
    """
    return _ascii_to_phred(s, 33)


def ascii_to_phred64(s):
    """Convert ascii string to Phred quality score with ASCII offset of 64.
    Illumina-specific ASCII offset of 64. This is used by Illumina in CASAVA
    versions prior to 1.8.0, and in Illumina internal formats (e.g.,
    export.txt files).
    """
    return _ascii_to_phred(s, 64)


def _drop_id_marker(s):
    """Drop the first character and decode bytes to text"""
    id_ = s[1:]
    try:
        return str(id_.decode('utf-8'))
    except AttributeError:
        return id_


def parse_fastq(data, strict=False, enforce_qual_range=True, phred_offset=None):
    r"""yields label, seq, and qual from a fastq file.
    
    Parameters
    ----------
    data : open file object or str
        An open fastq file (opened in binary mode) or a path to it.
    strict : bool, optional
        Defaults to ``False``. If strict is true a FastqParse error will be
        raised if the seq and qual labels dont' match.
    enforce_qual_range : bool, optional
        Defaults to ``True``. If ``True``, an exception will be raised if a
        quality score outside the range [0, 62] is detected
    phred_offset : {33, 64, None}, optional
        What Phred offset to use when converting qual score symbols to integers
        If ``None`` quality is returned as string.
    Returns
    -------
    label, seq, qual : (str, bytes, np.array)
        yields the label, sequence and quality for each entry
    Examples
    --------
    Assume we have a fastq formatted file with the following contents::
        @seq1
        AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
        +
        ````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
        @seq2
        TATGTATATATAACATATACATATATACATACATA
        +
        ]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb
    We can use the following code:
    >>> from StringIO import StringIO
    >>> from skbio.parse.sequences import parse_fastq
    >>> fastq_f = StringIO('@seq1\n'
    ...                     'AACACCAAACTTCTCCACCACGTGAGCTACAAAAG\n'
    ...                     '+\n'
    ...                     '````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF\n'
    ...                     '@seq2\n'
    ...                     'TATGTATATATAACATATACATATATACATACATA\n'
    ...                     '+\n'
    ...                     ']KZ[PY]_[YY^```ac^\\\`bT``c`\\aT``bbb\n')
    >>> for label, seq, qual in parse_fastq(fastq_f, phred_offset=64):
    ...     print(label)
    ...     print(seq)
    ...     print(qual)
    seq1
    AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
    [32 32 32 32 25 30 20 29 32 29 35 30 35 33 34 35 33 35 35 32 30 12 34 30 35
     35 25 20 28 20 28 25 28 23  6]
    seq2
    TATGTATATATAACATATACATATATACATACATA
    [29 11 26 27 16 25 29 31 27 25 25 30 32 32 32 33 35 30 28 28 32 34 20 32 32
     35 32 28 33 20 32 32 34 34 34]
    """

    if phred_offset == 33:
        phred_f = ascii_to_phred33
    elif phred_offset == 64:
        phred_f = ascii_to_phred64
    else:
        phred_f = None

    with gzip.open(data, mode='rt') as fi:
        iters = [iter(fi)] * 4
        for seqid, seq, qualid, qual in zip_longest(*iters):
            seqid = seqid.strip()
            # If the file simply ended in a blankline, do not error
            if seqid is '':
                continue
            # Error if an incomplete record is found
            # Note: seqid cannot be None, because if all 4 values were None,
            # then the loop condition would be false, and we could not have
            # gotten to this point
            if seq is None or qualid is None or qual is None:
                raise ValueError("Incomplete FASTQ record found at end "
                                      "of file")

            seq = seq.strip()
            qualid = qualid.strip()
            qual = qual.strip()

            seqid = _drop_id_marker(seqid)

            try:
                seq = str(seq.decode("utf-8"))
            except AttributeError:
                pass

            qualid = _drop_id_marker(qualid)
            if strict:
                if seqid != qualid:
                    raise ValueError('ID mismatch: {} != {}'.format(
                        seqid, qualid))

            # bounds based on illumina limits, see:
            # http://nar.oxfordjournals.org/content/38/6/1767/T1.expansion.html
            if phred_f is not None:
                qual = phred_f(qual)
            if (enforce_qual_range
                and (phred_f is not None)
                and ((qual < 0).any() or (qual > 62).any())):
                raise ValueError("Failed qual conversion for seq id: %s. "
                                      "This may be because you passed an "
                                      "incorrect value for phred_offset." %
                                      seqid)

            yield (seqid, seq, qual)

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
