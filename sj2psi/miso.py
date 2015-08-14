import sj2psi
import glob

MINIMUM_READS = 10


def intify_and_make_intron(junction):
    """Convert start and stop strings to ints and shorten the interval

    Interval is shortened to match the introns of SJ.out.tab files from STAR

    Parameters
    ----------
    junction : tuple
        (chrom, start, stop, strand) tuple of strings, e.g.
        ('chr1', '100', '200', '-')

    Returns
    -------
    intron : tuple
        (chrom, start, stop, strand) tuple of string, int, int, string.
        Adds 1 to original start and subtracts 1 from original stop

    >>> intify_and_make_intron(('chr1', '100', '200', '-'))
    ('chr1', 101, 199, '-')
    """
    chrom, start, stop, strand = junction
    start = int(start) + 1
    stop = int(stop) - 1
    return chrom, start, stop, strand

def extract_start_stop(exon):
    """Extract start and stop from a chr1,100-200,- tuple

    Parameters
    ----------
    exon : tuple
        (chrom, startstop, strand) tuple of strings, e.g.
        ('chr1', '100-200', '-')

    Returns
    -------
    extracted : tuple
        (chrom, start, stop, strand) tuple of strings, with start and stop
        extracted

    >>> extract_start_stop(('chr1', '100-200', '-'))
    ('chr1', '100', '200', '-')
    """
    chrom, startstop, strand = exon
    start, stop = startstop.split('-')
    return chrom, start, stop, strand

def fix_ri_exons(exons):
    """Convert RI event names to be consistent with other event types

    MISO RI event names are of the form: 'chr1:100-200:+@chr1:250-350:+', which
    use a '-' in between the start and stop rather than a ':' like SE or MXE
    events. An example SE event is:
    'chr1:100:200:+@chr1:300:400:+@chr1:500:600:+'

    This is a small utility function to extract the start and stop of RI event
    names so they are easier to work with.

    Parameters
    ----------
    exons : tuple
        (chrom, startstop, strand) tuple of strings, e.g.
        ('chr1', '100-200', '-')

    Returns
    -------
    converted : tuple
        A tuple of two (chrom, start, stop, strand) tuples of
        (str, int, int, str), with start and stop extracted

    >>> fix_ri_exons((('chr1', '100-200', '+'), ('chr1', '300-400', '+')))
    (('chr1', '100', '200', '+'), ('chr1', '300', '400', '+'))
    """
    exons = map(extract_start_stop, exons)
    return exons


def validate_event(event_name, sj, splice_type):
    """Validate a MISO SE or MXE splicing event with SJ.out.tab from STAR"""
    exons = map(lambda x: x.split(':'), event_name.split('@'))
    chrom = exons[0][0]
    strand = exons[0][-1]

    if splice_type == 'SE':
        exon1, exon2, exon3 = exons

        if strand == '+':
            # Exon1, exon3 junction (excluded, Psi~0 junction) is
            # the end of exon1 and the start of exon3
            exon1_exon3 = chrom, exon1[2], exon3[1], strand
            exon1_exon2 = chrom, exon1[2], exon2[1], strand
            exon2_exon3 = chrom, exon2[2], exon3[1], strand
        else:
            exon1_exon3 = chrom, exon3[2], exon1[1], strand
            exon1_exon2 = chrom, exon2[2], exon1[1], strand
            exon2_exon3 = chrom, exon3[2], exon2[1], strand

        isoform1_junctions = [exon1_exon3,]
        isoform2_junctions = [exon1_exon2, exon2_exon3]
        invalid_junctions = []

    elif splice_type == 'MXE':
        exon1, exon2, exon3, exon4 = exons

        if strand == '+':
            # Exon1, exon3 junction (excluded, Psi~0 junction) is
            # the end of exon1 and the start of exon3
            exon1_exon3 = chrom, exon1[2], exon3[1], strand
            exon1_exon2 = chrom, exon1[2], exon2[1], strand
            exon2_exon3 = chrom, exon2[2], exon3[1], strand
            exon1_exon4 = chrom, exon1[2], exon4[1], strand
            exon2_exon4 = chrom, exon2[2], exon4[1], strand
            exon3_exon4 = chrom, exon3[2], exon4[1], strand
        else:
            exon1_exon3 = chrom, exon3[2], exon1[1], strand
            exon1_exon2 = chrom, exon2[2], exon1[1], strand
            exon2_exon3 = chrom, exon3[2], exon2[1], strand
            exon1_exon4 = chrom, exon4[2], exon1[1], strand
            exon2_exon4 = chrom, exon4[2], exon2[1], strand
            exon3_exon4 = chrom, exon4[2], exon3[1], strand

        isoform1_junctions = [exon1_exon3, exon3_exon4]
        isoform2_junctions = [exon1_exon2, exon2_exon4]
        invalid_junctions = [exon1_exon4, exon2_exon3]
    elif splice_type == 'RI':
        exon1, exon2 = fix_ri_exons(exons)
#         print exon1, exon2

        if strand == '+':
            exon1_exon2 = chrom, exon1[2], exon2[1], strand
        else:
            exon1_exon2 = chrom, exon2[2], exon1[1], strand
        isoform1_junctions = [exon1_exon2]
        isoform2_junctions = []
        invalid_junctions = []

    isoform1_junctions = map(intify_and_make_intron, isoform1_junctions)
    isoform2_junctions = map(intify_and_make_intron, isoform2_junctions)
    invalid_junctions = map(intify_and_make_intron, invalid_junctions)

    try:
        # At least 10 reads on each and every isoform1 junction
        valid_isoform1 = all(sj.loc[junction].unique_junction_reads
                             >= MINIMUM_READS
                             for junction in isoform1_junctions)
    except KeyError:
        valid_isoform1 = False
    try:
        # At least 10 reads on each and every isoform2 junction
        valid_isoform2 = all(sj.loc[junction].unique_junction_reads
                             >= MINIMUM_READS
                             for junction in isoform2_junctions)
    except KeyError:
        valid_isoform2 = False
    try:
        # Absolutely no reads on "invalid" junctions
        invalid_junctions = any(sj.loc[junction].unique_junction_reads
                                > 0 for junction in invalid_junctions)
    except KeyError:
        invalid_junctions = False

    if not invalid_junctions and (valid_isoform1 or valid_isoform2):
        return True
    else:
        return False
