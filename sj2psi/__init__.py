import pandas as pd

__version__ = '0.2.2'

COLUMN_NAMES = ('chrom', 'intron_start', 'intron_stop', 'strand',
                'intron_motif', 'annotated',
                'unique_junction_reads', 'multimap_junction_reads',
                'max_overhang')
NEG_STRAND_INTRON_MOTIF = {'CT/AC': 'GT/AG',
                           'CT/GC': 'GC/AG',
                           'GT/AT': 'AT/AC',
                           'non-canonical': 'non-canonical'}

def int_to_intron_motif(n):
    if n == 0:
        return 'non-canonical'
    if n == 1:
        return 'GT/AG'
    if n == 2:
        return 'CT/AC'
    if n == 3:
        return 'GC/AG'
    if n == 4:
        return 'CT/GC'
    if n == 5:
        return 'AT/AC'
    if n == 6:
        return 'GT/AT'


def read_sj_out_tab(filename):
    """Read an SJ.out.tab file as produced by the RNA-STAR aligner into a
    pandas Dataframe

    Parameters
    ----------
    filename : str of filename or file handle
        Filename of the SJ.out.tab file you want to read in

    Returns
    -------
    sj : pandas.DataFrame
        Dataframe of splice junctions with the columns,
        ('chrom', 'intron_start', 'intron_stop', 'strand',
        'intron_motif', 'annotated', 'unique_junction_reads',
        'multimap_junction_reads', 'max_overhang')

    """
    sj = pd.read_table(filename, header=None, names=COLUMN_NAMES, sep='\s+')
    sj.intron_motif = sj.intron_motif.map(int_to_intron_motif)

    # Convert integer strand to symbol
    # Use index-based replacement because it's 100x faster
    rows = sj.strand == 1
    sj.loc[rows, 'strand'] = '+'
    rows = sj.strand == 2
    sj.loc[rows, 'strand'] = '-'

    # Translate negative strand intron motifs
    rows = sj.strand == '-'
    sj.loc[rows, 'intron_motif'] = sj.intron_motif[rows].map(
        lambda x: NEG_STRAND_INTRON_MOTIF[x])
    sj.annotated = sj.annotated.astype(bool)

    # Add intron location
    sj['intron_location'] = sj.chrom.astype(str) + ':' \
        + sj.intron_start.astype(str) + '-' \
        + sj.intron_stop.astype(str) + ':' \
        + sj.strand.astype(str)

    return sj


def chr_start_stop_to_sj_ind(chr_start_stop, sj):
    """Transform a 'chr1:100-200' string into index range of sj dataframe

    Parameters
    ----------
    chr_start_stop : str
        Genome location string of the format chr:start-stop
    sj : pandas.DataFrame
        Dataframe of splice junctions as created by read_sj_out_tab

    Returns
    -------
    ind : pandas.Series (bool)
        Boolean series which can be used to index the sj

    """
    chrom, startstop = chr_start_stop.replace(',', '').split(':')
    start, stop = map(int, startstop.split('-'))
    return (sj.chrom == chrom) & (start < sj.intron_start) \
        & (sj.intron_stop < stop)

def _full_index(sample_id_col):
    """Return all columns necessary to uniquely specify a junction"""
    return [sample_id_col, 'chrom', 'intron_start', 'intron_stop', 'strand']


def add_possible_donors_acceptors(sj, sample_id_col='sample_id',
                                  reads_col='total_junction_reads'):
    """Add other observed end (start) sites for shared junction starts (ends)

    For samples that share a donor (acceptor), add the acceptor (donor) to
    other samples, with zero reads. This ***explicitly*** adds all possible
    `donor, acceptor` pairs to each sample, so the other samples will properly
    get $\Psi$ scores of 1 for the `donor, acceptor` pairs they have not seen.

    **From:**

        Sample 1:
        chr1:100-200: 25 reads
        [   ]-------------[   ]
        Sample 2:
        chr1:100:250: 40 reads
        [   ]----------------------[    ]

    **To:**

        Sample 1:
        chr1:100-200: 25 reads
        [   ]-------------[   ]
        chr1:100:250: 0 reads
        [   ]----------------------[    ]
        Sample 2:
        chr1:100-200: 0 reads
        [   ]-------------[   ]
        chr1:100:250: 40 reads
        [   ]----------------------[    ]

    A one-line summary that does not use variable names or the
    function name.

    Several sentences providing an extended description. Refer to
    variables using back-ticks, e.g. `var`.

    Parameters
    ----------
    sj : pandas.DataFrame
        A table of splice junctions, with the required column names:
        "chrom", "intron_start", "intron_stop", "strand" and the
        `sample_id_col` and `reads_col` specified below
    sample_id_col : str, optional
        Name of the column containing sample ids. Default: "sample_id"
    reads_col : str, optional
        Name of the column containing read counts.
        Default: "total_junction_reads"

    Returns
    -------
    sj_appended : pandas.DataFrame
        A table of splice junctions with all observed splicing donor and
        acceptors for each sample in the dataset

    Note
    ----
    This function takes a long time and will use a lot of memory. Please plan
    accordingly.
    """
    dfs = []

    full_index = _full_index(sample_id_col)

    for intron in ('intron_start', 'intron_stop'):
        for name, df in sj.groupby(['chrom', intron, 'strand']):
            df = df.set_index(full_index)
            df = df[reads_col].unstack(
                [sample_id_col, 'chrom', intron, 'strand'])
            df = df.fillna(0)
            df = df.stack(
                [sample_id_col, 'chrom', intron, 'strand'])
            df = df.reset_index()
            df = df.sort(sample_id_col)
            df = df.rename(columns={0: reads_col})
            dfs.append(df)
    sj_appended = pd.concat(dfs, ignore_index=True)

    # Observed junctions appear twice, so remove them
    sj_appended = sj_appended.drop_duplicates()

    # Reorder columns
    sj_appended = sj_appended.reindex(columns=full_index + [reads_col])

    # Sort by sample id, chrom, start, stop, strand
    sj_appended = sj_appended.sort(full_index)
    return sj_appended


def calculate_psis(sj, sample_id_col='sample_id',
                   reads_col='total_junction_reads'):
    """Get percent spliced-in (psi) scores for 3' and 5' splice sites

    This method of getting psi scores is 5x faster than using
    `lambda x: x/x.sum()`

    Parameters
    ----------
    sj : pandas.DataFrame
        A table of splice junctions, with the required column names:
        "chrom", "intron_start", "intron_stop", "strand" and the
        `sample_id_col` and `reads_col` specified below
    sample_id_col : str, optional
        Name of the column containing sample ids. Default: "sample_id"
    reads_col : str, optional
        Name of the column containing read counts.
        Default: "total_junction_reads"

    Returns
    -------
    sj_psis : pandas.DataFrame
        A table of splice junctions with Psi scores calculated for both 3'
        and 5' splice sites
    """
    splice_sites = {r'$\Psi_5$': 'intron_start', r'$\Psi_3$': 'intron_stop'}

    full_index = _full_index(sample_id_col)

    for psi, ss in splice_sites.items():
        sj = sj.sort(full_index)
        groupby = ['sample_id', 'chrom', ss, 'strand']
        reads = sj.groupby(groupby)[reads_col].sum()
        reads.name = reads_col

        # sort_index() is necessary for the rows of sj_appended
        # to be exactly the right order and match with psi_scores
        sj = sj.set_index(groupby).sort_index()
        psi_scores = sj[reads_col].divide(reads).fillna(0)
        sj[psi] = psi_scores.values
        sj = sj.reset_index()
    return sj


def get_psis(sj, min_unique=5, min_multimap=10):
    """Calculate Percent spliced-in (Psi) scores of each junction

    As described in Pervouchine et al, Bioinformatics (2013)
    [doi: 10.1093/bioinformatics/bts678], we will take the approach of asking,
    how often is this donor site (5' splice site) used with this acceptor
    site (3' splice site), compared to ALL OTHER acceptors?

    Same goes for acceptor sites. How often is this acceptor site, used with
    this donor site, compared to ALL OTHER donors?

    To illustrate, check out this example. Each "-" represents 10 bp

    Splice junction fig     genome location     number of reads
    [  ]--------[    ]        chr1:100-180        90
    [  ]----------[  ]        chr1:100-200        10
    [     ]-------[  ]        chr1:130-200        40

    For the 5' splice site chr1:100, we have 90+10 = 100 total reads. Thus the
    "psi5" for chr1:100-180 is 90/100 = 0.9, and 0.1 for chr:100-200.

    For the 3' splice site chr1:200, we have 10+40 = 50 total reads. Thus the
    "psi3" for chr1:100-200 is 10/50 = 0.2, and 0.8 for chr:130-200.

    What's left is the uninteresting splice sites of chr1:180 and chr1:130,
    both of which didn't have any variance and were always used. Thus psi3
    for chr1:180 is 1.0, and psi5 for chr1:130 is 1.0 as well.

    Parameters
    ----------
    sj : pandas.DataFrame
        A splice junction dataframe as created by read_sj_out_tab, specifically
        with the columns,
        ('chrom', 'intron_start', 'intron_stop', 'strand',
        'intron_motif', 'annotated', 'unique_junction_reads',
        'multimap_junction_reads', 'max_overhang')
    min_unique : int, optional
        Minimum number of unique reads per junction. Default 5.
    min_multimap : int, optional
        Minimum number of multimapping reads per junction. Default 10

    Returns
    -------
    sj_with_psi : pandas.DataFrame
        The original dataframe, now with the columns psi5 and psi3 for
        percent spliced-in scores of each junction.

    >>> import pandas as pd
    >>> data = {'chrom': ['chr1', 'chr1', 'chr1'],
    ... 'intron_start':[100, 100, 130], 'intron_stop':[100, 200, 200],
    ... 'unique_junction_reads':[90, 10, 40],
    ... 'multimap_junction_reads':[0, 0, 0]}
    >>> sj = pd.DataFrame(data)
    >>> get_psis(sj) # doctest: +NORMALIZE_WHITESPACE
      chrom  intron_start  intron_stop  multimap_junction_reads  \\
    0  chr1              100             100                        0
    1  chr1              100             200                        0
    2  chr1              130             200                        0
    <BLANKLINE>
       unique_junction_reads  multimap_junction_reads_filtered  \\
    0                     90                                 0
    1                     10                                 0
    2                     40                                 0
    <BLANKLINE>
       unique_junction_reads_filtered  total_filtered_reads  psi5  psi3
    0                              90                    90   0.9   1.0
    1                              10                    10   0.1   0.2
    2                              40                    40   1.0   0.8
    <BLANKLINE>
    [3 rows x 10 columns]
    """
    sj['multimap_junction_reads_filtered'] = sj.multimap_junction_reads[
        sj.multimap_junction_reads >= min_multimap]
    sj['unique_junction_reads_filtered'] = sj.unique_junction_reads[
        sj.unique_junction_reads >= min_unique]
    sj['total_filtered_reads'] = sj.multimap_junction_reads_filtered.add(
        sj.unique_junction_reads_filtered)
    sj.total_filtered_reads = sj.total_filtered_reads.astype('float')

    # Calculate psi scores as in Pervouchine et al, Bioinformatics (2013)
    # doi: 10.1093/bioinformatics/bts678
    psi5_groupby = ['chrom', 'intron_start']
    psi3_groupby = ['chrom', 'intron_stop']

    groupbys = {'psi5': psi5_groupby, 'psi3': psi3_groupby}
    for name, groupby in groupbys.items():
        denominator = '{0}_denominator'.format(name)
        s = sj.groupby(groupby).total_filtered_reads.sum()
        s.name = denominator
        sj.set_index(groupby, inplace=True, drop=False)
        sj = sj.join(s)
        sj[name] = sj.total_filtered_reads / sj[denominator]
        sj.reset_index(inplace=True, drop=True)

    return sj
