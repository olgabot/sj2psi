import pandas as pd
import pandas.util.testing as pdt
import pytest
from StringIO import StringIO


@pytest.fixture
def example_sj_out_tab(tmpdir):
    s = """chr1    14830   14969   2       2       1       0       1       39
chr1    15039   15795   2       2       1       0       1       10
chr1    135359  135680  2       4       0       0       1       22
chr1    146510  155766  2       2       1       19      20      43
chr1    153544  155766  2       2       0       8       15      41
chr1    155832  164262  2       2       1       61      3       46
chr1    156087  156200  2       2       0       1       14      44
chr1    329977  334128  1       1       1       0       2       14
chr1    569184  569583  2       2       0       0       1       17
chr1    655581  659737  2       2       1       0       2       14
chr1    661725  662046  2       4       0       0       1       22
chr1    668587  671992  1       1       0       0       4       28
"""
    df = pd.read_table(StringIO(s), header=None, sep='\s+')
    filename = '{}/SJ.out.tab'.format(tmpdir)
    df.to_csv(filename, index=False, header=False, sep='\t')
    return filename


@pytest.fixture(params=['chr1:146000-655000',
                        pytest.mark.xfail('chr2:146000-655000')])
def chr_start_stop(request):
    return request.param


@pytest.fixture(params=[(0, 'non-canonical'), (1, 'GT/AG'), (2, 'CT/AC'),
                        (3, 'GC/AG'), (4, 'CT/GC'), (5, 'AT/AC'),
                        (6, 'GT/AT')])
def int_motif(request):
    return request.param


def test_int_to_intron_motif(int_motif):
    from sj2psi import int_to_intron_motif
    n, motif = int_motif
    pdt.assert_equal(int_to_intron_motif(n), motif)


def test_read_sj_out_tab(example_sj_out_tab):
    from sj2psi import read_sj_out_tab, COLUMN_NAMES, int_to_intron_motif

    test_output = read_sj_out_tab(example_sj_out_tab)

    true_output = pd.read_table(example_sj_out_tab, header=None,
                                names=COLUMN_NAMES, sep='\s+')
    true_output.intron_motif = true_output.intron_motif.map(
        int_to_intron_motif)
    true_output.annotated = true_output.annotated.astype(bool)

    pdt.assert_frame_equal(test_output, true_output)


@pytest.fixture(params=[None, 10])
def min_unique(request):
    return request.param


@pytest.fixture(params=[None, 20])
def min_multimap(request):
    return request.param


@pytest.fixture
def sj(example_sj_out_tab):
    from sj2psi import read_sj_out_tab
    return read_sj_out_tab(example_sj_out_tab)


def test_get_psis(sj, min_multimap, min_unique):
    from sj2psi import get_psis

    kwargs = {}
    kwargs['min_unique'] = 10 if min_unique is None else min_unique
    kwargs['min_multimap'] = 10 if min_multimap is None else min_multimap

    test_output = get_psis(sj.copy(), **kwargs)

    true_output = sj.copy()
    true_output['multimap_junction_reads_filtered'] = \
        true_output.multimap_junction_reads[
            true_output.multimap_junction_reads >= kwargs['min_multimap']]
    true_output['unique_junction_reads_filtered'] = \
        true_output.unique_junction_reads[
            true_output.unique_junction_reads >= kwargs['min_unique']]
    true_output['total_filtered_reads'] = \
        true_output.multimap_junction_reads_filtered.add(
            true_output.unique_junction_reads_filtered)
    true_output.total_filtered_reads = \
        true_output.total_filtered_reads.astype('float')

    # Calculate psi scores as in Pervouchine et al, Bioinformatics (2013)
    # doi: 10.1093/bioinformatics/bts678
    psi5_groupby = ['chrom', 'first_bp_intron']
    psi3_groupby = ['chrom', 'last_bp_intron']

    groupbys = {'psi5': psi5_groupby, 'psi3': psi3_groupby}
    for name, groupby in groupbys.iteritems():
        denominator = '{}_denominator'.format(name)
        s = true_output.groupby(groupby).total_filtered_reads.sum()
        s.name = denominator
        true_output.set_index(groupby, inplace=True, drop=False)
        true_output = true_output.join(s)
        true_output[name] = true_output.total_filtered_reads \
            / true_output[denominator]
        true_output.reset_index(inplace=True, drop=True)
    pdt.assert_frame_equal(test_output, true_output)


def test_chr_start_stop_to_sj_ind(chr_start_stop, sj):
    from sj2psi import chr_start_stop_to_sj_ind
    test_output = chr_start_stop_to_sj_ind(chr_start_stop, sj)

    chrom, startstop = chr_start_stop.replace(',', '').split(':')
    start, stop = map(int, startstop.split('-'))
    true_output = (sj.chrom == chrom) & (start < sj.first_bp_intron) \
        & (sj.last_bp_intron < stop)
    pdt.assert_array_equal(test_output, true_output)
