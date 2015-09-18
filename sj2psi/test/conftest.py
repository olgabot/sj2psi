import pandas as pd
import pytest
import six

# @pytest.fixture(params=['positive', 'negative'])
# def strand(request):
#     if request.param == 'positive':
#         return '+'
#     elif request.param == 'negative':
#         return '-'
#
# @pytest.fixture()
# def exon_locations(strand):
#     startstop = {'exon{}_alt': (50, 75),    # Exon 1 alt
#                  'exon{}': (100, 200),  # Exon 1
#                  'exon{}_a3ss': (250, 400),  # Exon 2, Alt 3' splice site
#                  'exon{}': (300, 400),  # Exon 2
#                  (300, 450),  # Exon 2, Alt 5' splice site
#                  (500, 600),  # Exon 3
#                  (700, 800),  # Exon 4
#                  (850, 900),  # Exon 4 alt
#     }


@pytest.fixture
def sj_out_tab(tmpdir):
    s = """chr1    76   299   1       2       1       0       1       39
chr1    201   299   1       1       1       0       1       10
chr1    201  249  1       1       0       0       1       22
chr1    201  799  1       1       1       19      20      43
chr1    201  799  1       1       0       8       15      41
chr1    155832  164262  1       1       1       61      3       46
chr1    156087  156200  1       1       0       1       14      44
chr1    329977  334128  1       1       1       0       2       14
chr1    569184  569583  1       1       0       0       1       17
chr1    655581  659737  1       1       1       0       2       14
chr1    661725  662046  1       1       0       0       1       22
chr1    668587  671992  1       1       0       0       4       28
"""
    df = pd.read_table(six.StringIO(s), header=None, sep='\s+')
    filename = '{0}/SJ.out.tab'.format(tmpdir)
    df.to_csv(filename, index=False, header=False, sep='\t')
    return filename