from graphlite import connect, V
import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import pytest

@pytest.fixture
def sj_exons():
    data = {'junction_location': ['chr1:76-299:+',  # Exon1alt-Exon2 junction
                                  'chr1:201-299:+',  # Exon1-Exon2 junction
                                  'chr1:201-249:+',  # Exon1-Exon2a3ss junction
                                  'chr1:201-499:+',  # Exon1-Exon3 junction
                                  'chr1:201-799:+',  # Exon1-Exon4 junction
                                  'chr1:401-499:+',  # Exon2-Exon3 junction
                                  'chr1:451-499:+',  # Exon2a5ss-Exon3 junction
                                  'chr1:401-699:+',  # Exon2-Exon4 junction
                                  'chr1:401-699:+',  # Exon3-Exon4 junction
                                  'chr1:401-849:+'],  # Exon3-Exon4alt junction
            # Upstream exon
            'exon_5p': ['exon:chr1:50-75:+',  # Exon1alt-Exon2 junction
                        'exon:chr1:100-200:+',  # Exon1-Exon2 junction
                        'exon:chr1:100-200:+',  # Exon1-Exon2a3ss junction
                        'exon:chr1:100-200:+',  # Exon1-Exon3 junction
                        'exon:chr1:100-200:+',  # Exon1-Exon4 junction
                        'exon:chr1:300-400:+,exon:chr1:250-400:+',
                        # Exon2-Exon3 junction
                        'exon:chr1:300-450:+',  # Exon2a5ss-Exon3 junction
                        'exon:chr1:300-400:+,exon:chr1:250-400:+',
                        # Exon2-Exon4 junction
                        'exon:chr1:500-600:+',  # Exon3-Exon4 junction
                        'exon:chr1:500-600:+'  # Exon3-Exon4alt junction
                        ],
            # Downstream exon
            'exon_3p': ['exon:chr1:300-400:+,exon:chr1:300-450:+',

                        'exon:chr1:300-400:+,exon:chr1:300-450:+',  #
                        'exon:chr1:250-400:+',  # Exon1-Exon2a3ss junction
                        'exon:chr1:500-600:+',  # Exon1-Exon3 junction
                        'exon:chr1:700-800:+',  # Exon1-Exon4 junction
                        'exon:chr1:500-600:+',  # Exon2-Exon3 junction
                        'exon:chr1:500-600:+',  # Exon2a5ss-Exon3 junction
                        'exon:chr1:500-600:+',  # Exon2-Exon4 junction
                        'exon:chr1:500-600:+',  # Exon3-Exon4 junction
                        'exon:chr1:850-900:+',  # Exon3-Exon4alt junction
                        ]
            }


    df = pd.DataFrame(data)
    return df

class TestAnnotator(object):

    @pytest.fixture
    def junction_exons(self, sj_exons):
        from sj2psi.annotate import Annotator
        return Annotator.make_junction_exon_table(sj_exons)

    @pytest.fixture
    def annotator(self, junction_exons):
        from sj2psi.annotate import Annotator
        return Annotator(junction_exons)

    def test_init(self, junction_exons):
        from sj2psi.annotate import Annotator

        test = Annotator(junction_exons)
        pdt.assert_frame_equal(test.junction_exons, junction_exons)
        assert test.db is None
        all_exons = junction_exons.exon.unique()
        all_junctions = junction_exons.junction_location.unique()
        items = np.concatenate([all_exons, all_junctions])
        int_to_item = pd.Series(items)
        item_to_int = pd.Series(dict((v, k) for k, v in
                                     int_to_item.iteritems()))

        pdt.assert_array_equal(test.all_exons, all_exons)
        pdt.assert_array_equal(test.all_junctions, all_junctions)
        pdt.assert_array_equal(test.items, items)
        pdt.assert_dict_equal(test.int_to_item, item_to_int)
        pdt.assert_dict_equal(test.item_to_int, item_to_int)

    def test_from_sj_exons(self, sj_exons, annotator):
        from sj2psi.annotate import Annotator

        test = Annotator.from_sj_exons(sj_exons)



@pytest.fixture
def graph():
    graph = connect(":memory:", graphs=['upstream', 'downstream'])

    items = ['exon:chr1:50-75:+',    # Exon 1 alt
             'exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:250-400:+',  # Exon 2, Alt 3' splice site
             'exon:chr1:300-400:+',  # Exon 2
             'exon:chr1:300-450:+',  # Exon 2, Alt 5' splice site
             'exon:chr1:500-600:+',  # Exon 3
             'exon:chr1:700-800:+',  # Exon 4
             'exon:chr1:850-900:+',  # Exon 4 alt

             'chr1:76-299:+',   # Exon1alt-Exon2 junction
             'chr1:201-299:+',  # Exon1-Exon2 junction
             'chr1:201-249:+',  # Exon1-Exon2a3ss junction
             'chr1:201-499:+',  # Exon1-Exon3 junction
             'chr1:201-799:+',  # Exon1-Exon4 junction
             'chr1:401-499:+',  # Exon2-Exon3 junction
             'chr1:451-499:+',  # Exon2a5ss-Exon3 junction
             'chr1:401-699:+',  # Exon2-Exon4 junction
             'chr1:401-699:+',  # Exon3-Exon4 junction
             'chr1:401-849:+',  # Exon3-Exon4alt junction
        ]
    int_to_item = pd.Series(items)
    item_to_int = pd.Series(dict((v, k) for k, v in int_to_item.iteritems()))

    with graph.transaction() as tr:
        tr.store(V(item_to_int['intron:exon1-exon3']).downstream(item_to_int['exon1']))
        tr.store(V(item_to_int['intron:exon1-exon3']).upstream(item_to_int['exon3']))
        tr.store(V(item_to_int['intron:exon1-exon2']).downstream(item_to_int['exon1']))
        tr.store(V(item_to_int['intron:exon1-exon2']).upstream(item_to_int['exon2']))
        tr.store(V(item_to_int['intron:exon1-exon2']).upstream(item_to_int['exon2a5ss']))
        tr.store(V(item_to_int['intron:exon2-exon3']).downstream(item_to_int['exon2']))
        tr.store(V(item_to_int['intron:exon2-exon3']).upstream(item_to_int['exon3']))
        tr.store(V(item_to_int['intron:exon2a5ss-exon3']).downstream(item_to_int['exon2a5ss']))
        tr.store(V(item_to_int['intron:exon2a5ss-exon3']).upstream(item_to_int['exon3']))

        tr.store(V(item_to_int['exon1']).upstream(item_to_int['intron:exon1-exon3']))
        tr.store(V(item_to_int['exon3']).downstream(item_to_int['intron:exon1-exon3']))
        tr.store(V(item_to_int['exon1']).upstream(item_to_int['intron:exon1-exon2']))
        tr.store(V(item_to_int['exon2']).downstream(item_to_int['intron:exon1-exon2']))
        tr.store(V(item_to_int['exon2a5ss']).downstream(item_to_int['intron:exon1-exon2']))
        tr.store(V(item_to_int['exon2']).upstream(item_to_int['intron:exon2-exon3']))
        tr.store(V(item_to_int['exon3']).downstream(item_to_int['intron:exon2-exon3']))
        tr.store(V(item_to_int['exon2a5ss']).upstream(item_to_int['intron:exon2a5ss-exon3']))
        tr.store(V(item_to_int['exon3']).downstream(item_to_int['intron:exon2a5ss-exon3']))




def test_se(graph):
    true = {('exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:300-400:+',  # Exon 2
             'exon:chr1:500-600:+'):  # Exon 3
                ('chr1:201-299:+',  # Exon1-Exon2 junction
                 'chr1:201-499:+',  # Exon1-Exon3 junction
                 'chr1:401-499:+'),  # Exon2-Exon3 junction
            ('exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:250-400:+',  # Exon 2, Alt 3' splice site
             'exon:chr1:500-600:+'):  # Exon 3
                ('chr1:201-249:+',  # Exon1-Exon2a3ss junction
                 'chr1:201-499:+',  # Exon1-Exon3 junction
                 'chr1:401-499:+'),  # Exon2-Exon3 junction
            ('exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:300-450:+',  # Exon 2, Alt 5' splice site
             'exon:chr1:500-600:+'):  # Exon 3
                ('chr1:201-299:+',  # Exon1-Exon2 junction
                 'chr1:201-499:+',  # Exon1-Exon3 junction
                 'chr1:451-499:+'),  # Exon2a5ss-Exon3 junction
            ('exon:chr1:300-400:+',  # Exon 2
             'exon:chr1:500-600:+',  # Exon 3
             'exon:chr1:700-800:+'):  # Exon 4
                ('chr1:401-499:+',  # Exon2-Exon3 junction
                 'chr1:401-699:+',  # Exon2-Exon4 junction
                 'chr1:401-699:+')}  # Exon3-Exon4 junction


def test_mxe(graph):
    true = {('exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:300-400:+',  # Exon 2
             'exon:chr1:500-600:+',  # Exon 3
             'exon:chr1:700-800:+'):  # Exon 4
                ('chr1:201-299:+',  # Exon1-Exon2 junction
                 'chr1:201-499:+',  # Exon1-Exon3 junction
                 'chr1:401-499:+',  # Exon2-Exon3 junction
                 'chr1:401-699:+'),  # Exon3-Exon4 junction

            ('exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:250-400:+',  # Exon 2, Alt 3' splice site
             'exon:chr1:500-600:+',  # Exon 3
             'exon:chr1:700-800:+'):  # Exon 4
                ('chr1:201-249:+',  # Exon1-Exon2a3ss junction
                 'chr1:201-499:+',  # Exon1-Exon3 junction
                 'chr1:401-499:+',  # Exon2-Exon3 junction
                 'chr1:401-699:+'),  # Exon3-Exon4 junction

            ('exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:300-450:+',  # Exon 2, Alt 5' splice site
             'exon:chr1:500-600:+',  # Exon 3
             'exon:chr1:700-800:+'):  # Exon 4
                ('chr1:201-299:+',  # Exon1-Exon2 junction
                 'chr1:201-499:+',  # Exon1-Exon3 junction
                 'chr1:451-499:+',  # Exon2a5ss-Exon3 junction
                 'chr1:401-699:+')}  # Exon3-Exon4 junction

def test_a5ss(graph):
    true = {('exon:chr1:300-400:+',  # Exon 2
             'exon:chr1:300-450:+',  # Exon 2, Alt 5' splice site
             'exon:chr1:500-600:+'):  # Exon 3
                ('chr1:401-499:+',  # Exon2-Exon3 junction
                 'chr1:451-499:+')}  # Exon2a5ss-Exon3 junction

def test_a3ss(graph):
    true = {('exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:250-400:+',  # Exon 2, Alt 3' splice site
             'exon:chr1:300-400:+'):  # Exon 2
                ('chr1:201-299:+',  # Exon1-Exon2 junction
                 'chr1:201-249:+')}  # Exon1-Exon2a3ss junction



def test_afe(graph):
    true = {('exon:chr1:50-75:+',  # Exon 1 alt
             'exon:chr1:100-200:+',  # Exon 1
             'exon:chr1:250-400:+'):
                ('chr1:76-299:+',   # Exon1alt-Exon2 junction
                 'chr1:201-299:+')}  # Exon1-Exon2 junction

def test_ale(graph):
    true = {('exon:chr1:500-600:+',  # Exon 3
             'exon:chr1:700-800:+',  # Exon 4
             'exon:chr1:850-900:+'):  # Exon 4 alt
                ('chr1:401-699:+',  # Exon3-Exon4 junction
                 'chr1:401-849:+')}  # Exon3-Exon4alt junction
