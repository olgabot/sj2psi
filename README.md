sj2psi
======

[[![Build Status](https://travis-ci.org/olgabot/sj2psi.svg)](https://travis-ci.org/olgabot/sj2psi)](https://travis-ci.org/olgabot/sj2psi)[[![Coverage Status](https://coveralls.io/repos/olgabot/sj2psi/badge.svg)](https://coveralls.io/r/olgabot/sj2psi)
](https://coveralls.io/r/olgabot/sj2psi)

DOI: 10.5281/zenodo.9885

Annotation-free estimation of percent spliced in of a junction. This
will convert [RNA-STAR aligner](http://bioinformatics.oxfordjournals.org/content/29/1/15.long) "SJ.out.tab" files to "Percent spliced-in"
(Psi)
scores. Here's an example of an `SJ.out.tab` file:

```
chr1    30040   30563   1       1       1       0       131     45
chr1    30668   30975   1       1       1       0       123     46
chr1    146510  155766  2       2       1       50      92      46
chr1    155832  694346  2       2       0       6       14      26
chr1    317782  322037  1       1       1       0       2       3
chr1    320939  321031  1       1       1       2       4       31
chr1    322229  324287  1       1       1       0       3       21
chr1    322229  324438  1       1       0       0       5       40
chr1    324346  324438  1       1       1       0       2       13
chr1    324711  325802  2       4       0       0       1       30
chr1    663813  664904  2       4       0       0       1       38
chr1    665185  667396  2       2       0       0       2       37
chr1    665185  670802  2       2       0       0       4       40
chr1    667588  682074  2       2       0       0       44      43
chr1    668594  670802  2       2       0       0       4       40
chr1    670994  682074  2       2       0       0       47      43
```

As described in [Pervouchine et al, Bioinformatics (2013)](http://bioinformatics.oxfordjournals.org/content/29/2/273.long), we will take the approach of asking,
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

    >>> import pandas as pd
    >>> data = {'chrom': ['chr1', 'chr1', 'chr1'],
    ... 'first_bp_intron':[100, 100, 130], 'last_bp_intron':[180, 200, 200],
    ... 'unique_junction_reads':[90, 10, 40],
    ... 'multimap_junction_reads':[0, 0, 0]}
    >>> sj = pd.DataFrame(data)
    >>> get_psis(sj, min_multimap=0)
      chrom  first_bp_intron  last_bp_intron  multimap_junction_reads  \
    0  chr1              100             180                        0
    1  chr1              100             200                        0
    2  chr1              130             200                        0
    <BLANKLINE>
       unique_junction_reads  multimap_junction_reads_filtered  \
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
