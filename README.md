primerDesign
============

A tool to design highly specific PCR primers for the validation of genomic alterations, including structural variants. Supported variant types are: SNVs/SNPs, deletions, tandem duplications, inversion, and translocations.

The tool was applied in the 1000 genomes project paper "An integrated map of structural variation in 2,504 human genomes" by Sudmant et al. (Sep 2015; http://dx.doi.org/10.1038/nature15394).


Prerequisites
-------------

To run primerDesign, you need the following software/tools:
* Primer3 (http://sourceforge.net/projects/primer3/files/primer3)
* NCBI BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST)
* Python 2.7 with packages:
  * xlrd (https://pypi.python.org/pypi/xlrd)
  * xlwt (https://pypi.python.org/pypi/xlwt)


Furthermore, you need to generate a BLAST database for your reference genome, e.g.:
`makeblastdb -in REF_GENOME.fa -dbtype nucl -parse_seqids`


NB: primerDesign has only been tested on Linux systems.


Input
-----

#### Configuration file

In the configuration file, you specify the location of the reference genome, of the Primer3 paramter files, and of the executables primer3_core, as blastn, as well as blastdbcmd. Furthermore, you can define the paraemters used for the primer design.


#### List of genomic variants

The list of genomic variants to design PCR primers for can be provided as an Excel xls file or as a tab-delimited text file. The columns are:
* chromosome of variant start
* position of variant start
* chromosome of variant end (only necessary for translocations, for all other variants you can put ".")
* position of variant end
* variant type 
* comment (if there is no comment, put ".")

Variant type can be one of the following
* 'snv' (single nucleotide variant/polymorphism)
* 'del' (deletion)
* 'dup' (tandem duplication)
* 'inv3to3', 'inv5to5', 'invAlt', 'invRef' (inversion with indication of orientation. For 'invAlt' and 'invRef', one primer pair for each inversion breakpoint is designed)
* 'trans3to3', 'trans3to5', 'trans5to3', 'trans5to5' (translocation with indication of orientation)


Running primerDesign
----------------------

`pyhon primerDesign.py CONFIG_FILE VARIANT_LIST [RESULT_FILE]`

If no result file is specified, the results will be printed to stdout.


Output
------

The pipeline produces two list of primer pairs. The first one contains primer pairs where at least one of the two primers is unique. And in the second list, both primers of a pair are unique (however, they might be further away from the break point or of lower quality).

The columns of the output are the same as the input, plus:
Primer3 quality, primer1, primer2, uniqueness primer1, uniqueness primer2, temp primer1, temp primer2, distance to breakpoint primer1, distance to breakpoint primer2, primer1 sequencable, primer2 sequencable, primer1 coordinate, primer2 coordinate, expected band size variant, expected band size reference.



