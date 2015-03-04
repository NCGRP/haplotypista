haplotypista generates a new data matrix by combining a given number of SNPs into a single haplotype block.
	Each unique haplotype block is given a new, integer-coded, allelic state.  Also calculates block length.

To compile:  use "make"
Usage: haplotypista -i inputfile -o outputfile -l logfile -b blocklengthstart blocklengthend -m missingdatachar
where, 
-b specifies a range of blocklengths to consider
blocklength = length of haplotype block in number of adjacent SNPs to be combined
-m specifies the missing data character used in the input file
-p specifies the ploidy, 1 = haploid, 2 = diploid, etc.

Examples: ./haplotypista -i hexin.txt -o hexout.txt -l hexlog.txt -b 2 4 -m ? -p 1
          ./haplotypista -i AtExample.txt -o AtExout.txt -l AtExlog.txt -b 5 8 -m ? -p 1
          ./haplotypista -i PopulusExample.txt -o PopExout.txt -l PopExlog.txt -b 1 4 -m ? -p 2


Input file format:
First row--space delimited list assigning SNP to chromosome
Second row--position of SNP along chromosome
Third + rows --individual genotypes, sample name is first column

Example input file with 20 SNPs (Unix line breaks), haploid data:
1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 3
42 456 6032 6142 10054 11529 79 876 1024 1125 12058 3 24 53 657 1001 1200 5654 1000254 1000256
a 0 1 0 0 0 1 1 1 0 0 1 1 1 0 1 0 1 0 1
b 1 1 0 0 0 1 0 1 0 0 ? 0 0 0 1 1 1 0 0
c 0 1 1 1 0 1 0 1 0 1 1 0 0 0 1 0 ? 0 0

Example input file with 10 SNPs, phased diploid data, 2 rows per individual:
1 1 1 1 1 2 2 3 3 3
10003096 10003352 10003834 10004403 10004622 45319876 45754486 9689083 9689915 9737708
ALAA-20-1 C G A T ? G C A C G
ALAA-20-1 A T C T ? G C A C G
ALAA-20-2 A T C T A A C A C G
ALAA-20-2 A T A T A G C A C G
ALAA-20-3 A T C T A G T A C G
ALAA-20-3 A T C T A G C A C G


Output:
Produces a series of output data sets with unique haplotypes recoded as unique integers.
	Output data sets receive the suffix ".bX" where X is the block length, i.e. the number of
	contiguous SNPs used to define the allelic state.  If a string of SNPs contains the 
	missing data character defined with -m, the resulting haplotype is recoded as missing 
	data, receiving the designation "-9999" in the output file.  Row 1 of the output data set 
	indicates the chromosome where the haplotype lies.  Row 2 is the length of the 
	haplotype block, using the units in row 2 of the input file.
Also produces a log file containing summary statistics for the output data sets.  Column
	headers are as follows:
	b = haplotype block length
	chromosome = chromosome number ("0" indicates all chromosomes combined)
	n loci = number of loci used for allele count statistics
	Mean allele count = average number of alleles across newly coded loci
	SD allele count = standard deviation of allele count across newly coded loci
	n haplotype length = number of recoded haplotypes used for haplotype length statistics
	Mean haplotype length = average length of a newly coded locus
	SD haplotype length = standard deviation in length across newly coded loci
