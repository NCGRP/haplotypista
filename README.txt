Haplotypista generates a new data matrix by combining a given number of SNPs into a single haplotype block.
	Each unique haplotype block is given a new, integer-coded, allelic state.  Also calculates physical
	length of resulting blocks.

To compile:  use "make"
Usage: haplotypista -i inputfile -o outputfile -b minblocklength maxblocklength -m missingdatachar -p ploidy
      [-l logfile ] [-g genomicpositions] [-v popnamemap] [-x minbp maxbp]

Mandatory command line flags:
-i    input file 
-o    output file root
-b    range of haplotype block lengths to consider, where blocklength = number of adjacent
      SNPs to be combined. -b 0 0 combines all SNPs in each fragment into a single haplotype block
-m    missing data character used in the input file
-p    ploidy, 1 = haploid, 2 = diploid, etc.

Optional command line flags:
-l    specify a log file name other than the default
-g    genomic positions to use. Provide a comma-delimited list of the form:
      1.75000:1.1000000,14.8697509:14.8697509
      which specifies bp 75000-1000000 (inclusive) of named fragment 1 and bp 8697509 of named fragment 14.
      Additionally, a line break delimited list of positions may be piped to haplotypista,
      without invoking the -g option.  If both are supplied, the piped list will be used.
-v    write m+ input files (.var and .dat) for recoded matrices, required argument is a path
      to a file containing a line break delimited list that maps sample names to populations.
-x    filter haplotype blocks by physical length, include in output only those blocks >= minbp and
      <= maxbp.


Examples: ./haplotypista -i hexin.txt -o hexout -b 2 4 -m ? -p 1
          ./haplotypista -i AtExample.txt -o AtExout -b 5 8 -m ? -p 1
          ./haplotypista -i AtExample.txt -o AtExout1 -b 5 8 -m ? -p 1 -x 1 10000
          ./haplotypista -i PopulusExample.txt -o PopExout -b 1 4 -m ? -p 2 
          ./haplotypista -i PopulusExample.txt -o PopEx2out -b 1 4 -m ? -p 2 -g 1.75000:1.1000000,14.8697509:14.8697509 -v Poppopid.txt
          echo "2.1:2.10000000" | ./haplotypista -i AtExample.txt -o AtEx2out -b 5 8 -m ? -p 1 -v Atpopid.txt

Input file format:
First row--space delimited list assigning SNP to named fragment
Second row--position of SNP along fragment
Third row--amino acid category of each SNP, non-genic (0), synonymous (1),
  non-synonymous (2), unknown (9)
Fourth + rows --individual genotypes, sample name is first column

Example input file with 20 SNPs, haploid data, all loci non-genic:
1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 3
42 456 6032 6142 10054 11529 79 876 1024 1125 12058 3 24 53 657 1001 1200 5654 1000254 1000256
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
a 0 1 0 0 0 1 1 1 0 0 1 1 1 0 1 0 1 0 1 0
b 1 1 0 0 0 1 0 1 0 0 ? 0 0 0 1 1 1 0 0 0
c 0 1 1 1 0 1 0 1 0 1 1 0 0 0 1 0 ? 0 0 1

Example input file with 10 SNPs, phased diploid data (2 rows per individual):
1 1 1 1 1 2 2 3 3 3
10003096 10003352 10003834 10004403 10004622 45319876 45754486 9689083 9689915 9737708
0 0 2 0 0 1 0 0 1 0
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
	data, receiving the designation "-9" in the output file.  
	Row 1--named fragment where the haplotype lies.
	Row 2--length of the haplotype block, using the units in row 2 of the input file.
	Row 3--midpoint of the haplotype block on the fragment.
	Row 4--counts of non-genic, synonymous, and non-synonymous substitutions. For example,
		for blocklength = 5 in above example, 4:0:1 3:2:0.  In this case the first block 
		contains 4 non-genic SNPs and 1 non-synonymous SNP. Block 2 contains 3 non-genic, 
		2 synonymous SNPs and 0 non-synonymous SNPs.
	Row 5--position of first SNP in the haplotype block
	Row 6--position of last SNP in the haplotype block
	
Also produces a log file containing summary statistics for the output data sets.  Column
	headers are as follows:
	b = haplotype block length ("max" indicates a single haplotype block was computed for all SNPs in each fragment)
	fragment = fragment name ("all" indicates all fragments combined)
	n loci = number of loci used for allele count statistics
	Mean allele count = average number of alleles across recoded loci
	SD allele count = standard deviation of allele count across recoded loci
	n haplotype length = number of recoded haplotypes used for haplotype length statistics
	Mean haplotype length = average length of a recoded locus
	SD haplotype length = standard deviation in length across recoded loci
