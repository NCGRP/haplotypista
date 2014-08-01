haplotypista generates a new data matrix by combining a given number of SNPs into a single haplotype block.
	Each unique haplotype block is given a new, integer-coded, allelic state.  Also calculates block length.
	
Input file format:
First row--space delimited list assigning SNP to chromosome
Second row--position of SNP along chromosome
Third + rows --individual genotypes, sample name is first column

Example input file with 9 SNPs (Unix line breaks):
1 1 1 2 2 3 3 3 3
42 456 6032 79 876 3 24 53 657
a 0 1 0 0 0 1 1 1 0
b 1 1 0 0 0 1 0 1 0
c 0 1 1 1 0 1 0 1 0

To compile:  use "make"
Usage: haplotypista -i inputfile -o outputfile -l logfile -h blocklength 
where, 
blocklength = length of haplotype block in number of adjacent SNPs to be combined

Examples: haplotypista -i hin.txt -o hout.txt -l hlog.txt -h 2
