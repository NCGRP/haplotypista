haplotypista generates a new data matrix by combining a given number of SNPs into a single haplotype block.
	Each unique haplotype block is given a new, integer-coded, allelic state.  Also calculates block length.
	
Input file format:
First row--space delimited list assigning SNP to chromosome
Second row--position of SNP along chromosome
Third + rows --individual genotypes, sample name is first column

Example input file with 20 SNPs (Unix line breaks):
1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 3
42 456 6032 6142 10054 11529 79 876 1024 1125 12058 3 24 53 657 1001 1200 5654 1000254 1000256
a 0 1 0 0 0 1 1 1 0 0 1 1 1 0 1 0 1 0 1
b 1 1 0 0 0 1 0 1 0 0 1 0 0 0 1 1 1 0 0
c 0 1 1 1 0 1 0 1 0 1 1 0 0 0 1 0 1 0 0

To compile:  use "make"
Usage: haplotypista -i inputfile -o outputfile -l logfile -b blocklengthstart blocklengthend 
where, 
-b specifies a range of blocklengths to consider
blocklength = length of haplotype block in number of adjacent SNPs to be combined
-m specifies the missing data character

Examples: ./haplotypista -i hin.txt -o hout.txt -l hlog.txt -b 2 4 -m ?

