#include "h.hpp"


/*
To compile:  use "make"
Usage: haplotypista -i inputfile -o outputfile -l logfile -b blocklengthstart blocklengthend 
where, 
-b specifies a range of blocklengths to consider
blocklength = length of haplotype block in number of adjacent SNPs to be combined
-m specifies the missing data character

Examples: ./haplotypista -i hin.txt -o hout.txt -l hlog.txt -b 2 4 -m ?
*/



/***************FUNCTIONS*****************/


//split a string on whitespace
vector<string> split(string const &input) 
{ 
    istringstream buffer(input);
    vector<std::string> ret;

    copy(std::istream_iterator<string>(buffer), 
              std::istream_iterator<string>(),
              back_inserter(ret));
    return ret;
}

vector<unsigned long long> MyVecToULL(vector<std::string> pos, unsigned long long& maxpos)
{
	vector<unsigned long long> ULLvec(pos.size());
	for (unsigned long i=0;i<pos.size();++i) 
	{
		ULLvec[i] = strtoull(pos[i].c_str(), NULL, 10);
		if (ULLvec[i] > maxpos) maxpos = ULLvec[i]; //calculate the maximum position value while going thru everything, update as reference
	}
	return ULLvec;
}
			
int MyGetAlleleCount(unsigned long long h, vector<vector<int> > hapvecint)
{
	int curr = 0;
	int max = curr;
	for (unsigned int j=0;j<hapvecint.size();++j) //indiv
	{
		curr = hapvecint[j][h]; //get allele name for all individuals for current haplotype, the max (+1) will be the number of alleles
		if (curr > max) max = curr;
	}
	return (max + 1); //number of alleles = max + 1 (since 0 is an allele)
}

double MyMean(vector<int> v)
{
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double m = sum / v.size();
	return m;
}

double MyStdev(double m, vector<int> v)
{
	std::vector<double> diff(v.size());
	std::transform(v.begin(), v.end(), diff.begin(),
				   std::bind2nd(std::minus<double>(), m));
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double stdev = sqrt( sq_sum / (v.size() - 1) );
	return stdev;
}

double MyMeanULL(vector<unsigned long long> v)
{
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double m = sum / v.size();
	return m;
}

double MyStdevULL(double m, vector<unsigned long long> v)
{
	std::vector<double> diff(v.size());
	std::transform(v.begin(), v.end(), diff.begin(),
				   std::bind2nd(std::minus<double>(), m));
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double stdev = sqrt( sq_sum / (v.size() - 1) );
	return stdev;
}

//quickly reads a large .dat file into a memory buffer
char * MyBigRead(char* DatFilePath)
{
	FILE * pFile;
	unsigned long long lSize;
	char * buffer;
	size_t result;

	pFile = fopen ( DatFilePath , "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);

	// allocate memory to contain the whole file:
	buffer = (char*) malloc (sizeof(char)*lSize);
	if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

	// copy the file into the buffer:
	result = fread (buffer,1,lSize,pFile);
	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

	fclose (pFile);
	return buffer;
}

vector<vector<std::string> > MyReadInfile(char* InFilePath)
{

	//start the clock
	cout << "Reading data...\n";
	time_t startm,endm;
	time (&startm);

	//read the whole file into a buffer using fread
	char * buffer;
	buffer = MyBigRead(InFilePath);
	stringstream s(buffer); //put giant char array into a stream
	
	strcpy(buffer,""); //clear char*
	
	//read buffer into a vector, one line per item
	vector<std::string> bufvec;
	std::string foo;
	while (getline(s, foo))//get line from s, put in foo, consecutively
	{
		bufvec.push_back(foo);  
	}
	
	s.str(""); //clear stringstream
		
	//sort vector so that individuals from the same population form consecutive elements
	//std::sort(bufvec.begin(), bufvec.end()); //no need to use fancy sort, lexicographic should be fine
		
	//split lines of bufvec into 2d vector
	vector<vector<std::string> > bufvec2d(bufvec.size()); 
	for (unsigned long i=0;i<bufvec.size();++i) bufvec2d[i] = split(bufvec[i]);
	vector<std::string>().swap(bufvec); //clear bufvec
	
	//stop the clock
	time (&endm);
	double dif = difftime (endm,startm);
	if (dif==1) cout << "  " << dif << " second.\n";	
	else cout << "  " << dif << " seconds.\n";
	
	return bufvec2d;
}



//convert alleles to integer coding to save memory, vector access order 2
vector<vector<int> > MyRecode(vector<vector<std::string> > hapvec, std::string missingchar)
{
	cout << "Recoding data...\n";
	time_t startm, endm;
	time (&startm);

	//initialize hapvecint, the recoded data
	vector<vector<int> > hapvecint( hapvec.size(), vector<int>((hapvec[0].size())-1) ); //declare and size vector to hold new integer coded alleles
	
	//unsigned int iz = ColKeyToAllAlleleByPopList.size();
	for (unsigned int i=1;i<hapvec[0].size();++i) //go thru each locus
	{
		vector<std::string> AllelesEncountered; //will contain the unique set of alleles at the locus
		for (unsigned int k=0;k<hapvec.size();++k) //go thru all individuals
		{
			int ColIndex = i;
			std::string a = hapvec[k][ColIndex];
			
			//search for missing substring
			std::size_t found = a.find(missingchar);
			if (found!=std::string::npos) hapvecint[k][ColIndex-1] = -9999; //add the missing data value
			//if (a == "9999") hapvecint[k][ColIndex] = -9999; //add the missing data value
			else
			{
				int AlleleInt; //the new, integerized, name of the allele
				std::vector<std::string>::iterator itr = std::find(AllelesEncountered.begin(), AllelesEncountered.end(), a);
				if (itr != AllelesEncountered.end()) //the allele has been found before
				{
			//cout << "old allele a="<<a<<"\n"; 

					AlleleInt = itr - AllelesEncountered.begin(); //convert itr to index, the index is the integerized allele name
					hapvecint[k][ColIndex-1] = AlleleInt; //add the new name
				}
				else // you have a new allele
				{
			//cout << "new allele a="<<a<<"\n"; 
					AllelesEncountered.push_back(a); //add new allele to list of those encountered
					AlleleInt = AllelesEncountered.size() - 1;  //calculate integerized allele name, starts at 0
					hapvecint[k][ColIndex-1] = AlleleInt; //ColIndex-1 since hapvecint has no column for indiv id's
				}
			}
		}
	}
	//stop the clock
	time (&endm);
	double dif = difftime (endm,startm);
	if (dif==1) cout << "  " << dif << " second.\n";	
	else cout << "  " << dif << " seconds.\n";

	return hapvecint;
}
	
bool fileExists(const char *fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

/***************MAIN*****************/

int main( int argc, char* argv[] )
{
	//parse the command line for options
	char* InFilePath;
	char* OutFilePath;
	char* LogFilePath;
	unsigned long bstart;
	unsigned long bend; //start and end for range of block lengths
	std::string missingchar;
	for (int i=0;i<argc;i++)
	{
		if ( string(argv[i]) == "-i" ) 
    	{
        	InFilePath = argv[i+1];
		}

		if ( string(argv[i]) == "-o" ) 
    	{
        	OutFilePath = argv[i+1];
 		}
		
		if ( string(argv[i]) == "-l" ) 
    	{
		 	LogFilePath = argv[i+1];
		}

		if ( string(argv[i]) == "-b" ) 
    	{
		 	bstart = strtoul( argv[i+1], NULL, 10);
		 	bend = strtoul( argv[i+2], NULL, 10);
		}

		if ( string(argv[i]) == "-m" ) 
    	{
		 	missingchar = argv[i+1];
		}
	}
	
	//test whether all required files specified on the command line exist
	vector<std::string> BadFiles;
	if (fileExists(InFilePath) == 0) 
	{
		std::string bf = "InFilePath = ";
		bf += InFilePath;
		BadFiles.push_back(bf);
	}
	
	if (BadFiles.size() > 0)
	{
		{
			cout << "\nThe following variables appear to contain misspecified paths:\n";
			for (unsigned int i=0;i<BadFiles.size();++i)
			{
				cout << "  " << BadFiles[i] << "\n";
			}
			cout << "\nPlease check the command line.  Quitting...\n\n";
		}
		exit (EXIT_FAILURE);
	}
	
	
	//read the input file
	vector<vector<std::string> > bufvec2d = MyReadInfile(InFilePath);
	cout << bufvec2d.size() - 2 << " individuals, " << bufvec2d[0].size() << " SNPs\n";
	
	//open log file
	ofstream logger;
	logger.open(LogFilePath);
	logger.close(); //quick open close done to clear any existing file each time program is run
	logger.open(LogFilePath, ios::out | ios::app); //open file in append mode
	logger << "b	chromosome	n loci	Mean allele count	SD allele count	n haplotype length	Mean haplotype length	SD haplotype length\n";
	
	//COMBINE ADJACENT ALLELES INTO HAPLOTYPES
	cout << "Compressing haplotypes...\n";
	vector<vector<std::string> > hapvec; //initialize combined haplotype vector
	vector<std::string> SNPsizevec; //will contain size of newly fused SNP region
	vector<std::string> chrvec; //will contain chromosomal location of newly fused region
	
	vector<std::string> chr = bufvec2d[0]; //get list of chromosome designation for each SNP
	vector<std::string> pos = bufvec2d[1]; //get list of positions of SNPs on chromosome
	unsigned long long maxpos = 0; //updated as reference in MyVecToUll
	vector<unsigned long long> posull = MyVecToULL(pos, maxpos); //convert string vector to unsigned long long vector
	vector<std::string> ind(bufvec2d.size()-2); //get list of individual id's
	for (unsigned int i=2;i<bufvec2d.size();++i) ind[i-2] = bufvec2d[i][0];
	
	//cycle through range of blocklengths, backwards, since larger b values run faster
	for (unsigned int b=bend;b>=bstart;--b)
	{
		vector<vector<std::string> >().swap(hapvec); //clear hapvec
		vector<vector<std::string> > hapvec( (bufvec2d.size() - 2) ); //size hapvec, # indiv is same as bufvec2d minus two header lines
		vector<std::string>().swap(SNPsizevec); //clear SNPsizevec
		vector<std::string>().swap(chrvec); //clear chrvec

		//cycle through SNPs, one individual at a time
		for (unsigned long i=2;i<bufvec2d.size();++i) //start at third row
		{
			cout << bufvec2d[i][0] << ", individual "<<i-1<<"/"<<bufvec2d.size() - 2<<" with blocklength "<<b<<"/"<<bend<<"\n";
		
			vector<std::string> currindiv = bufvec2d[i]; //for convenience get the current indiv SNPs as a separate vector
			hapvec[i-2].push_back(currindiv[0]); //add the indiv id to hapvec, noting correction to index i for starting on third line of bufvec2d
		
			//cycle through SNPs
			unsigned long j = 1;
			while (j<currindiv.size()) //j = allele index (first column is sample name)
			{
				string startchr = chr[j-1];//get the chromosome of the starting SNP for the fusion
				//unsigned long long startpos = strtoull( pos[j-1].c_str(), NULL, 10 ); //get start position, convert string to unsigned long long
				unsigned long long startpos = posull[j-1]; //get start position
				unsigned long k = 0;
				std::string newallele;
				while ( (k < b) && (j<currindiv.size()) ) //b=blocklength
				{
					if (chr[j-1] == startchr) //verify that current SNP is on the same chromosome as the starting SNP for the fusion
					{
						newallele += currindiv[j];
						++k;
						++j;
					}
					else break; //leave the while loop if you encounter a new chromosome, do not ++j.
								//A truncated fusion product may remain!!
				}
		
				//calculate length of SNP region fused, add to vector SNPsizevec
				//do this only for individual #1, since it is the same for all
				if (i == 2)
				{
					//unsigned long long endpos = strtoull( pos[j-2].c_str(), NULL, 10 );//determine the end position, it is the previous SNP, (j-1)-1 = j-2
					unsigned long long endpos = posull[j-2];//determine the end position, it is the previous SNP, (j-1)-1 = j-2
					unsigned long long SNPlen = endpos - startpos + 1;  //+1 to include the SNP position on both ends
					if ( SNPlen > maxpos ) 
					{
						cout << "An error has occurred. endpos="<<endpos<<", startpos="<<startpos<<", SNPlen="<<SNPlen<<"\n  SNP order may not be consecutive.\n";
					}
			
					stringstream ss;
					ss << SNPlen;
					string str = ss.str();
					SNPsizevec.push_back(str);
			
					//make note of the chromosomal location of the new fusion product
					chrvec.push_back(startchr);
				}
			
				//add the new fused allele to the data for this individual
				hapvec[i-2].push_back(newallele);
			
			}//snp
		}//individual
		
		
		
		//translate concatenated binary codes into integer allele calls
		vector<vector<int> > hapvecint = MyRecode(hapvec, missingchar);
		
		//WRITE LOG FILE
		cout << "Calculating statistics...\n";
		time_t startm, endm;
		time (&startm);

		vector<int> allelecounts(hapvecint[0].size());
		vector<unsigned long long> sizevec(hapvecint[0].size());
		vector<unsigned long long> SNPsizevecULL = MyVecToULL(SNPsizevec, maxpos); //convert string vector to unsigned long long vector


		//calculate mean allele counts and haplotype lengths for entire data set, log
		//get max value for each haplotype in hapvecintj. number of alleles = max + 1 (since 0 is an allele)
		for (unsigned int i=0;i<hapvecint[0].size();++i) //haplotypes
		{
			allelecounts[i] = MyGetAlleleCount(i, hapvecint); //get count of alleles at haplotype h	
		}
		unsigned long nac = allelecounts.size();
		double meanac = MyMean(allelecounts);
		double sdac = MyStdev(meanac, allelecounts);
		unsigned long nhl = sizevec.size();
		double meanhl;
		double sdhl;
		if (b == 1)
		{
			meanhl = 1;
			sdhl = 0;
		}
		else
		{
			meanhl = MyMeanULL(SNPsizevecULL);
			sdhl = MyStdevULL(meanhl, SNPsizevecULL);
		}
		logger << b << "\t" << "0" << "\t" << nac << "\t" << meanac << "\t" << sdac << "\t" << nhl << "\t" << meanhl << "\t" << sdhl << "\n"; //chrname = 0 means all chromosomes

		//calculate chromosome specific means, log. have to use chrvec, not chr, because chr has all SNPs, not fused haplotypes
		unsigned long long h = 0; //h indexes thru haplotypes
		while (h != hapvecint[0].size())
		{
			std::string chrname = chrvec[h]; //get the current chromosome name
			vector<int>().swap(allelecounts); //clear allelecounts
			allelecounts.resize( count(chrvec.begin(), chrvec.end(), chrname) ); //resize vector to number of haplotypes from chromosome chrname
			vector<unsigned long long>().swap(sizevec); //clear sizevec
			sizevec.resize(allelecounts.size()); //resize sizevec to same value that was calculated for allelecounts (faster than doing the count() again)
			
			unsigned int z = 0;
			while (chrvec[h] == chrname) //while current haplotype h is on chromosome chrname
			{
				allelecounts[z] = MyGetAlleleCount(h, hapvecint); //get count of alleles at haplotype h
				sizevec[z] = SNPsizevecULL[h]; //place size of haplotype h on chromosome chrname in a vector
				++h;
				if (h == chrvec.size()) break; //get out of loop after last chromosome
				++z;
			}
			
			//calculate mean allele counts and haplotype lengths for current chromosome chrname, log
			nac = allelecounts.size();
			meanac = MyMean(allelecounts);
			sdac = MyStdev(meanac, allelecounts);
			nhl = sizevec.size();
			if (b==1)
			{
				meanhl = 1;
				sdhl = 0;
			}
			else
			{
				meanhl = MyMeanULL(sizevec);
				sdhl = MyStdevULL(meanhl, sizevec);
			}
			logger << b << "\t" << chrname << "\t" << nac << "\t" << meanac << "\t" << sdac << "\t" << nhl << "\t" << meanhl << "\t" << sdhl << "\n";
		}
		
		//stop the clock
		time (&endm);
		double dif = difftime (endm,startm);
		if (dif==1) cout << "  " << dif << " second.\n";	
		else cout << "  " << dif << " seconds.\n";


		//WRITE NEW DATA SET
		//create the write path by appending blocklength
		std::string OutFilePathS = OutFilePath; //char* to std::string
		OutFilePathS += ".b";
		stringstream ss;
		ss << b;
		OutFilePathS += ss.str();
		
		//write the new data matrix to the output file
		ofstream output;
		output.open(OutFilePathS.c_str()); //convert std::string to const char* via c_str()
		output.close(); //quick open close done to clear any existing file each time program is run
		output.open(OutFilePathS.c_str(), ios::out | ios::app); //open file in append mode
	
		//write chromosomal location
		for (unsigned int i=0;i<chrvec.size();++i)
		{
			if (i == chrvec.size() - 1) output << chrvec[i];
			else output << chrvec[i] << " ";
		}
		output << "\n";
	
		//write fused SNP size
		for (unsigned int i=0;i<SNPsizevec.size();++i)
		{
			if (i == SNPsizevec.size() - 1) output << SNPsizevec[i];
			else output << SNPsizevec[i] << " ";
		}
		output << "\n";
	
		//write new SNPs
		for (unsigned long i=0;i<hapvecint.size();++i)
		{
			stringstream ss;
			ss << ind[i];
			output << ss.str() << " ";
			for (unsigned long long j=0;j<hapvecint[i].size();++j)
			{
				output << hapvecint[i][j] << " ";
			}
			output << "\n";
		}
		
		//wrap up write step
		output.close();

	}//blocklength range		
	
	//wrap up write to log
	logger.close();
	
	return 0;
}
