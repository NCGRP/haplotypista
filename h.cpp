#include "h.hpp"


/*
To compile:  use "make"
Usage: see README.txt
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

//determine the size of a file on disk
std::ifstream::pos_type filesize(const char* filename)
{
	std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
	return in.tellg(); 
}

vector<vector<std::string> > MyReadInfile(char* InFilePath)
{

	//start the clock
	cout << "Reading data...\n";
	time_t startm,endm;
	time (&startm);

	char * buffer = MyBigRead(InFilePath);
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
	vector<vector<int> > hapvecint( hapvec.size(), vector<int>( (hapvec[0].size())-1) ); //declare and size vector to hold new integer coded alleles
	
	for (unsigned int i=1;i<hapvec[0].size();++i) //go thru each locus, start at 1 to skip row label
	{
		//cout << "loc=" << i << "\n";
		
		vector<std::string> AllelesEncountered; //will contain the unique set of alleles at the locus
		for (unsigned int k=0;k<hapvec.size();++k) //go thru all individuals
		{
			int ColIndex = i;
			std::string a = hapvec[k][ColIndex];
			
			//search for missing substring
			std::size_t found = a.find(missingchar);
			if (found!=std::string::npos) hapvecint[k][ColIndex-1] = -9; //add the missing data value #this puts a missing value for a haplotype that contains any missing nucleotide. this might not be the smartest way to deal with this. too conservative? --PR 1/27/14
			else
			{
				int AlleleInt; //the new, integerized, name of the allele
				std::vector<std::string>::iterator itr = std::find(AllelesEncountered.begin(), AllelesEncountered.end(), a);
				if (itr != AllelesEncountered.end()) //the allele has been found before
				{
			//cout << "old allele a="<<a<<"\n"; 

					AlleleInt = itr - AllelesEncountered.begin(); //convert itr to index, the index is the integerized allele name
					hapvecint[k][ColIndex-1] = AlleleInt; //add the integerized allele name
				}
				else // you have a new allele
				{
			//cout << "new allele a="<< a << " " << "AllelesEncountered.size()=" << AllelesEncountered.size() << "\n"; 
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

vector<Location> MyProcessLocations(std::string in, std::string& DoReduce)
{
	//change status of switch
	DoReduce = "yes";

	//initialize
	std::string a;
	vector<Location> locs;
	std::replace( in.begin(), in.end(), ',', ' '); //replace all commas with space
	vector<string> ranges = split(in); //split space delimited string into vector with entry like 1.2355:1.6952 to specify a single range
	
	for (unsigned long long i=0;i<ranges.size();++i)
	{
		a = ranges[i];
		std::replace( a.begin(), a.end(), ':', ' ' );
		std::replace( a.begin(), a.end(), '.', ' ' );
		vector<string> l = split(a); //split space delimited string into a vector with position 0 and 2 specify chromosome, positions 1 and 3 specify beginloc, endloc
		
		//convert vector<string> to vector<unsigned long long>, add values to Location struct, push_back onto return variable
		vector<unsigned long long> ll;
		for (unsigned int j=0;j<l.size();++j)
		{
			ll.push_back( stoull(l[j]) );	
		}
		if ( ll[0] != ll[2] ) 
		{
			cout << "Range specification (-g) can not span chromosomes . Quitting...\n";
			exit (EXIT_FAILURE);
		}
		Location loc { std::make_pair(ll[0],i), ll[1], ll[3] };
//		Location loc { ll[0], ll[1], ll[3] };
		locs.push_back(loc); 
	}
	return locs;
}

vector<std::string> MyDoReduce(vector<vector<std::string> >& bufvec2d, vector<Location> locs)
{
	std::string chr;
	std::string rangeid;
	unsigned long long beginloc;
	unsigned long long endloc;
	vector<unsigned long long> incvec; //will hold indexes of all columns to be included from rows containing genotypic data
	vector<std::string> rangevec; //will hold rangeid from MyProcessLocations for each column, original fragment name in
	                              //bufvec2d[0] will be converted to this rangeid to treat each user defined fragment (-g)
	                              //independently so that multiple ranges per chromosome and overlapping ranges are allowed
	vector<std::string> chrwrite; //will hold the original fragment names, not the supplanting range ids, so that final output
	                              //can contain sensible locations
		
	for (unsigned long long l=0;l<locs.size();++l)
	{
		
		chr = to_string(locs[l].chr.first); //get the chromosome for current location
		rangeid = to_string(locs[l].chr.second); //get the range id for the current range as string
		beginloc = locs[l].beginloc;
		endloc = locs[l].endloc;
		//cout << "l=" << l << " chr=" << chr << " rangeid=" << rangeid << " beginloc=" << beginloc << " endloc=" << endloc << "\n";
		
		//extract index values to locate columns in bufvec2d that correspond to regions of interest
		//if statement tests: chr is correct, position is between range specified in locs[l]
		for (unsigned long long i=0;i<bufvec2d[0].size();++i)
		{
			if ( (bufvec2d[0][i] == chr) && (stoull(bufvec2d[1][i]) >= beginloc) && (stoull(bufvec2d[1][i]) <= endloc) )
			{
				incvec.push_back(i+1);
				rangevec.push_back(rangeid); //list of rangeids to supplant original fragment names
				chrwrite.push_back(chr); //original fragment name preserved in chrwrite
				//cout << "i=" << i << " rangeid=" << rangeid << "\n";
				//cout << "i+1=" << i+1 << "\n";	
				//cout << " " << stoull(bufvec2d[0][i]) << "\n";
				//cout << " " << stoull(bufvec2d[1][i]) << "\n";
			}
		}
	}
	
	//extract relevant data
	std::string pb;
	vector<vector<std::string> > tempvec2d(bufvec2d.size()); //initialize temporary storage with same row number as bufvec2d
	//add sample names to rows 4 thru end
	for (unsigned int j=3;j<bufvec2d.size();++j)
	{
		pb = bufvec2d[j][0];
		tempvec2d[j].push_back(pb);
	}
	//iterate through rows, extracting columns from each row corresponding to index values
	for (unsigned int j=0;j<bufvec2d.size();++j)
	{
		//iterate over list of indexed
		for (unsigned int i=0;i<incvec.size();++i)
		{
			//treat 1st three rows differently, use indexed position -1
			if ( j == 0 )
			{
				pb = rangevec[i]; //replace original fragment name with range id
				tempvec2d[j].push_back(pb);
			}
			else if ( j == 1 || j == 2 )
			{
				pb = bufvec2d[j][incvec[i]-1];
				tempvec2d[j].push_back(pb);
			}
			else
			{
				pb = bufvec2d[j][incvec[i]];
				tempvec2d[j].push_back(pb);
			}
		}
	}

	//overwrite bufvec2d for update by reference
	vector<vector<std::string> >().swap(bufvec2d); //clear bufvec2d prior to filling with new information
	bufvec2d = tempvec2d;

	return chrwrite;
}

vector<std::string> MyProcessPopFile(char* PopFilePath, int ploidy)
{	
	//read file into a vector, one line per item
	std::ifstream popfile(PopFilePath);
	vector<std::string> popvec;
	std::string foo;
	while (getline(popfile, foo))//get line from s, put in foo, consecutively
	{
		popvec.push_back(foo);  
	}

	//place into set to determine unique pop names
	std::set<std::string> popuniq;
	for (unsigned int i=0;i<popvec.size();++i)
	{
		popuniq.insert(popvec[i]); 
	}

	//label each pop name with consecutive values for each ind
	vector<std::string> rowlabs(popvec.size());
	std::string ins = "no"; //match indicator
	for(set<string>::const_iterator it = popuniq.begin(); it != popuniq.end(); it++)
	{
		unsigned int x = 1;
		unsigned int i = 0;
		
		//cout << "i=" << i << ", x=" << x << "\n";

		while ( i < (popvec.size()) )
		{
			//test whether current vector element is the same as the unique pop id from set popuniq
			if ( popvec[i] == *it )
			{
				ins = "yes";
				//repeat ploidy times if above true
				for (int p=0;p<ploidy;++p)
				{
					//add row label for p-th consecutive row
					rowlabs[i] = (popvec[i] + " " + to_string(x));
					++i;
		//cout << "i=" << i << ", x=" << x << "\n";
				}
			}
			
			if ( ins == "no") ++i;  //index popvec if no match to set iterator
			if ( ins == "yes" ) 
			{
				ins = "no"; //reset match indicator
				++x; //index ind label only after labeling all rows corresponding to the same ind
			}

		}
	}
	
	return rowlabs;
}	

int MyWriteNormal(std::string OutFilePathS, vector<std::string> chrvec, vector<std::string> chrwritevec, vector<std::string> SNPsizevec, vector<std::string> midptvec, vector<std::string> aacatvec, vector<std::string> blockstartvec, vector<std::string> blockendvec, vector<std::string> ind, vector<vector<int> > hapvecint, std::string DoReduce)
{
	cout << "Writing standard output...\n";
	OutFilePathS += ".txt";
	//write the new data matrix to the output file
	ofstream output;
	output.open(OutFilePathS.c_str()); //convert std::string to const char* via c_str()
	output.close(); //quick open close done to clear any existing file each time program is run
	output.open(OutFilePathS.c_str(), ios::out | ios::app); //open file in append mode

	//write chromosomal location
	vector<std::string> cvec;
	if ( DoReduce == "yes" ) cvec = chrwritevec;
	else cvec = chrvec;
	for (unsigned int i=0;i<cvec.size();++i)
	{
		if (i == cvec.size() - 1) output << cvec[i];
		else output << cvec[i] << " ";
	}
	output << "\n";

	//write fused SNP size
	for (unsigned int i=0;i<SNPsizevec.size();++i)
	{
		if (i == SNPsizevec.size() - 1) output << SNPsizevec[i];
		else output << SNPsizevec[i] << " ";
	}
	output << "\n";

	//write haplotype block midpoint
	for (unsigned int i=0;i<midptvec.size();++i)
	{
		if (i == midptvec.size() - 1) output << midptvec[i];
		else output << midptvec[i] << " ";
	}
	output << "\n";

	//write amino acid category codes
	for (unsigned int i=0;i<aacatvec.size();++i)
	{
		if (i == aacatvec.size() -1) output << aacatvec[i];
		else output << aacatvec[i] << " ";
	}
	output << "\n";
	
	//write haplotype block startpoint
	for (unsigned int i=0;i<blockstartvec.size();++i)
	{
		if (i == blockstartvec.size() - 1) output << blockstartvec[i];
		else output << blockstartvec[i] << " ";
	}
	output << "\n";
	
	//write haplotype block endpoint
	for (unsigned int i=0;i<blockendvec.size();++i)
	{
		if (i == blockendvec.size() - 1) output << blockendvec[i];
		else output << blockendvec[i] << " ";
	}
	output << "\n";
	
	//write alleles
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
	
	return 0;
}

int MyWriteM(std::string OutFilePathS, vector<std::string> chrvec, vector<std::string> chrwritevec, vector<std::string> midptvec, vector<std::string> ind, vector<vector<int> > hapvecint, vector<std::string> rowlabs, std::string DoReduce)
{
	cout << "Writing M+ output...\n";

	//modify file name for m+ .dat and .var format
	std::string datfile = OutFilePathS + ".dat";
	std::string varfile = OutFilePathS + ".var";
	
	//write the new data matrix to dat file
	ofstream output;
	output.open(datfile.c_str()); //convert std::string to const char* via c_str()
	output.close(); //quick open close done to clear any existing file each time program is run
	output.open(datfile.c_str(), ios::out | ios::app); //open file in append mode

	std::string linei;
	std::string datai;
	std::string foo;
	for (unsigned int i=0;i<hapvecint.size();++i)
	{
		datai.clear();
		for (unsigned int j=0;j<hapvecint[i].size();++j)
		{
			//change coding of missing data value -9
			if (hapvecint[i][j] == -9) foo = "9999";
			else foo = to_string(hapvecint[i][j]);
			
			//add string value to output variable
			if (j == hapvecint[i].size() - 1) datai += foo;
			else datai += (foo + " ");
		}
		//cout << "datai=" << datai << "\n";
		
		linei = rowlabs[i] + " " + ind[i] + " " + datai;
		//cout << "linei=" << linei << "\n";
		output << linei << "\n";
	}

	//wrap up write dat file
	output << '\0'; //add null terminator, for some this reason this is necessary for output that exceeds some size unknown to me
	output.close();
	

	//create locus names
	std::string loc;
	vector<std::string> locnames;
	vector<std::string> cvec; //holds the original fragment names to write
	if ( DoReduce == "yes" ) cvec = chrwritevec;
	else cvec = chrvec;
	unsigned long long j = 1;
	for (unsigned int i=0;i<cvec.size();++i)
	{
		loc = cvec[i] + "." + midptvec[i] + "." + to_string(j);
		locnames.push_back(loc);
		++j;
	}

	//write the var file
	output.open(varfile.c_str()); //convert std::string to const char* via c_str()
	output.close(); //quick open close done to clear any existing file each time program is run
	output.open(varfile.c_str(), ios::out | ios::app); //open file in append mode
	
	output << "code 0\nindividu 0\nSample 1 0 0 1 5\n"; //write var file header
	for (unsigned int i=0;i<locnames.size();++i)
	{
		output << locnames[i] << " 2 1 0 1 5\n";
	}

	//wrap up write var file
	output << '\0'; //add null terminator, probably not necessary here because var file will seldom exceed size at which null terminator is no longer added, but having two at the end can't hurt
	output.close();

	return 0;
}


/***************MAIN*****************/

int main( int argc, char* argv[] )
{
	//parse the command line for options
	char* InFilePath;
	char* OutFilePath;
	char* PopFilePath;
	std::string LogFilePath;
	unsigned long bstart;
	unsigned long bend; //start and end for range of block lengths
	std::string missingchar;
	int ploidy;
	vector <Location> locs; //holds ranges for genomic regions of interest
	std::string DoReduce = "no";
	std::string MakeM = "no";
	for (int i=0;i<argc;i++)
	{
		if ( string(argv[i]) == "-i" ) 
    	{
        	InFilePath = argv[i+1];
		}

		if ( string(argv[i]) == "-o" ) 
    	{
        	OutFilePath = argv[i+1];
        	LogFilePath = std::string(OutFilePath) + "log.txt";
 		}
		
/*		if ( string(argv[i]) == "-l" ) 
    	{
		 	LogFilePath = argv[i+1];
		}
*/
		if ( string(argv[i]) == "-b" ) 
    	{
		 	bstart = strtoul( argv[i+1], NULL, 10);
		 	bend = strtoul( argv[i+2], NULL, 10);
		}

		if ( string(argv[i]) == "-m" ) 
    	{
		 	missingchar = argv[i+1];
		}

		if ( string(argv[i]) == "-p" ) 
    	{
		 	ploidy = atoi(argv[i+1]);
		}

		if ( string(argv[i]) == "-g" ) 
    	{
		 	//process command line supplied genes of interest
		 	locs = MyProcessLocations(argv[i+1], DoReduce);
		}
		
		if ( string(argv[i]) == "-v" ) 
    	{
			MakeM = "yes";
			PopFilePath = argv[i+1];
		}

	}
	
	//process genes of interest that are piped in
	//test whether there is a piped list. does program see a terminal or a file/pipe as the source of the stdin?
	if (isatty(fileno(stdin))); //puts("stdin is connected to a terminal"); //do nothing
  	else
  	{
    	//puts("stdin is NOT connected to a terminal");
    	//stdin is from a pipe, |, or file, <
		std::string pipedin;
		std::string s;
		while (std::getline(std::cin, s))
		{
			pipedin += s;
			pipedin += ",";
		}
		pipedin.pop_back();//remove hanging comma
		locs = MyProcessLocations(pipedin, DoReduce);
		cout << "\n" << "pipedin=" << pipedin << "\n";
	
	}
	
	/* write out command line options and genes of interest specification
		for (int i=0;i<argc;i++)
		{
			cout << argv[i] << "\n";	
		}
	
		for (int i=0;i<locs.size();++i)
		{
			cout << i << " " << locs[i].chr.first << " " << locs[i].chr.second << " " << locs[i].beginloc << " " << locs[i].endloc << "\n";
		}
	*/
	
	//test whether all required files specified on the command line exist
	vector<std::string> BadFiles;
	std::string bf;
	if (fileExists(InFilePath) == 0) 
	{
		bf = "InFilePath = ";
		bf += InFilePath;
		BadFiles.push_back(bf);
	}
	if (MakeM == "yes")
		{
		if (fileExists(PopFilePath) == 0) 
		{
			bf = "PopFilePath = ";
			bf += PopFilePath;
			BadFiles.push_back(bf);
		}
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

	//test that header lines are equal length
	if ( ( bufvec2d[0].size() != bufvec2d[1].size() ) || ( bufvec2d[0].size() != bufvec2d[2].size() ) || ( bufvec2d[1].size() != bufvec2d[2].size() ) )
	{
		cout << "\nThe first three lines of the input file do not all contain the same number of elements.  Quitting...\n\n";
		exit (EXIT_FAILURE);
	}
	
	//print out stats on input file
	cout << "Data set contains:\n";
	cout << "  " << (bufvec2d.size() - 3) / ploidy << " individuals, " << bufvec2d.size() - 3 << " haplotypes, " << bufvec2d[0].size() << " SNPs, ploidy = " << ploidy << "N\n";
	
	//read the population specification file, process it
	vector<std::string> rowlabs;
	if (MakeM == "yes") rowlabs = MyProcessPopFile(PopFilePath, ploidy);
	
	//reduce data set to regions of interest
	//when DoReduce="yes", chr and chrvec contain range ids for fragments, not original fragment names.
	//when DoReduce="yes", chrwrite and chrwritevec contain the original fragment names that correspond to the -g selected regions
	vector<std::string> chrwrite;
	vector<std::string> chrwritevec;
	if (DoReduce == "yes")
	{
		chrwrite = MyDoReduce(bufvec2d, locs); //updates bufvec2d as reference
		
		//print out stats on reduced file
		cout << "Reducing data set to genomic positions defined by -g...\n";
		cout << "Reduced data set contains:\n";
		cout << "  " << (bufvec2d.size() - 3) / ploidy << " individuals, " << bufvec2d.size() - 3 << " haplotypes, " << bufvec2d[0].size() << " SNPs, ploidy = " << ploidy << "N\n";
	}
	
	/*print out the data set	
		for (unsigned int i=0;i<bufvec2d.size();++i)
		{
			for (unsigned int j=0;j<bufvec2d[i].size();++j)
			{
				cout << bufvec2d[i][j] << " ";
			}
			cout << "\n";
		}
	*/ 
	
	//open log file
	ofstream logger;
	logger.open(LogFilePath.c_str());
	logger.close(); //quick open close done to clear any existing file each time program is run
	logger.open(LogFilePath.c_str(), ios::out | ios::app); //open file in append mode
	logger << "b	chromosome	n loci	Mean allele count	SD allele count	n haplotype length	Mean haplotype length	SD haplotype length\n";
	
	//COMBINE ADJACENT ALLELES INTO HAPLOTYPES
	cout << "Compressing haplotypes...\n";
	vector<vector<std::string> > hapvec; //initialize combined haplotype vector
	vector<std::string> SNPsizevec; //will contain size of newly fused SNP region
	vector<std::string> chrvec; //will contain chromosomal location of newly fused region
	vector<std::string> midptvec; //will contain the nucleotide position of the midpoint of the newly fused region
	vector<std::string> aacatvec; //will contain the counts of non-genic, synonymous, and non-synonymous sites contained in the fused region in the form {4 1 0, 3 0 2, 5 0 0}
	vector<std::string> blockstartvec; //will contain the nucleotide position of the first SNP of the haplotype block
	vector<std::string> blockendvec; //will contain the nucleotide position of the last SNP of the haplotype block

	vector<std::string> chr = bufvec2d[0]; //get list of chromosome designation for each SNP, this is rangeid when DoReduce="yes"
	vector<std::string> pos = bufvec2d[1]; //get list of positions of SNPs on chromosome
	vector<std::string> aacat = bufvec2d[2]; //get list of amino acid category of each SNP
	
	unsigned long long maxpos = 0; //updated as reference in MyVecToUll
	vector<unsigned long long> posull = MyVecToULL(pos, maxpos); //convert string vector to unsigned long long vector
	vector<std::string> ind(bufvec2d.size()-3); //get list of individual id's, not used until write step
	for (unsigned int i=3;i<bufvec2d.size();++i) ind[i-3] = bufvec2d[i][0];
	
	//cycle through range of blocklengths, backwards, since larger b values run faster
	for (unsigned int b=bend;b>=bstart;--b)
	{
		vector<vector<std::string> >().swap(hapvec); //clear hapvec
		vector<vector<std::string> > hapvec( (bufvec2d.size() - 3) ); //size hapvec, # indiv is same as bufvec2d.size() minus three header lines
		vector<std::string>().swap(SNPsizevec); //clear SNPsizevec
		vector<std::string>().swap(midptvec); //clear midptvec
		vector<std::string>().swap(chrvec); //clear chrvec
		vector<std::string>().swap(chrwritevec); //clear chrwritevec
		vector<std::string>().swap(aacatvec); //clear aacatvec
		vector<std::string>().swap(blockstartvec); //clear blockstartvec
		vector<std::string>().swap(blockendvec); //clear blockendvec
		
		//cycle through SNPs, one individual at a time
		for (unsigned long i=3;i<bufvec2d.size();++i) //start at fourth row
		{
			cout << bufvec2d[i][0] << ", haplotype "<<i-2<<"/"<<bufvec2d.size() - 3<<" with blocklength "<<b<<"/"<<bend<<"\n";
		
			vector<std::string> currindiv = bufvec2d[i]; //for convenience get the current indiv SNPs as a separate vector
			hapvec[i-3].push_back(currindiv[0]); //add the indiv id to hapvec, noting correction to index i for starting on fourth line of bufvec2d
		
			//cycle through SNPs
			unsigned long j = 1; //j = SNP position index (starts at 1 because first column is sample name)
			unsigned long m = 0; //m = amino acid code index (starts at 0)
			while (j<currindiv.size()) 
			{
				string startchr = chr[j-1];//get the chromosome of the starting SNP for the fusion, =rangeid when DoReduce="yes"
				string startchrwrite; //initialize variable only useful MyReduce="yes"
				if (DoReduce == "yes") startchrwrite = chrwrite[j-1]; //get the original fragment name of the starting SNP for the fusion
				unsigned long long startpos = posull[j-1]; //get start position
				unsigned long k = 0;
				std::string newallele;
				std::string newaacats;
				while ( (k < b) && (j<currindiv.size()) ) //b=blocklength
				{
					if (chr[j-1] == startchr) //verify that current SNP is on the same chromosome as the starting SNP for the fusion
					{
						newallele += currindiv[j]; //concatenate the allele calls for the current block
						if (i ==3)
						{
							newaacats += aacat[m]; //concatenate the amino acid category codes for the current block, only need to do this when working on first individual since it will be the same for all
						}
						++k;
						++j;
						++m;
					}
					else break; //leave the while loop if you encounter a new chromosome, do not ++j.
								//A truncated fusion product may remain!!
				}
		
				
				
					//calculate length of SNP region fused, the midpoint of haplotype block, and the frequency of non-genic, non-synonymous, and synonymous SNPs in the block
					//add to vectors
					//do this only for individual #1, since it is the same for all
					if (i == 3)
					{
						//calculate length
						unsigned long long endpos = posull[j-2];//determine the end position, it is the previous SNP, (j-1)-1 = j-2
						unsigned long long SNPlen = endpos - startpos + 1;  //+1 to include the SNP position on both ends
						if ( SNPlen > maxpos ) 
						{
							cout << "An error has occurred. endpos="<<endpos<<", startpos="<<startpos<<", SNPlen="<<SNPlen<<"\n  SNP order may not be consecutive.\n";
						}
					
						//calculate midpt
						unsigned long long midpt = startpos+(SNPlen/2);
					
						//calculate amino acid category code frequencies
						unsigned long long ng = std::count(newaacats.begin(), newaacats.end(), '0'); //number of non-genic SNPs contained in the haplotype block
						unsigned long long syn = std::count(newaacats.begin(), newaacats.end(), '1'); //number of synonymous SNPs contained in the haplotype block
						unsigned long long ns = std::count(newaacats.begin(), newaacats.end(), '2'); //number of non-synonymous SNPs contained in the haplotype block
					
						//test whether the number of SNPs included in the haplotype is equal to the blocklength, i.e. is the haplotype truncated?
						//if the haplotype is not truncated, include it in the output
						if (newallele.length() == b)
						{
							//cout << startchr << "\t" << j << "\t" << startpos << "\t" << endpos << "\t" << SNPlen << "\t" << midpt << "\n";
			
							//add the block length to SNPsizevec
							stringstream ss;
							ss << SNPlen;
							string str = ss.str();
							SNPsizevec.push_back(str);
					
							//add the midpoint of the haplotype block to midptvec
							stringstream sd;
							sd << midpt;
							str = sd.str();
							midptvec.push_back(str);
					
							//add the amino acid category codes to aacatvec
							stringstream se;
							se << ng << ":" << syn << ":" << ns;
							str = se.str();
							aacatvec.push_back(str);
							
							//add the position of the first SNP in block to blockstartvec
							stringstream sf;
							sf << startpos;
							str = sf.str();
							blockstartvec.push_back(str);
							
							//add the position of the last SNP in block to blockendvec
							stringstream sg;
							sg << endpos;
							str = sg.str();
							blockendvec.push_back(str);
							
							//make note of the chromosomal location of the new fusion product
							chrvec.push_back(startchr);
							if (DoReduce == "yes") chrwritevec.push_back(startchrwrite);
						}
					}
			
					//add the new fused allele to the data for this individual
					//but only if the fusion product is not truncated
					//cout << "b=" << b << "newaacats.length()=" << newaacats.length() << "\t" << "newallele.length()=" << newallele.length() << "\n";
					if (newallele.length() == b)
					{
						hapvec[i-3].push_back(newallele);
					}
			
				

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


		//WRITE NEW DATA SETS
		//create the write path by appending blocklength
		std::string OutFilePathS = OutFilePath; //char* to std::string
		OutFilePathS += ".b";
		stringstream ss;
		ss << b;
		OutFilePathS += ss.str();
		
		MyWriteNormal(OutFilePathS, chrvec, chrwritevec, SNPsizevec, midptvec, aacatvec, blockstartvec, blockendvec, ind, hapvecint, DoReduce);
		if (MakeM == "yes") MyWriteM(OutFilePathS, chrvec, chrwritevec, midptvec, ind, hapvecint, rowlabs, DoReduce);
		
	}//blocklength range		
	
	//wrap up write to log
	logger.close();
	
	return 0;
}
