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

/*double MyStdev(double m, vector<int> v)
{
	double accum = 0.0;
	std::for_each ( v.begin(), v.end(), [&](const double d) {
		accum += (d - m) * (d - m);
	});
	double stdev = sqrt(accum / (v.size()-1));
	return stdev;
}
*/

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

/*double MyStdevULL(double m, vector<unsigned long long> v)
{
	double accum = 0.0;
	std::for_each 
	(	v.begin(), 
		v.end(),
		[&](const double d) {accum += (d - m) * (d - m);}
	);
	double stdev = sqrt(accum / (v.size()-1));
	return stdev;
}
*/
/*
vector<std::string> MySetKernel(char* KerFilePath)
{
	//declare variables
	std::string foo;
	vector<std::string> KernelAccessionList;
	
	//add mandatory accessions to list
	std::ifstream infile;
	infile.open(KerFilePath);
        while( !infile.eof() ) // To get all the lines.
       {
	        std::getline (infile, foo); // Saves the line in foo.
	        if (foo == "") break; //get out of while when foo string is empty
	        KernelAccessionList.push_back(foo);
	    }

	return KernelAccessionList;
}


//removes duplicate values from a list but does not alter sort order
vector<std::string> unsortedRemoveDuplicates(vector<std::string> numbers)
{
	vector<std::string> uniqvec;
	std::string b;
	for (unsigned int i=0;i<numbers.size();++i)
	{
		b = numbers[i];
		if (std::find( uniqvec.begin(), uniqvec.end(), b) == uniqvec.end() ) uniqvec.push_back(b); //number has not been seen, add it
	}
	return uniqvec;
}

//tabulates the ploidy for each locus
vector<int> GetPloidy(vector<std::string> AllLociNameList, vector<std::string> UniqLociNameList)
{
	int n;
	std::string b;
	vector<int> PloidyList;
	for (unsigned int i=0;i<UniqLociNameList.size();++i)
	{
		b = UniqLociNameList[i];
		n = std::count (AllLociNameList.begin(), AllLociNameList.end(), b);
		PloidyList.push_back(n);
	}
	
	return PloidyList;
}

int MyProcessVarFile(char* VarFilePath, vector<int>& AllColumnIDList, vector<std::string>& AllLociNameList, vector<int>& ActiveColumnIDList, vector<std::string>& ActiveLociNameList, vector<int>& TargetColumnIDList, vector<std::string>& TargetLociNameList, vector<vector<int> >& ColKeyToAllAlleleByPopList, vector<int>& ReferenceOrTargetKey, vector<int>& PloidyList, vector<std::string>& UniqLociNameList)
{
    //declare variables
    std::string foo;
    vector<std::string> foovector;
    unsigned int k;
    
    int i=0; // i is the row number
	std::ifstream infile;
	infile.open(VarFilePath);
    while( !infile.eof() ) // To get all the lines.
    {
	    std::getline (infile, foo); // Saves the line in foo.
	        
	    //split foo on whitespace
		foovector = split(foo);

		//identify active columns with qualitative data, classify those as reference or target
		if (foovector[1] == "2") //column holds qualitative data
		{
			AllColumnIDList.push_back(i);
			AllLociNameList.push_back(foovector[0]);
				
			if ((foovector[2] == "1") && (foovector[3] == "0")) //reference variable
			{
				ActiveColumnIDList.push_back(i);
				ActiveLociNameList.push_back(foovector[0]);
			}
			else if ((foovector[2] == "0") && (foovector[3] == "1"))  //target variable
			{
				TargetColumnIDList.push_back(i);
				TargetLociNameList.push_back(foovector[0]);
			}
		}
	     
		i++;
		foovector.clear();  //zero vector foovector
    }
	
	infile.close();
	
	
	//make a key showing which loci are target and which are reference
	//this will be used later to sort out the AllAlleleByPopList
	std::string uniqloc;
	std::string currloc;
	
	//find unique locus names, retains sort order
	UniqLociNameList = unsortedRemoveDuplicates(AllLociNameList);

	//tabulate the ploidy for each locus
	PloidyList = GetPloidy(AllLociNameList, UniqLociNameList);
	
	//define 2d vector that is key to columns, size to number of unique loci
	ColKeyToAllAlleleByPopList.resize( UniqLociNameList.size() ); //size to number of loci
	int b;
	for (unsigned int i=0;i<UniqLociNameList.size();++i)
	{
		uniqloc = UniqLociNameList[i];
		for (k=0;k<AllLociNameList.size();++k)
		{
			currloc = AllLociNameList[k];
			if (currloc == uniqloc)
			{
				b = AllColumnIDList[k]; //get the corresponding value in the columnID list
				ColKeyToAllAlleleByPopList[i].push_back(b);	
			}	
		}		
	}
	
	//define a 1d vector that describes whether the loci in the ColKey are reference(0) or target(1)
	double rt;
	ReferenceOrTargetKey.resize( ColKeyToAllAlleleByPopList.size() ); //sized to same length as key
	for (unsigned int i=0;i<ColKeyToAllAlleleByPopList.size();++i)
	{
		rt=0;
		//test whether all elements are categorized as reference or as target, if not, raise error
		for (k=0;k<ColKeyToAllAlleleByPopList[i].size();++k)
		{
			b=ColKeyToAllAlleleByPopList[i][k];
			if(std::find(ActiveColumnIDList.begin(), ActiveColumnIDList.end(), b) != ActiveColumnIDList.end())
			{
   			 	//column is reference
   			 	rt=rt+0;
			} 
			else if(std::find(TargetColumnIDList.begin(), TargetColumnIDList.end(), b) != TargetColumnIDList.end())
			{
    			//column is target
    			rt=rt+1;
			}
		}
		
		//test whether columns in key represent reference or target loci
		if (rt == 0) ReferenceOrTargetKey[i] = 0; //it is a reference
		else if (  rt/ColKeyToAllAlleleByPopList[i].size() == 1 ) ReferenceOrTargetKey[i] = 1; //it is a target
		else 
		{
			cout << "ERROR:  Some loci are described as both reference and target.  Please check your var file. Quitting...\n\n";
			exit (EXIT_FAILURE);
		}
	}
	return 0;
}
*/

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

/*
//adds an allele to the appropriate locus in the ActiveAlleleList
void MyUpdateActiveAlleleList(vector<vector<std::string> >& ActiveAlleleList, int CurrItemIndex, std::string NewAllele) //updates ActiveAlleleList by reference
{	
	vector<std::string> OldAlleles;
	
	//add to the current array being built or to the appropriate existing array
	OldAlleles = ActiveAlleleList[CurrItemIndex]; //pull out the list of alleles currently held for this locus
	OldAlleles.push_back (NewAllele); //add the new allele to the existing list of alleles for the current locus
	ActiveAlleleList[CurrItemIndex] = OldAlleles; //update ActiveAlleleList by reference
}

//reduces the master vector of all alleles into subsets containing reference or target loci only
int MyReduceToRef(vector<vector<vector<int> > > AllAlleleByPopList, vector<int> ReferenceOrTargetKey, vector<vector<vector<int> > >& ActiveAlleleByPopList, vector<vector<vector<int> > >& TargetAlleleByPopList)
{
	unsigned int r, t, i, j;
	vector<int> b;
	
	//resize level 1 of vectors, they contain the same number of populations as AllAlleleByPopList
	ActiveAlleleByPopList.resize(AllAlleleByPopList.size());
	TargetAlleleByPopList.resize(AllAlleleByPopList.size());
	
	//reserve level 2 of vectors, they contain the number of loci in ReferenceOrTargetKey
	r = std::count(ReferenceOrTargetKey.begin(), ReferenceOrTargetKey.end(), 0);
	t = std::count(ReferenceOrTargetKey.begin(), ReferenceOrTargetKey.end(), 1);
	for (i=0;i<AllAlleleByPopList.size();++i)
	{
		ActiveAlleleByPopList[i].reserve(r);
		TargetAlleleByPopList[i].reserve(t);
	}

	
	for (i=0;i<AllAlleleByPopList.size();++i) //iterate thru populations
	{
		for (j=0;j<AllAlleleByPopList[i].size();++j) //iterate thru loci
		{
			b.clear();
			b = AllAlleleByPopList[i][j];
			
			if (ReferenceOrTargetKey[j] == 0) //it is a reference locus
			{
				ActiveAlleleByPopList[i].push_back(b);
			}
			
			else if (ReferenceOrTargetKey[j] == 1) //it is a target locus
			{
				TargetAlleleByPopList[i].push_back(b);
			}
		}
	}
	//ActiveAlleleByPopList & TargetAlleleByPopList have now been updated
	return 0;
}
*/

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
	
	/*
	//stop the clock
	time (&endm);
	dif = difftime (endm,startm);
	if (procid == 0) 
	{
		if (dif==1) cout << "    " << dif << " second.\n";	
		else cout << "    " << dif << " seconds.\n";	
	}

	
	if (procid == 0) cout << "  Building data structures...\n";
	time (&startm);
	
	//simplify bufvec2d by retaining only the first column, with POPID, then swap clear to save memory
	vector<string> PopIDvec;
	for (i=0;i<bufvec2d.size();++i)
		PopIDvec.push_back(bufvec2d[i][0]);
	vector<vector<std::string> >().swap(bufvec2d); //clear bufvec2d

	//break up vector into a 3d vector by population:  { { {pop1ind1elems},{pop1ind2elems},...}, { {pop2ind1elems},{pop2ind2elems},...} }
	vector<vector<vector<int> > > ByPop3d;
	for (i=0;i<bufvec2dint.size();++i)
	{
		//NewPopID = bufvec2d[i][0];  //access string in bufvec2d
		NewPopID = PopIDvec[i];
		IndivPerPop.push_back(NewPopID); //add the pop ID to a list to calc pop sizes later

		if (NewPopID != OldPopID) //then create a new population in ByPop3D
		{
			ByPop3d.resize(ByPop3d.size() + 1);
			
			//add the new population name to the AccessionNameList
			FullAccessionNameList.push_back(NewPopID);
		}
		
		//push row of integer elements on current line, onto last item of ByPop3d, which might be a new population as added just above
		//remember that the first three columns are 0 in bufvec2dint, so will also be 0 in ByPop3d
		ByPop3d[ByPop3d.size()-1].push_back(bufvec2dint[i]);
		
		OldPopID = NewPopID;
	}
	vector<vector<int> >().swap(bufvec2dint); //clear bufvec2dint, the integerized data is in ByPop3d

	//stop the clock
	time (&endm);
	dif = difftime (endm,startm);
	if (procid == 0) 
	{
		if (dif==1) cout << "    " << dif << " second.\n";	
		else cout << "    " << dif << " seconds.\n";	
	}

	//print out ByPop3d
	for (i=0;i<ByPop3d.size();++i)
	{
		cout << "Pop" << i << "\n";
		for (j=0;j<ByPop3d[i].size();++j)
		{
			cout << " Ind" << j << "  ";
			for (k=0;k<ByPop3d[i][j].size();++k)
			{
				cout << ByPop3d[i][j][k] << ",";
			}
			cout << "\n";
		}
	
	}

	if (procid == 0) cout << "  Condensing data...\n";
	time (&startm);

	//resize AllAlleleByPopListSet
	AllAlleleByPopListSet.resize(ByPop3d.size());//resize number of populations
	for (i=0;i<AllAlleleByPopListSet.size();++i)
	{
		AllAlleleByPopListSet[i].resize(ColKeyToAllAlleleByPopList.size()); //resize number of loci
																		   //the index in ColKey is the locus index in AllAllelesByPopList level 2
																		   //the value of ColKey is the index of the allele in ByPop3d level 3
	}
		

	//calculate size of AllAlleles
	AllAlleles.reserve(row*AllColumnIDList.size());
	
	//condense alleles by locus in AllAllelesByPopList, within each population
	int AlleleIndex;
	for (i=0;i<ByPop3d.size();++i) //go thru pops
	{
		for (j=0;j<ByPop3d[i].size();++j) //go thru indivs
		{
			for (k=0;k<ColKeyToAllAlleleByPopList.size();++k) //go through each locus
			{
				for (l=0;l<ColKeyToAllAlleleByPopList[k].size();++l) //assign columns to loci
				{
					AlleleIndex = ColKeyToAllAlleleByPopList[k][l];
					int NewAllele = ByPop3d[i][j][AlleleIndex]; //get the allele in the specified column
					AllAlleles.push_back(NewAllele); //add the allele to the list of all alleles, missing data included
					if (NewAllele != -9999) //exclude missing data
						AllAlleleByPopListSet[i][k].insert(NewAllele); //add the allele to the set of unique alleles at locus k, pop i	
				}
			}
		}
	}
	
	vector<vector<vector<int> > >().swap(ByPop3d); //clear ByPop3d
			
	//stop the clock
	time (&endm);
	dif = difftime (endm,startm);
	if (procid == 0) 
	{
		if (dif==1) cout << "    " << dif << " second.\n";	
		else cout << "    " << dif << " seconds.\n";	
	}


	return 0;
}
*/

/*//removes duplicate alleles and missing data (9999) from the supplied vector
vector<std::string> MyFilterDuplicates(vector<std::string> ListToFilter)
{
	//remove duplicates
	sort( ListToFilter.begin(), ListToFilter.end() );
	ListToFilter.erase( std::unique( ListToFilter.begin(), ListToFilter.end() ), ListToFilter.end() );
	
	//remove missing data
	ListToFilter.erase( std::remove( ListToFilter.begin(), ListToFilter.end(), "9999" ), ListToFilter.end() );
	
	return ListToFilter;
}

//removes duplicate numbers from the supplied vector
vector<std::string> MyFilterDuplicatesII(vector<std::string> ListToFilter)
{
	//remove duplicates
	sort( ListToFilter.begin(), ListToFilter.end() );
	ListToFilter.erase( std::unique( ListToFilter.begin(), ListToFilter.end() ), ListToFilter.end() );
	
	return ListToFilter;
}
*/

/*
//returns maximum number of alleles possible at each locus for active and target
vector<int> MyGetMaxs(vector<vector<vector<int> > > ActiveAlleleByPopList)
{
	unsigned int i, j, k;
	vector<int> ActiveMaxAllelesList;
	vector<int> CurrLoc;
	set<int> NewSet;
	
	for (i=0;i<ActiveAlleleByPopList[0].size();++i)
	{
		NewSet.clear();
		for (j=0;j<ActiveAlleleByPopList.size();j++)
		{
			CurrLoc = ActiveAlleleByPopList[j][i]; //you are traversing the locus 'column' of the 3d grid
												   //get locus i for population j
			for (k=0;k<CurrLoc.size();++k)
			{
				NewSet.insert(CurrLoc[k]);	//place alleles into set to eliminate redundancies
			}
		}
		ActiveMaxAllelesList.push_back(NewSet.size()); //the set size after adding all populations for locus i is the maximum number of alleles
	}
	return ActiveMaxAllelesList; 
}
*/

bool fileExists(const char *fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

/*
vector<int> MyRemoveTargetAlleles(vector<int> AllAlleles, vector<int> AllColumnIDList, vector<int> TargetColumnIDList)
{
	unsigned int checker, j;
	unsigned int i=0;
	vector<int> AllRefAlleles;
	
	while (i<AllAlleles.size())
	{
		j=0;
		while (j<AllColumnIDList.size()) //this routine follows along, mapping column ids to target or, implicitly, reference
		{
			checker = AllColumnIDList[j];		
			if(std::find(TargetColumnIDList.begin(), TargetColumnIDList.end(), checker) != TargetColumnIDList.end()) 
			{
				// TargetColumnIDList contains checker
				++j;
				++i;
			} 
			else 
			{
				//TargetColumnIDList does not contain checker
				AllRefAlleles.push_back(AllAlleles[i]);
				++j;
				++i;
			}
		}
	}
	return AllRefAlleles;
}	

//places a continuous string of alleles into a 2d vector with samples as rows
void MyMakeRefAllelesIntoRows(vector<int> AllRefAlleles, vector<std::string> ActiveLociNameList, vector<vector<int> >& RefAllelesIntoRows)
{
	unsigned int j; //counts items per row
	unsigned int i=0, k=0; //k is the row number
	while (i<AllRefAlleles.size())
	{
		j=0;
		while (j<ActiveLociNameList.size()) //keeps track of the number of alleles traversed before setting the next row
		{
			RefAllelesIntoRows[k][j] = AllRefAlleles[i];
			++i;
			++j;
		}
		++k;
	}
}	

//remove columns of data to make 2d vector of alleles by locus 
void MyMakeRefAllelesByLocus(vector<vector<int> > RefAllelesIntoRows, vector<string> ActiveLociNameList, vector<std::pair<std::string, vector<int> > >& RefAllelesByLocus)
{
	vector<std::pair<std::string, vector<int> > > lola;
	std::pair<std::string, vector<int> > la; //locus name, allele list pair
	
	unsigned int i, j, k;
	int locindex;
	std::string locname, b;
	vector<std::string> foo;
	for (i=0;i<ActiveLociNameList.size();++i)
	{
		locname = ActiveLociNameList[i];
		//test whether the list of locus names/alleles already contains this locus
		locindex = -1; //default value will cause an error if not explicitly set
		for (j=0;j<lola.size();++j)
		{
			b = lola[j].first;
			if (locname == b) //there is already data for that locus
			{
				locindex = j;
			}
		}
		if (locindex == -1)  //locname not found, add a pair
		{
			locindex = lola.size(); 
			lola.push_back(la);//add an empty pair onto the vector of loci
			lola[locindex].first = locname;
		}
		
		//add the column of data defined by locname to the appropriate pair.second, defined by locindex in lola
		for (k=0;k<RefAllelesIntoRows.size();++k)
		{
			int al = RefAllelesIntoRows[k][i];
			if (al != -9999) lola[locindex].second.push_back(al); //ignore missing data
		}
	}
	
	//update RefAllelesByLocus
	RefAllelesByLocus = lola;
	
	vector<std::pair<std::string, vector<int> > >().swap(lola); //clear lola

	
	//print out each locus name followed by all the alleles found within it
	for (i=0;i<lola.size();++i)
	{
		cout << lola[i].first << "\n";
		foo = lola[i].second;
		for (j=0;j<foo.size();++j)
		{
			cout << " " << foo[j];
		}
		cout << "\n";
	}
	
}

//calculates allele frequencies for all alleles at all loci, updates vector of struct Alfreq, which contains the relational data
void MyCalculateAlleleFrequencies(vector<std::pair<std::string, vector<int> > > RefAllelesByLocus, vector<Alfreq>& AlleleFrequencies)
{
	unsigned int i, j, z;
	double freq;
	std::string b;
	vector<int> AllAlleles;
	vector<int> UniqAlleles;
	set<int> AlleleSet;

	Alfreq laf; //locusname, allelenames, frequencies
	static const struct Alfreq emptylaf; //this will be used to zero struct between loops
	
	for (i=0;i<RefAllelesByLocus.size();++i)
	{
		//empty containers
		laf = emptylaf;
		vector<int>().swap(UniqAlleles); //clear UniqAlleles
		AlleleSet.clear(); //clear AlleleSet
		
		//get locus name, add to struct
		b = RefAllelesByLocus[i].first;
		laf.locusname = b;
		
		//compress list of all alleles at this locus into unique alleles
		AllAlleles = RefAllelesByLocus[i].second; //get the vector of all alleles at locus i
		for (j=0;j<AllAlleles.size();++j) 
		{
			AlleleSet.insert( AllAlleles[j] ); //filter out redundant alleles by dumping the vector into a set
		}
		UniqAlleles.assign(AlleleSet.begin(), AlleleSet.end()); //assign the unique alleles to a vector
		
		//for each allele in UniqAlleles count the number of occurrences in AllAlleles, calc frequency, add allele name and freq to struct
		for (j=0;j<UniqAlleles.size();++j)
		{
			int al = UniqAlleles[j];
			laf.allelenames.push_back(al); //add allele name to struct
			z = count(AllAlleles.begin(), AllAlleles.end(), al);
			freq = double(z)/double(AllAlleles.size());
			laf.frequencies.push_back(freq); //add allele frequency to struct
		}
		
		AlleleFrequencies.push_back(laf);
	}
}
*/

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
	logger << "b	chromosome	n allele count	Mean allele count	SD allele count	n haplotype length	Mean haplotype length	SD haplotype length\n";
	
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
			cout << "Individual "<<i-1<<"/"<<bufvec2d.size() - 2<<" with blocklength "<<b<<"/"<<bend<<"\n";
		
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
