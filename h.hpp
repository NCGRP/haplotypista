#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <math.h>
#include <numeric>
#include <set>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <sstream>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>
using namespace std;


/***************STRUCTS*****************/

//holds ranges for genomic regions of interest
struct Location
{
  std::pair <unsigned long long, unsigned long long> chr; // without -g...fragment name from input file : fragment name from input file
                                                          // with -g .....fragment name from input file : unique name for -g defined range
//  unsigned long long chr; //chromosome or fragment name
  unsigned long long beginloc;
  unsigned long long endloc;
};

/***************CLASSES*****************/


/***************VARIABLES*****************/


/***************SHARED FUNCTIONS WITHIN FILES*****************/


/***************SHARED FUNCTIONS BETWEEN FILES*****************/
