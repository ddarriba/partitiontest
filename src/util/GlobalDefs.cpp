#include "GlobalDefs.h"
#include "Utilities.h"

using namespace std;

namespace partest {

#ifdef PTHREADS
	int number_of_threads = 1;
#endif

	bitMask protModelsMask = Utilities::binaryPow(8*sizeof(bitMask)-1)-1;
	bool ckpAvailable = false;
	string ckpPath;
	string ckpStartingTree = "starting_tree";
}
