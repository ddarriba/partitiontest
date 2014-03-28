#include "GlobalDefs.h"
#include "Utilities.h"

using namespace std;

namespace partest {

	bitMask protModelsMask = Utilities::binaryPow(8*sizeof(bitMask))-1;;
	bool ckpAvailable = false;
	string ckpPath;
}
