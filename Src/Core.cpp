#include "Core.h"

using namespace Manifest_Math;


//DO NOT USE
// DO NOT KEEP SCROLLING
// 
// 
// 
// NOT STAY AWAY
// STOP
// 
// STOP
// I TOLD YOU
// 
// STOP!!!
// 
//!!WARNING!! FUNCTION IS SLOWER(and less accurate...) THAN 1.0F/STD::SQRT(X) included for why the fuck not
const MFfloat Manifest_Math::Isqrt(MFfloat&& x2)//const removed for ability to reinterpret cast without needing temp varbiable
{//l variable is needed to be able ot take address for return cast
	auto l = 0x5f3759df - (*reinterpret_cast<long*>(&x2) >> 1);//x guess	
	return *reinterpret_cast<MFfloat*>(&l) * (1.5f - ((x2 * 0.5) * *reinterpret_cast<MFfloat*>(&l) * *reinterpret_cast<MFfloat*>(&l)));//newtonian
}