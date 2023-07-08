#include "Simplex.h"

using namespace Manifest_Math;

MFbool Manifest_Math::PositiveHalfSpace(const MFvec3& direction, const MFvec3& AO)
{
	//DLOG(36, "D*AO: " << Dot(direction, AO));
	return Dot(direction, AO) > 0;
}