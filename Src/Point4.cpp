#include "Point4.h"

using namespace Manifest_Math;

MFpoint3 Manifest_Math::Normalize(const MFpoint4& componentPoint)
{
	return	{ static_cast<MFvec3>(componentPoint) / componentPoint.w };
}