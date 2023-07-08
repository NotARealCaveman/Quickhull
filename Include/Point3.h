#pragma once
#include "Vector3.h"

//used when defining a point in 3d space
//4th component is assumed to be 1, not required to be stored in memory
namespace Manifest_Math
{
	struct MFpoint3 : public MFvec3
	{
		MFpoint3() = default;
		MFpoint3(const MFfloat& x, const MFfloat& y, const MFfloat& z) :MFvec3{ x,y,z } {};
		MFpoint3(const MFfloat& u) : MFvec3{ u } {};
		MFpoint3(const MFvec3& v) : MFvec3{ v } {};
	};		
	//returns the point on a line segment closet to an arbitary non collinear point
	MFpoint3 ClosestPointfromlineSegment(const MFpoint3& q, const MFpoint3& start, const MFpoint3& end);
	MFfloat PointDistanceFromLineSegmentSquared(const MFpoint3& q, const MFpoint3& start, const MFpoint3&);
}
