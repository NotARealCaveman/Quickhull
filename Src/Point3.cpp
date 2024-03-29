#include "Point3.h"

using namespace Manifest_Math;


//given the line segment AB and point Q,  point p can be found by
//P(t)=A+t(B-A) where t is the projecttion of c onto AB
//since a line segment is used and not inifite length, t must be clamped
MFpoint3 Manifest_Math::ClosestPointfromlineSegment(const MFpoint3& q, const MFpoint3& start, const MFpoint3& end)
{
	const MFvec3 ab = end - start;
	MFfloat t = Dot(q - start, ab) / Dot(ab, ab);
	t = Max(t, 0.0f);
	t = Min(t, 1.0f);

	return start + ab * t;
}

MFfloat Manifest_Math::PointDistanceFromLineSegmentSquared(const MFpoint3& q, const MFpoint3& 
start, const MFpoint3& end)
{
	const MFvec3 ab = end - start, ac = q - start, bc = q - end;
	const MFfloat proj = Dot(ac, ab);
	const MFfloat projLength = Dot(ac, ac);
	const MFfloat segmentLength = Dot(ab, ab);
	const MFfloat segmentEnd = Dot(bc, bc);
	if (proj <= 0.0f)//not on segment
		return projLength;
	else if (proj >= segmentLength)//past segment
		return segmentEnd;
	//projection lies on segment
	 return Dot(ac, ac) - proj * proj / segmentLength;
}