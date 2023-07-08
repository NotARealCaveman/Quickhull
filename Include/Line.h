#pragma once
#include "Point3.h"
#include "Transform.h"

namespace Manifest_Math
{	
	//implicit line - Plücker Coordinates
	//l{v|m}, where v is the direction of the line and m is the moment
	//assume two points on a line p1,p2 then v=p1-p2, m = p1 x p2
	//if l is normalized, meaning all 6 l/||v||, then m is the distance from the origin to the line
	//if l is not normalized, then ||m||/||v|| is the perpendicular distance to the origin
	struct MFline
	{
		MFline() = default;
		MFline(const MFfloat& vx, const MFfloat& vy, const MFfloat& vz, const MFfloat& mx, const MFfloat& my, const MFfloat& mz);
		MFline(const MFvec3& _v, const MFvec3& _m);

		MFvec3 v;
		MFvec3 m;
	};	
	//transforms a line from some coordSys to sys A
	MFline Transform(const MFline& l, const MFtransform& hA);
	//transforms a line from some coordSys to sys A taking advantage of properties of orthogonal M3
	MFline TransformOrthonormal(const MFline& l, const MFtransform& hA);
	//let L(t)=p+tv, where p = point, and tv = scalar direction of line  from points p
	//returns distance from point q to line l->v
	const MFfloat DistancePointLine(const MFpoint3& q, const MFpoint3& l, const MFvec3& v);
	//N denotes line direction is unit vector
	const MFfloat NDistancePointLine(const MFpoint3& q, const MFpoint3& l, const MFvec3& v);
	const MFfloat DistanceLineLine(const MFpoint3& p1, const MFvec3 v1, const MFpoint3& p2, const MFvec3& v2);	
}
