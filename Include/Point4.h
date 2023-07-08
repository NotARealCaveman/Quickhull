#pragma once
#include "Vector4.h"
#include "Point3.h"

//used when defining a point in 3d space
//unlike Point3D w is not implicitly 1, thus P(p|w)
//P=pw
//pNormal = P/w
//division by w moves vector Pxyzw to plane where w=1
namespace Manifest_Math
{
	struct MFpoint4 : public MFvec4
	{
		MFpoint4() = default;
		MFpoint4(const MFfloat& _x, const MFfloat& _y, const MFfloat& _z, const MFfloat& _w) : MFvec4{ _x,_y,_z,_w } {};
		MFpoint4(const MFvec3& p, const MFfloat& w) : MFvec4{ p,w } {};
		MFpoint4(const MFfloat& u) : MFvec4{ {u},1.0f } {};
		MFpoint4(const MFvec4& v) : MFvec4{ v } {};		
	};
	inline MFpoint4 operator+(const MFpoint4& componentPoint, const MFpoint4& inPoint)
	{
		auto&& interpolate = MFpoint4{ static_cast<MFvec3>(componentPoint) + static_cast<MFvec3>(inPoint),componentPoint.w + inPoint.w };
		
		return interpolate / interpolate.w;
	}
	inline MFpoint4 operator-(const MFpoint4& componentPoint, const MFpoint4& inPoint)
	{
		auto&& interpolate = MFpoint4{ static_cast<MFvec3>(componentPoint) - static_cast<MFvec3>(inPoint),componentPoint.w - inPoint.w };

		//division by 0 is fine for floats
		//it is defined as the point at +/-inf for x >0,<0 respectively.
		//for the division 0/0 it is defined as -nan(ind)
		//Not-A-Number(indeterminate) -> result of invalid opperation 
		return interpolate / interpolate.w;		
	}
	//normalized a P4 to P3|w=1
	MFpoint3 Normalize(const MFpoint4& componentPoint);
}
