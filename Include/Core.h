#pragma once
#include <iostream>
#include <iomanip>

#include "Typenames.h"

using namespace Manifest_Utility;

namespace Manifest_Math
{

	constexpr MFfloat Pi = 3.1415927f;
	constexpr MFfloat Pi2 = Pi * 2.0f;
	constexpr MFfloat EpsilonMax = 1.0f - FLT_EPSILON;	

	inline const MFfloat Radians(const MFfloat& degrees) { return (Pi / 180.0f) * degrees; };
	inline const MFfloat Degrees(const MFfloat& radians) { return (180.0f / Pi) * radians; };
//calcluates Cotanget as 1.0f/tangent(t)
	inline const MFfloat Cot(const MFfloat& t) { return 1.0f / tanf(t); };
	inline const MFfloat& Min(const MFfloat& a, const MFfloat& b) { return a < b ? a : b; };
//returns a if greater than clamping value
	inline const MFfloat ClampMin(const MFfloat& a, const MFfloat& clamp) { return a < clamp ? clamp : a; };	
	inline const MFfloat& Max(const MFfloat& a, const MFfloat& b) { return a > b ? a : b; };
//returns a if less than clamping vlaue
	inline const MFfloat ClampMax(const MFfloat& a, const MFfloat& clamp) { return a > clamp ? clamp : a; };
//clamps value between desired range
	inline const MFfloat ClampRange(const MFfloat& lower, const MFfloat& a, const MFfloat& upper)
	{
		return Min(Max(lower, a),upper);
	}
//once i figure it out
	inline const MFfloat Abs(const MFfloat& a) { return  fabsf(a); };	
	inline const MFfloat Sqrt(const MFfloat& a) { return sqrtf(a); };
	inline const MFfloat Sin(const MFfloat& a) { return sinf(a); };
	inline const MFfloat Cos(const MFfloat& a) { return cosf(a); };
}
























namespace Manifest_Math
{
	const MFfloat Isqrt(MFfloat&& x2);
}

