#include "Vector2.h"

using namespace Manifest_Math;

MFvec2::MFvec2(const MFfloat& u)
	:x(u), y(u) {};

MFvec2::MFvec2(const MFfloat& _x, const MFfloat& _y)
	: x(_x), y(_y) {};

MFvec2 Manifest_Math::Normalize(const MFvec2& inVec)
{
	return inVec / std::sqrtf(inVec.x * inVec.x + inVec.y * inVec.y);
}

MFfloat Manifest_Math::Dot(const MFvec2& a, const MFvec2& b)
{
	return  a.x * b.x + a.y * b.y;
}