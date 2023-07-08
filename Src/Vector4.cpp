#include "Vector4.h"
using namespace Manifest_Math;

//vector math
MFvec4 Manifest_Math::operator*(const MFvec4& componentVec, const MFfloat& scale)
{
	return MFvec4{ componentVec.x * scale,componentVec.y * scale,componentVec.z * scale, componentVec.w * scale };
};
MFvec4 Manifest_Math::operator/(const MFvec4& componentVec, const MFfloat& scale)
{
	auto s = 1.0f / scale;
	return MFvec4{ componentVec.x * s,componentVec.y * s,componentVec.z * s, componentVec.w * s };
};
MFvec4 Manifest_Math::operator+(const MFvec4& componentVec, const MFvec4& inVec)
{
	return MFvec4{ componentVec.x + inVec.x, componentVec.y + inVec.y, componentVec.z + inVec.z, componentVec.w + inVec.w };
};
MFvec4 Manifest_Math::operator-(const MFvec4& componentVec, const MFvec4& inVec)
{
	return MFvec4{ componentVec.x - inVec.x, componentVec.y - inVec.y, componentVec.z - inVec.z, componentVec.w - componentVec.w };
};
MFvec4 Manifest_Math::operator-(const MFvec4& inverseVec)
{
	return MFvec4{ -inverseVec.x, -inverseVec.y, -inverseVec.z , -inverseVec.w };
};
//vector manipulation -- see definitions in "vector3.h"
MFfloat Manifest_Math::Magnitude(const MFvec4& componentVec)
{
	return std::sqrt(componentVec.x * componentVec.x + componentVec.y * componentVec.y + componentVec.z * componentVec.z + componentVec.w + componentVec.w);
}
MFvec4 Manifest_Math::Normalize(const MFvec4& componentVec) { return componentVec / Magnitude(componentVec); };

//Dot product
MFfloat Manifest_Math::Dot(const MFvec4& componentVec, const MFvec4& inVec)
{
	return componentVec.x * inVec.x + componentVec.y * inVec.y + componentVec.z * inVec.z + componentVec.w * inVec.w;
}
//3D vector Projection
MFvec4 Manifest_Math::Project(const MFvec4& componentVec, const MFvec4& inVec)
{
	return inVec * (Dot(componentVec, inVec) / Dot(inVec, inVec));
}
//3D vector Rejection
MFvec4 Manifest_Math::Reject(const MFvec4& componentVec, const MFvec4& inVec)
{
	return componentVec - inVec * (Dot(componentVec, inVec) / Dot(inVec, inVec));
}

MFvec4 Manifest_Math::Ceil(const MFvec4& a)
{
	return MFvec4
	{
		std::ceilf(a.x),
		std::ceilf(a.y),
		std::ceilf(a.z),
		std::ceilf(a.w)
	};
}
MFvec4 Manifest_Math::Floor(const MFvec4& a)
{
	return MFvec4
	{
		std::floorf(a.x),
		std::floorf(a.y),
		std::floorf(a.z),
		std::floorf(a.w)
	};
}