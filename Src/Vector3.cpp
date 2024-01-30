#include "Vector3.h"

using namespace Manifest_Math;

void Manifest_Math::swap(MFvec3& a, MFvec3& b)
{
	const MFvec3 t{ b };
	b = a;
	a = t;
}

//vector math	
MFfloat Manifest_Math::Magnitude(const MFvec3& componentVec)
{
	return std::sqrtf(componentVec.x * componentVec.x + componentVec.y * componentVec.y + componentVec.z * componentVec.z);
};
MFfloat Manifest_Math::MagnitudeSquared(const MFvec3& componentVec)
{
	return componentVec.x * componentVec.x + componentVec.y * componentVec.y + componentVec.z * componentVec.z;
};

//rot
MFvec3 Manifest_Math::Normalize(const MFvec3& componentVec)
{
	return componentVec / std::sqrt(componentVec.x * componentVec.x + componentVec.y * componentVec.y + componentVec.z * componentVec.z);
};

MFvec3 Manifest_Math::Cross(const MFvec3& componentVec, const MFvec3& inVec)
{
	return MFvec3
	{
		componentVec.y * inVec.z - componentVec.z * inVec.y ,
		componentVec.z * inVec.x - componentVec.x * inVec.z ,
		componentVec.x * inVec.y - componentVec.y * inVec.x
	};
};

MFfloat Manifest_Math::Dot(const MFvec3& componentVec, const MFvec3& inVec)
{
	return componentVec.x * inVec.x + componentVec.y * inVec.y + componentVec.z * inVec.z;
};

MFfloat Manifest_Math::DotTheta(const MFvec3& componentVec, const MFvec3& inVec)
{
	return std::acosf(ClampRange(-1.0f,(componentVec.x * inVec.x + componentVec.y * inVec.y + componentVec.z * inVec.z),1.0f));
};

MFvec3 Manifest_Math::Project(const MFvec3& componentVec, const MFvec3& inVec)
{
	return inVec * (Dot(componentVec, inVec) / Dot(inVec, inVec));
};

MFvec3 Manifest_Math::Reject(const MFvec3& componentVec, const MFvec3& inVec)
{
	return componentVec - inVec * (Dot(componentVec, inVec) / Dot(inVec, inVec));
}

MFvec3 Manifest_Math::Slerp(const MFvec3& v, const MFvec3& a, const MFfloat& t)
{
	return v * cosf(t) + (a * Dot(v, a) * (1 - cosf(t)) + (Cross(a, v) * sinf(t)));
}

//converts an axis-angle to euler angles
//euler angles are expresed in radians
//roll/bank-rotation about z axis
//pitch/attitude-rotation about x axis
//yaw/heading-rotation about y axis
MFvec3 Manifest_Math::AxisAngleToEuler(const MFvec3& axis, const MFfloat& theta)
{
	MFvec3 euler;//pitch-yaw-roll

	const MFfloat s = sinf(theta);
	const MFfloat c = cosf(theta);
	const MFfloat t = 1.0f - c;
	//axis is assumed to be normalized, if not normalized, normalized before passing in
	if ((axis.x * axis.y * t + axis.z * s) > 0.998f) { // north pole singularity detected
		euler.y = 2 * atan2f(axis.x * sinf(theta / 2), cosf(theta / 2));
		euler.x = Pi / 2;
		euler.z = 0;

		return euler;
	}
	if ((axis.x * axis.y * t + axis.z * s) < -0.998f) { // south pole singularity detected
		euler.y = -2.0f * atan2f(axis.x * sinf(theta / 2.0f), cosf(theta / 2.0f));
		euler.x = -Pi / 2.0f;
		euler.z = 0.0f;

		return euler;
	}

	//no singularity dected
	euler.y = atan2f(axis.y * s - axis.x * axis.z * t, 1.0f - (axis.y * axis.y + axis.z * axis.z) * t);
	euler.x = asinf(axis.x * axis.y * t + axis.z * s);
	euler.z = atan2f(axis.x * s - axis.y * axis.z * t, 1.0f - (axis.x * axis.x + axis.z * axis.z) * t);

	return euler;
}


MFvec3 Manifest_Math::Max(const MFvec3& a, const MFvec3& b)
{
	return MFvec3
	{
		std::max(a.x,b.x),
		std::max(a.y,b.y),
		std::max(a.z,b.z)
	};
}
MFvec3 Manifest_Math::Min(const MFvec3& a, const MFvec3& b)
{
	return MFvec3
	{
		std::min(a.x,b.x),
		std::min(a.y,b.y),
		std::min(a.z,b.z)
	};
}
MFvec3 Manifest_Math::Ceil(const MFvec3& a)
{
	return MFvec3
	{
		std::ceilf(a.x),
		std::ceilf(a.y),
		std::ceilf(a.z)
	};
}
MFvec3 Manifest_Math::Floor(const MFvec3& a)
{
	return MFvec3
	{
		std::floorf(a.x),
		std::floorf(a.y),
		std::floorf(a.z)
	};
}
//not to be confused with the applications of the dot product
//the Square function returns the squared components of a vector, i.e selfwise dot product
MFvec3 Manifest_Math::Square(const MFvec3& a)
{
	return MFvec3{ a.x * a.x, a.y * a.y,a.z * a.z };
}
MFvec3 Manifest_Math::ComponentMultiply(const MFvec3& a, const MFvec3& b)
{
	return MFvec3{ a.x * b.x, a.y * b.y,a.z * b.z };
}
//returns the absolute value of each axis in vec
MFvec3 Manifest_Math::VAbs(const MFvec3& a)
{
	return MFvec3{ std::fabsf(a.x), std::fabsf(a.y), std::fabsf(a.z) };
}
