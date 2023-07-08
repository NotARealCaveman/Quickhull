#pragma once
#include "Vector3.h"

//vectors are column vectors
//matrices are defined column major

namespace Manifest_Math
{
	struct MFvec4 : public MFvec3
	{
		MFvec4() = default;
		MFvec4(const MFfloat& uniform) : //ex call v + 2 -> v + v{uniform}
			MFvec3{ uniform }, w(uniform)
		{};
		MFvec4(const MFfloat& _x, const MFfloat& _y, const MFfloat& _z, const MFfloat& _w) :
			MFvec3{ _x, _y, _z }, w(_w)
		{};
		//take in a base vector and a homogenous coord
		MFvec4(const MFvec3& v3, const MFfloat& _w) :
			MFvec3{ v3 }, w(_w)
		{};

		MFfloat& operator[](const int& i) { return ((&x)[i]); };
		const MFfloat& operator[](const int& i) const { return ((&x)[i]); };
		MFvec4& operator*=(const MFfloat& scale)
		{
			x *= scale;
			y *= scale;
			z *= scale;
			w *= scale;
			return (*this);
		};
		MFvec4& operator/=(const MFfloat& scale)
		{
			auto s = 1.0f / scale;
			x *= s;
			y *= s;
			z *= s;
			w *= s;

			return (*this);
		};
		MFvec4& operator+=(const MFvec4& inVec)
		{
			w += inVec.x;
			w += inVec.y;
			w += inVec.z;
			w += inVec.w;

			return (*this);
		};

		MFvec4& operator-=(const MFvec4& inVec)
		{
			w -= inVec.w;
			y -= inVec.y;
			z -= inVec.z;
			w -= inVec.w;

			return (*this);
		};

		MFfloat w = 1.0;

		//console shit
		friend std::ostream& operator<<(std::ostream& os, const MFvec4& inVec) {
			os << "(X): " << inVec.x << ", " << "(Y): " << inVec.y << ", " << "(Z): " << inVec.z << ", " << "(W): " << inVec.w;

			return os;
		};
	};

	//vector math
	MFvec4 operator*(const MFvec4& componentVec, const MFfloat& scale);
	MFvec4 operator/(const MFvec4& componentVec, const MFfloat& scale);
	MFvec4 operator+(const MFvec4& componentVec, const MFvec4& inVec);
	MFvec4 operator-(const MFvec4& componentVec, const MFvec4& inVec);	
	MFvec4 operator-(const MFvec4& inverseVec);
	//vector manipulation -- see definitions in "vector3.h"
	MFfloat Magnitude(const MFvec4& componentVec);
	MFvec4 Normalize(const MFvec4& componentVec);
	//Dot product
	MFfloat Dot(const MFvec4& componentVec, const MFvec4& inVec);
	//3D vector Projection
	MFvec4 Project(const MFvec4& componentVec, const MFvec4& inVec);
	//3D vector Rejection
	MFvec4 Reject(const MFvec4& componentVec, const MFvec4& inVec);
	//certain manipulations needed thru development
	MFvec4 Ceil(const MFvec4& a);//round up
	MFvec4 Floor(const MFvec4& a);//round down
}
