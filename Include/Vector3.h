#pragma once
#include <iostream>
#include <algorithm>

#include "DebugLogger.h"

#include "Core.h"

//vectors are column vectors
//matrices are defined column major
namespace Manifest_Math
{
	struct MFvec3
	{
		MFvec3() = default;
		MFvec3(const MFfloat uniform) : //ex call v + 2 -> v + v{uniform}
			x{ uniform }, y{ uniform }, z{ uniform }
		{};
		MFvec3(const MFfloat& _x, const MFfloat& _y, const MFfloat& _z) :
			x{ _x }, y{ _y }, z{ _z }
		{};

		MFfloat& operator[](const MFint32& i) { return ((&x)[i]); };
		const MFfloat& operator[](const MFint32& i) const { return ((&x)[i]); };
		MFvec3& operator*=(const MFfloat& scale)
		{
			x *= scale;
			y *= scale;
			z *= scale;
			return (*this);
		};		
		MFvec3& operator/=(const MFfloat& scale)
		{
			auto s = 1.0f / scale;
			x *= s;
			y *= s;
			z *= s;
			return (*this);
		};		
		MFvec3& operator+=(const MFvec3& inVec)
		{
			x += inVec.x;
			y += inVec.y;
			z += inVec.z;
			return (*this);
		};

		MFvec3& operator-=(const MFvec3& inVec)
		{
			x -= inVec.x;
			y -= inVec.y;
			z -= inVec.z;
			return (*this);
		};
		//used for map in isosurface
		MFbool operator<(const MFvec3& other) const
		{
			if (x != other.x)
				return x < other.x;
			if (y != other.y)
				return y < other.y;
			return z < other.z;
		};

		MFfloat* begin() 
		{
			return &x;
		}
		MFfloat* end() 
		{
			return (&z)+1;
		}


		const MFfloat* begin() const
		{
			return &x;
		}
		const MFfloat* end() const
		{
			return (&z) + 1;
		}

		MFfloat x;
		MFfloat y;
		MFfloat z;

		//console shit
		friend std::ostream& operator<<(std::ostream& os, const MFvec3& inVec) {
			os << /*&inVec.x <<*/ "(X): " <<  inVec.x << ", " << /*&inVec.y << */"(Y): " << inVec.y << ", " << /*&inVec.z <<*/ "(Z): " << inVec.z;

			return os;
		};
	};

	void swap(MFvec3& a, MFvec3& b);

	inline MFbool operator==(const MFvec3& a, const MFvec3& b)
	{
		if ((a.x == b.x) && (a.y == b.y) && (a.z == b.z))
			return true;
		else 
			return false;
	}
	/*inline MFbool operator!=(const MFvec3& a, const MFvec3& b)
	{
		return !(a == b);
	}*/
	//vector manipulation
	inline MFvec3 operator*(const MFvec3& componentVec, const MFfloat& scale)
	{
		return MFvec3{ componentVec.x * scale,componentVec.y * scale,componentVec.z * scale };
	};
	inline MFvec3 operator/(const MFvec3& componentVec, const MFfloat& scale)
	{
		auto s = 1.0f / scale;
		return MFvec3{ componentVec.x * s,componentVec.y * s,componentVec.z * s };
	};
	inline MFvec3 operator+(const MFvec3& componentVec, const MFvec3& inVec)
	{
		return MFvec3{ componentVec.x + inVec.x, componentVec.y + inVec.y, componentVec.z + inVec.z };
	};
	inline MFvec3 operator-(const MFvec3& componentVec, const MFvec3& inVec)
	{
		return MFvec3{ componentVec.x - inVec.x, componentVec.y - inVec.y, componentVec.z - inVec.z };
	};
	inline MFvec3 operator-(const MFvec3& inverseVec)
	{
		return MFvec3{ -inverseVec.x, -inverseVec.y, -inverseVec.z };
	};	

	inline MFvec3 Lerp(const MFvec3& v0, const MFvec3& v1, const MFfloat& t)
	{
		return v0 + (v1 - v0) * t;
	}
	//vector math	
	//length of vector
	MFfloat Magnitude(const MFvec3& componentVec);
	//squared length - saves the sqrt and allows for faster length comparisons against 2 sqrd vectors
	MFfloat MagnitudeSquared(const MFvec3& componentVec);
	//expands to /c\ =c/||c||
	//internally expands to c/sqrt(cX^2+cY^2+cZ^2)
	//which expands to 1/||c||cX,1/||c||cX,1/||c||cY,1/||c||cZ
	MFvec3 Normalize(const MFvec3& componentVec);
	//Cross Product(3D wedge product)
	//returns vector orthogonal to plane whice Va&vB lie
	MFvec3 Cross(const MFvec3& componentVec, const MFvec3& inVec);
	//Dot product, of two unit vectors, will be in range [-1,1]
	//x->-1, vectors point >90* apart
	//x=0 vectors are orthogonal
	//x->1, vectors point <90* together
	//Componentwise dot* product
	//if vA==vB then Dot(V,V) = V^2
	MFfloat Dot(const MFvec3& componentVec, const MFvec3& inVec);
	//returns theta from dot product
	//strips theta from cos(0) with cos-1
	//theta is returned as radians
	MFfloat DotTheta(const MFvec3& componentVec, const MFvec3& inVec);
	//trigonometric definition 
	//a*b = ||a||||b||cos(0), let 0 = theta
	//a and b are normalized
	//then a*b = cos(0) which 0=cos^-1(A*B)

	//3D vector Projection
	//a||b=(a*b)/b^2*b
	//gives a vector projA which is parallel to b
	//dot product returns "likeness" of two vectors, normalized to b if not already ( /v\=v/||v||, ||v||=v*v)
	//projected dot product then scales vB
	MFvec3 Project(const MFvec3& componentVec, const MFvec3& inVec);
	//3D vector Rejection
	//a_|_b = a - a||b
	//gives a vector rejA which is perpendicular to b
	//the combination of Projections and Rejctions can be use for orthogonalization
	//provides vectors which are orthogonal to a desired vector and can be used as the sides of a right triangle
	//the hypotenuse of this triangle is the vector which the projection/rejection is done on
	MFvec3 Reject(const MFvec3& componentVec, const MFvec3& inVec);
	//Spherical Linear Interpolation
	//rotates a vector n by some axis-angle x,0
	MFvec3 Slerp(const MFvec3& n, const MFvec3& a, const MFfloat& t);
	//returns the roll pitch and yaw of an axis angle
	MFvec3 AxisAngleToEuler(const MFvec3& a, const MFfloat& t);
	//manipulation needed thru development
	MFvec3 Max(const MFvec3& a, const MFvec3& b);//returns maxs of both
	MFvec3 Min(const MFvec3& a, const MFvec3& b);//returns mins of both
	MFvec3 Ceil(const MFvec3& a);//round down
	MFvec3 Floor(const MFvec3& a);//round up
	MFvec3 Square(const MFvec3& a);//self dot product
	MFvec3 ComponentMultiply(const MFvec3& a, const MFvec3& b);
	MFvec3 VAbs(const MFvec3& a);//absolute vector values
}