#pragma once
#include "Line.h"
#include "Transform.h"
#include "Point4.h"

namespace Manifest_Math
{
	//a plane can be defined as
	//f[n|d], such that n represents a normal vector to some point p on the plane
	//and w represents the point at which the normal values intersect the plane f
	struct MFplane
	{
		MFplane() = default;
		//takes in normal vector of plane
		MFplane(const MFvec3& n, const MFfloat d) :
			x(n.x), y(n.y), z(n.z), w(d)
		{};
		//takes in normal components of plane
		MFplane(const MFfloat nx, const MFfloat ny, const MFfloat nz, const MFfloat d) :
			x(nx), y(ny), z(nz), w(d)
		{};
		//returns the normal vector of the plane
		inline const MFvec3& Normal()const { return reinterpret_cast<const MFvec3&>(x); };
		//sets the plane given a normal and offset from the origin
		void Set(const MFvec3& normal, const MFfloat& d)
		{
			*reinterpret_cast<MFvec3*>(&x) = normal;
			w = d;
		}
		void Set(const MFfloat& _x, const MFfloat& _y, const MFfloat& _z, const MFfloat& _d)
		{
			x = _x;
			y = _y;
			z = _z;
			w = _d;
		}
		explicit operator MFvec4() const
		{
			return reinterpret_cast<const MFvec4&>(*this);
		}
		MFfloat& operator[](const MFint32& i) { return ((&x)[i]); };
		const MFfloat& operator[](const MFint32& i) const { return ((&x)[i]); };

		//normal vector of plane
		MFfloat x;
		MFfloat y;
		MFfloat z;
		//distance from origin, f*O = O
		MFfloat w;//also known as d

		//console shit
		friend std::ostream& operator<<(std::ostream& os, const MFplane& inPlane)
		{
			os << /*&inPlane <<*/ " " << inPlane.x << "X " << inPlane.y << "Y " << inPlane.z << "Z " <<inPlane.w <<"w";
			return os;
		};
	};	

	//returns the normalized plane f->fn
	MFplane Normalize(const MFplane& f);
	//returns the transformation of fA -> fH
	const MFplane operator*(const MFplane& f, const MFtransform& h);
	//returns dot product between the plane and direction, implicit w = 0		
	//if plane is normalized then f*p=||f->p||
	const MFfloat Dot(const MFplane& f, const MFvec3& v);
	//returns dot product between the plane and point, implicit w = 1	
	const MFfloat Dot(const MFplane& f, const MFpoint3& p);
	//returns a transform matrix which will reflect a point pacross a plane f	
	const MFtransform ReflectPlane(const MFplane& f);
	//tests for intersection of line and plane
	const MFfloat InteresectionLinePlane(const MFpoint3& p, const MFvec3& v, const MFplane& f, MFpoint3& q);
	//tests for intersection of ray and plane
	const MFfloat InteresectionRayPlane(const MFpoint3& p, const MFvec3& v, const MFplane& f, MFpoint3& q);
	//tests for intersection of ray and front of the plane
	const MFfloat InteresectionRayPlaneFront(const MFpoint3& p, const MFvec3& v, const MFplane& f, MFpoint3& q);
	//tests for intersection point of three planes
	const MFfloat IntersectionThreePlanes(const MFplane& fA, const MFplane& fB, const MFplane& fC, MFpoint3& p);
	//tests for intersection point of two palnes and f[v|0]
	const MFfloat IntersectionTwoPlanes(const MFplane& fA, const MFplane& fB, MFpoint3& p, MFvec3& v);	
	//returns point4 representing the point closest to the origin on the plane farthest from the origin along the normal n		
	MFpoint4 PlanePointOrigin(const MFplane& f);
	//returns the projected point P of point Q onto the plane F
	MFpoint3 ClosestPointfromPlane(const MFpoint3& q, const MFplane& f);
}


