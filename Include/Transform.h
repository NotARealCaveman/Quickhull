#pragma once
#include "Matrix4.h"
#include "Point3.h"

namespace Manifest_Math
{
	//inherits from MFmat4 to be used in places where expected
	//when used 4 row is assumed to be 0,0,0,1
	//4th column is used to store 
	//H[M t]
	// [0 1]
	struct MFtransform : public MFmat4
	{
		MFtransform() = default;
		//Ax,y,z,0 Bx,y,z,0	Cx,y,z,0 Tx,y,z,1
		MFtransform(const MFfloat& f00, const MFfloat& f01, const MFfloat& f02,
			const MFfloat& f10, const MFfloat& f11, const MFfloat& f12,
			const MFfloat& f20, const MFfloat& f21, const MFfloat& f22,
			const MFfloat& p30, const MFfloat& p31, const MFfloat& p32);
		//vectors are typically column major, vector ctor takes row major
		MFtransform(const MFvec3& a, const MFvec3& b, const MFvec3& c, const MFpoint3& p);
		//used for taking in 3 row vectors with implicit {0}1 4th row, translation is contained in V
		MFtransform(const MFvec4& a, const MFvec4& b, const MFvec4& c);
		MFvec3& operator[](const int& j) { return *reinterpret_cast<MFvec3*>(field[j]); };
		const MFvec3& operator[](const int& j) const { return *reinterpret_cast<const MFvec3*>(field[j]); };
		MFpoint3& GetTranslation() { return *reinterpret_cast<MFpoint3*>(field[3]); };
		const MFpoint3& GetTranslation() const { return *reinterpret_cast<const MFpoint3*>(field[3]); };
		void SetTranslation(const MFpoint3& translation);
	};

	//matrix manipulation
	MFtransform Inverse(const MFtransform& transform);
	MFmat3 NormalMatrix(const MFtransform& transform);//returns inverse tranpose 3x3 of the transform
	//transforms coord sys A to coord sys B
	MFtransform operator*(const MFtransform& hA, const MFtransform& hB);
	//represents the transformation of a point in somes local sys to sys A, h*p
	MFpoint3 operator*(const MFtransform& hA, const MFpoint3& p);
	//represents the transformation of a normal in sys B to sys A, v*h
	MFvec3 operator*(const MFvec3& n, const MFtransform& hB);
	//converts a vector created with a Y-up to a Z up
	MFpoint3 ConvertWorldUpZ(const MFpoint3& convert);
	//converts a vector created with a Z-up to a Y up
	MFpoint3 ConvertWorldUpY(const MFpoint3& convert);
	//xform identity - identity defined here  since xform->mat4 but not mat4->xform
	MFtransform Identity();
}