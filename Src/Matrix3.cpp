#include "Matrix3.h"

using namespace Manifest_Math;

//Matrix functions
MFmat3 Manifest_Math::Transpose(const MFmat3& mA)
{
	return MFmat3
	{
		mA(0,0),mA(0,1),mA(0,2),
		mA(1,0),mA(1,1),mA(1,2),
		mA(2,0),mA(2,1),mA(2,2),
	};
}
//returns the "magnitude" of Ma
//|M| is a scalar = hypervolume of nD parallelotope
//uses the Scalar Triple Product to calculate the determinate
//[a b c] aXb = volume of parallelogram AB, aXb* c = volume of parallelogram scaled to D3 by c, making it the volume of a parallelpiped
MFfloat Manifest_Math::Determinate(const MFmat3& mA)
{
	return Dot(Cross(mA[0], mA[1]), mA[2]);//aXb*c
}

//retruns the inverse of M3
//M^-1=Madj^T/|M|
//M^-1=1/|M|*Madj
//M^-1=1/([a,b,c])(MbXMc, McXMa, MaXMb)
MFmat3 Manifest_Math::Inverse(const MFmat3& mA)
{
	//|M| of M3 = Volume of the parallelpiped
	//grab column vectors of M and inverse determinate
	const auto& a = mA[0];
	const auto& b = mA[1];
	const auto& c = mA[2];	

	//create row vectors of C(M) and multiply by |M^-1|
	//Cofactor of M3 is the cross product of the alternate vectors of the row position
	//C(M) * |M|^-1 will build the Madj row vectors as M^-1
	auto r0 = Cross(b, c);
	auto r1 = Cross(c, a);
	auto r2 = Cross(a, b);
	//|M|^-1 in this case is 1/[a b c] 
	const auto& iDet = 1/Dot(r2,c);
	r0 *= iDet;
	r1 *= iDet;
	r2 *= iDet;

	//this allows vectors to be plugged in rather than component by component
	return MFmat3{ r0,r1,r2 };//calls row major vector constructor
}
///transformations
MFmat3 Manifest_Math::Rotation_X(const MFfloat& theta)
{
	const auto& c = std::cosf(theta);
	const auto& s = std::sinf(theta);

	return MFmat3
	{
		MFvec3{1,0,0},
		MFvec3{0,c,-s},
		MFvec3{0,-s,c}
	};
}
MFmat3 Manifest_Math::Rotation_Y(const MFfloat& theta)
{
	const auto& c = std::cosf(theta);
	const auto& s = std::sinf(theta);

	return MFmat3
	{
		MFvec3{c,0,s},
		MFvec3{0,1,0},
		MFvec3{-s,0,c}
	};
}
MFmat3 Manifest_Math::Rotation_Z(const MFfloat& theta)
{
	const auto& c = std::cosf(theta);
	const auto& s = std::sinf(theta);

	return MFmat3
	{
		MFvec3{c,-s,0},
		MFvec3{c,s,0},
		MFvec3{0,0,1}
	};
}
MFmat3 Manifest_Math::Rotation(const MFfloat& theta, const MFvec3& axis)
{
	const auto& c = std::cosf(theta);
	const auto& s = std::sinf(theta);
	const auto& d = 1.0f - c;//moves from V_|_a to v'-a

	const auto& dAx = axis.x * d;		
	const auto& AxAy = dAx * axis.y;
	const auto& dAy = axis.y * d;
	const auto& AxAz = dAx * axis.z;
	const auto& dAz = axis.z * d;
	const auto& AyAz = dAy * axis.z;
	
	return MFmat3
	{
		c+dAx*axis.x, AxAy + s * axis.z, AxAz - s * axis.y,
		AxAy-s*axis.z, c+ dAy * axis.y, AyAz + s * axis.x,
		AxAz + s * axis.y, AyAz - s * axis.x, c + dAz * axis.z
	};
}
MFmat3 Manifest_Math::Reflect(const MFvec3& axis)
{
	const auto& x = axis.x * -2.0f;
	const auto& AxAy = x * axis.y;
	const auto& y = axis.y * -2.0f;
	const auto& AxAz = x * axis.y;	
	const auto& z = axis.z * -2.0f;
	const auto& AyAz = y * axis.z;

	return MFmat3
	{
		x * axis.x + 1.0f, AxAy,AxAz,
		AxAy, y * axis.y + 1.0f, AyAz,
		AxAz, AyAz, z * axis.z + 1.0f
	};
}


MFmat3 Manifest_Math::Scale(const MFvec3& scale)
{
	return MFmat3
	{
		scale.x,0.0f,0.0f,
		0.0f,scale.y,0.0f,
		0.0f,0.0f, scale.z
	};
}

MFmat3 Manifest_Math::Scale_AxisAligned(const MFfloat& scale, const MFvec3& axis)
{
	const auto& s = scale - 1.0f;
	const auto& x = axis.x * s;
	const auto& AxAy = x * axis.y;
	const auto& y = axis.y * s;
	const auto& AxAz = x * axis.z;
	const auto& z = axis.z * s;
	const auto& AyAz = y * axis.z;

	return MFmat3
	{
		x * axis.x + 1.0f, AxAy, AxAz,
		AxAy, y * axis.y + 1.0f, AyAz,
		AxAz,AyAz, z * axis.z + 1.0f
	};
}

MFmat3 Manifest_Math::Skew(const MFfloat& theta, const MFvec3& axis, const MFvec3& axisPerpendicular)
{
	const auto& t = std::tan(theta);
	const MFfloat& x = axis.x * t;
	const MFfloat& y = axis.y * t;
	const MFfloat& z = axis.z * t;
	return MFmat3
	{
		x * axisPerpendicular.x + 1.0f, y * axisPerpendicular.x, z * axisPerpendicular.x,
		x * axisPerpendicular.y , y * axisPerpendicular.y + 1.0f, z * axisPerpendicular.y,
		x * axisPerpendicular.z , y * axisPerpendicular.z, z * axisPerpendicular.z + 1.0f,
	};
}

//returns a matrix such that if a = ssMa, ssMa*b=axb
MFmat3 Manifest_Math::SkewSymmetric(const MFvec3& a)
{
	return MFmat3
	{
		0.0f,a.z,-a.y,
		-a.z,0.0f,a.x,
		a.y,-a.x,0.0f
	};
}

MFmat3 Manifest_Math::Involution(const MFvec3& axis)
{
	const auto& x = axis.x * 2.0f;
	const auto& AxAy = x * axis.y;
	const auto& y = axis.y * 2.0f;
	const auto& AxAz = x * axis.z;
	const auto& z = axis.z * 2.0f;
	const auto& AyAz = y * axis.z;

	return MFmat3
	{
		x * axis.x - 1.0f, AxAy,AxAz,
		AxAy, y * axis.y - 1.0f, AyAz,
		AxAz, AyAz, z * axis.z - 1.0f
	};
}
//takes in a non orthonormal base with u1 set to the x axis. computes u2 and u3 internally depending on the desired v1 and v2 stored in the Y-Z columns
//u1 = v1|u2=v2-(u1_|_v2)|u3=v3-(u1_|_v3)-(u2_|_v3)
MFmat3 Manifest_Math::GramSchmidt(const MFvec3& a, const MFvec3& b, const MFvec3& c)
{
	MFvec3 u0 = a;
	MFvec3 u1 = Normalize(b - Project(b, u0));
	MFvec3 u2 = Normalize(c - Project(c, u0) - Project(c, u1));	
	return MFmat3{ u0.x,u0.y,u0.z,
	u1.x,u1.y,u1.z,
	 u2.x,u2.y,u2.z }; 	
}