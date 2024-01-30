#include "Transform.h"

using namespace Manifest_Math;

//Ax,y,z,0 Bx,y,z,0	Cx,y,z,0 Tx,y,z,1
MFtransform::MFtransform(const MFfloat& f00, const MFfloat& f01, const MFfloat& f02,
	const MFfloat& f10, const MFfloat& f11, const MFfloat& f12,
	const MFfloat& f20, const MFfloat& f21, const MFfloat& f22,
	const MFfloat& p30, const MFfloat& p31, const MFfloat& p32) :	
	MFmat4{ f00, f01, f02, 0.0f, 
			f10, f11, f12, 0.0f, 
			f20, f21, f22, 0.0f, 
			p30, p31, p32, 1.0f }
{};
//vectors are typically column major, vector ctor takes row major
MFtransform::MFtransform(const MFvec3& a, const MFvec3& b, const MFvec3& c, const MFpoint3& p) :
	MFmat4{	a.x,a.y,a.z,0.0f,
			b.x,b.y,b.z,0.0f,
			c.x,c.y,c.z,0.0f,
			p.x,p.y,p.z,1.0f }
{};
//used for taking in 3 row vectors with implicit {0}1 4th row, translation is contained in V
MFtransform::MFtransform(const MFvec4& a, const MFvec4& b, const MFvec4& c) : MFmat4{ a,b,c,{{0.0f},1.0f} }
{};



void MFtransform::SetTranslation(const MFpoint3& translation)
{
	field[3][0] = translation.x;
	field[3][1] = translation.y;
	field[3][2] = translation.z;
}


//retruns the inverse of M4
//see hAtrix4.h for details
//implied homogenous vector {0},1 ; calculation optimizations can be used
MFtransform Manifest_Math::Inverse(const MFtransform& hA)
{
	//grab 3D C-Vectors of M4x3
	const MFvec3& a = reinterpret_cast<const MFvec3&>(hA[0]);
	const MFvec3& b = reinterpret_cast<const MFvec3&>(hA[1]);
	const MFvec3& c = reinterpret_cast<const MFvec3&>(hA[2]);
	const MFvec3& d = reinterpret_cast<const MFvec3&>(hA[3]);

	//create 0 vector once, when crossed in 3d all vectors are parallel to 0,0,0
	//Cross product of parallel vectors is a 0 vector
	MFvec3 s = Cross(a, b);//returns vector perpendicular to vectors a&b
	MFvec3 t = Cross(c, d);//returns vector perpendicular to vectors c&d

	MFfloat iDet = 1.0f / Dot(s, c);
	//create inverse scalars
	s *= iDet;
	t *= iDet;
	const MFvec3 v = c * iDet;
	//using inverse scalars, calculate inverse row vectors
	const MFvec3 r0 = Cross(b, v);
	const MFvec3 r1 = Cross(v, a); 
	
	return MFtransform
	{
		r0.x,r1.x,s.x,
		r0.y,r1.y,s.y,
		r0.z,r1.z,s.z,
		-Dot(b,t),Dot(a,t),-Dot(d,s)
	};
}

//applies homogenous inverse optimzation with tranpose construction of row values into column vectors
//used to construct a normal matrix from a model matrix transform
MFmat3 Manifest_Math::NormalMatrix(const MFtransform& hA)
{
	//grab 3D C-Vectors of M4x3
	const MFvec3& a = reinterpret_cast<const MFvec3&>(hA[0]);
	const MFvec3& b = reinterpret_cast<const MFvec3&>(hA[1]);
	const MFvec3& c = reinterpret_cast<const MFvec3&>(hA[2]);
	const MFvec3& d = reinterpret_cast<const MFvec3&>(hA[3]);
	//4th row vector of M4 is assumed to be
	//const MFvec4& p{ 0.0,0.0,0.0,1.0f}; 


	MFvec3 s = Cross(a, b);//returns vector perpendicular to vectors a&b
	MFvec3 t = Cross(c, d);//returns vector perpendicular to vectors c&d	
	

	MFfloat iDet = 1.0f / Dot(s, c);// +Dot(t, u); t*0
	//create inverse scalars
	s *= iDet;
	t *= iDet;
	//0=u *= iDet;
	const MFvec3 v = c * iDet;
	//using inverse scalars, calculate inverse row vectors
	const MFvec3 r0 = Cross(b, v);// + t * p.y; +0
	const MFvec3 r1 = Cross(v, a);// - t * p.x; -0		

	//returns the tranpose of the inverse
	//vectors are row vectors, ctor expects column major coordinates
	return MFmat3
	{
		r0.x,r0.y,r0.z,
		r1.x,r1.y,r1.z,
		s.x,s.y,s.z 
	};
}

MFtransform Manifest_Math::operator*(const MFtransform& hA, const MFtransform& hB)
{
	return MFtransform
	{
		//Column 0
		hA(0, 0) * hB(0, 0) + hA(0, 1) * hB(1, 0) + hA(0, 2) * hB(2, 0) + hA(0, 3) * hB(3, 0),//Row 0
		hA(1, 0) * hB(0, 0) + hA(1, 1) * hB(1, 0) + hA(1, 2) * hB(2, 0) + hA(1, 3) * hB(3, 0), //Row 1
		hA(2, 0) * hB(0, 0) + hA(2, 1) * hB(1, 0) + hA(2, 2) * hB(2, 0) + hA(2, 3) * hB(3, 0), //Row 2		
		//Column 1
		hA(0, 0) * hB(0, 1) + hA(0, 1) * hB(1, 1) + hA(0, 2) * hB(2, 1) + hA(0, 3) * hB(3, 1),//Row 0
		hA(1, 0) * hB(0, 1) + hA(1, 1) * hB(1, 1) + hA(1, 2) * hB(2, 1) + hA(1, 3) * hB(3, 1), //Row 1
		hA(2, 0) * hB(0, 1) + hA(2, 1) * hB(1, 1) + hA(2, 2) * hB(2, 1) + hA(2, 3) * hB(3, 1), //Row 2		
		//Column 2
		hA(0, 0) * hB(0, 2) + hA(0, 1) * hB(1, 2) + hA(0, 2) * hB(2, 2) + hA(0, 3) * hB(3, 2),//Row 0
		hA(1, 0) * hB(0, 2) + hA(1, 1) * hB(1, 2) + hA(1, 2) * hB(2, 2) + hA(1, 3) * hB(3, 2), //Row 1
		hA(2, 0) * hB(0, 2) + hA(2, 1) * hB(1, 2) + hA(2, 2) * hB(2, 2) + hA(2, 3) * hB(3, 2), //Row 2		
		//Column 3
		hA(0, 0) * hB(0, 3) + hA(0, 1) * hB(1, 3) + hA(0, 2) * hB(2, 3) + hA(0, 3) * hB(3, 3), //Row 0
		hA(1, 0) * hB(0, 3) + hA(1, 1) * hB(1, 3) + hA(1, 2) * hB(2, 3) + hA(1, 3) * hB(3, 3), //Row 1
		hA(2, 0) * hB(0, 3) + hA(2, 1) * hB(1, 3) + hA(2, 2) * hB(2, 3) + hA(2, 3) * hB(3, 3) //Row 2		
	};
}


MFpoint3 Manifest_Math::operator*(const MFtransform& hA, const MFpoint3& p)
{
	return MFpoint3
	{
		hA(0,0) * p.x + hA(0,1) * p.y + hA(0,2) * p.z + hA(0,3),
		hA(1,0) * p.x + hA(1,1) * p.y + hA(1,2) * p.z + hA(1,3),
		hA(2,0) * p.x + hA(2,1) * p.y + hA(2,2) * p.z + hA(2,3)
	};
}

MFpoint3 Manifest_Math::operator*(const MFtransform& hA, const MFvec3& v)
{
	return MFpoint3
	{
		hA(0,0) * v.x + hA(0,1) * v.y + hA(0,2) * v.z,
		hA(1,0) * v.x + hA(1,1) * v.y + hA(1,2) * v.z,
		hA(2,0) * v.x + hA(2,1) * v.y + hA(2,2) * v.z
	};
}

MFvec3 Manifest_Math::operator*(const MFvec3& n, const MFtransform& hB)
{
	return MFvec3
	{
		n.x * hB(0,0) + n.y * hB(1,0) + n.z * hB(2,0),
		n.x * hB(0,1) + n.y * hB(1,1) + n.z * hB(2,1),
		n.x * hB(0,2) + n.y * hB(1,2) + n.z * hB(2,2)
	};
}

MFpoint3 Manifest_Math::ConvertWorldUpZ(const MFpoint3& convert)
{
	const auto zUp = MFtransform
	{
		{1.0f,0.0f,0.0f},
		{0.0f,0.0f,1.0f},
		{0.0f,-1.0f,0.0},
		{0.0f}
	};

	return zUp * convert;
}

MFpoint3 Manifest_Math::ConvertWorldUpY(const MFpoint3& convert)
{
	const auto zUp = MFtransform
	{
		{1.0f,0.0f,0.0f},
		{0.0f,0.0f,-1.0f},
		{0.0f,1.0f,0.0},
		{0.0f}
	};

	return zUp * convert;
}

MFtransform Manifest_Math::Identity()
{
	return MFtransform
	{
		1.0f,0.0f,0.0f,
		0.0f,1.0f,0.0f,
		0.0f,0.0f,1.0f,
		0.0f,0.0f,0.0f
	};
}