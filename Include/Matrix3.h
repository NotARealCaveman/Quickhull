#pragma once
#include "DebugLogger.h"

#include "Vector3.h"

namespace Manifest_Math
{
	//memory layout
	//rows contain vectors X Y Z
	//columns contain compXYZ
	struct MFmat3
	{
	public:
		//[C][R] column major access
		MFfloat field[3][3];
		/* memory layout
		* f[0][0-2] = vX = vXx vXy vXz
		* f[1][0-2] = vY = vYx vYy vYz
		* f[2][0-2] = vZ = vZx vZy vZz
		*/
	public:
		MFmat3() = default;
		//takes in matrix values as column major
		//stores them in memory layout above
		MFmat3(MFfloat f00, MFfloat f01, MFfloat f02,
			MFfloat f10, MFfloat f11, MFfloat f12,
			MFfloat f20, MFfloat f21, MFfloat f22) :
			/// [00][01][02]  [10][11][12]  [20][21][22]
			field{ f00,f01,f02,f10,f11,f12,f20,f21,f22 }
		{};
		//takes in row major vectors and stores as column major
		//row vUxyz = x == Col 0, y == Col 1 , z == Col 2
		//col vKabc = a == Row 0, b == Row 1 , c == Row 2
		MFmat3(const MFvec3& a, const MFvec3& b, const MFvec3& c) :
			/// [00][01][02]  [10][11][12]  [20][21][22]
			field{ a.x,b.x,c.x,a.y,b.y,c.y,a.z,b.z,c.z  }
		{};

		//accessors	
		//takes in row major coordinates for sake of reading
		// 1 4 7 =>EX (0,1) = 4 != [0,1] == [1,0]
		// 2 5 8 =>EX (1,2) = 8 != [1,2] == [2,1]
		// 3 6 9 =>EX (2,1) = 6 != [2,1] == [1,2]
		MFfloat& operator()(const int& row, const int& column)
		{
			return field[column][row];
		}
		const MFfloat& operator()(const int& row, const int& column)const
		{
			return field[column][row];
		}
		MFvec3& operator[](const int& colVec)
		{
			return *reinterpret_cast<MFvec3*>(field[colVec]);
		}
		const MFvec3& operator[](const int& colVec)const
		{
			return *reinterpret_cast<const MFvec3*>(field[colVec]);
		}
		const MFmat3& operator*=(const MFfloat& s)
		{
			*reinterpret_cast<MFvec3*>(field[0]) *= s;
			*reinterpret_cast<MFvec3*>(field[1]) *= s;
			*reinterpret_cast<MFvec3*>(field[2]) *= s;
			return *this;
		}		
		//console shit
		friend std::ostream& operator<<(std::ostream& os, const MFmat3& inMat)
		{
			os <<  &inMat.field << " \n";
			os << inMat.field[0][0] << " " << inMat.field[1][0] << " " << inMat.field[2][0] << " \n";
			os << inMat.field[0][1] << " " << inMat.field[1][1] << " " << inMat.field[2][1] << " \n";
			os << inMat.field[0][2] << " " << inMat.field[1][2] << " " << inMat.field[2][2] << " \n";

			return os;
		};
	};		

	MFmat3 Transpose(const MFmat3& mA);
	MFmat3 Inverse(const MFmat3& mA);
	///transforms
	MFfloat Determinate(const MFmat3& mA);	
//rotates x axis about yz plane
	MFmat3 Rotation_X(const MFfloat& theta);
	//rotates y axis about xz plane
	MFmat3 Rotation_Y(const MFfloat& theta);
//rotates z axis about xy plane
	MFmat3 Rotation_Z(const MFfloat& theta);
//rotates V about an axis
	MFmat3 Rotation(const MFfloat& theta, const MFvec3& axis);
//rotates V about a plane perpendicular to an axis 
	MFmat3 Reflect(const MFvec3& axis);
//scale V by some x,y,z
	MFmat3 Scale(const MFvec3& scale);
//scales V by some value s along an axis
	MFmat3 Scale_AxisAligned(const MFfloat& scale, const MFvec3& axis);
//skews V by the projection of V onto axisUp by theta(between V and VrejAxis) in the direction of the axis
	MFmat3 Skew(const MFfloat& theta, const MFvec3& axis, const MFvec3& axisPerpendicular);
	MFmat3 SkewSymmetric(const MFvec3& a);
//not used but included for completeness
//duplicating a refelection /rotating by 180 degrees is an involution matrix
	MFmat3 Involution(const MFvec3& axis);
//applies Gram-Schmidt process to orthogonalize a matrix - original matrix is overriden
	MFmat3 GramSchmidt(const MFvec3& a, const MFvec3& b, const MFvec3& c);

	inline MFmat3 operator*(const MFmat3& mA, const MFmat3& mB)
	{
		return MFmat3
		{
			mA(0,0) * mB(0,0) + mA(0,1) * mB(1,0) + mA(0,2) * mB(2,0), //M(0,0) = M[0][0]
			mA(1,0) * mB(0,0) + mA(1,1) * mB(1,0) + mA(1,2) * mB(2,0), //M(1,0) = M[0][1]
			mA(2,0) * mB(0,0) + mA(2,1) * mB(1,0) + mA(2,2) * mB(2,0), //M(2,0) = M[0][2]

			mA(0,0) * mB(0,1) + mA(0,1) * mB(1,1) + mA(0,2) * mB(2,1), //M(0,1) = M[1][0]
			mA(1,0) * mB(0,1) + mA(1,1) * mB(1,1) + mA(1,2) * mB(2,1), //M(1,1) = M[1][1]
			mA(2,0) * mB(0,1) + mA(2,1) * mB(1,1) + mA(2,2) * mB(2,1), //M(2,1) = M[1][2]	 

			mA(0,0) * mB(0,2) + mA(0,1) * mB(1,2) + mA(0,2) * mB(2,2), //M(0,2) = M[2][0]
			mA(1,0) * mB(0,2) + mA(1,1) * mB(1,2) + mA(1,2) * mB(2,2), //M(1,2) = M[2][1]
			mA(2,0) * mB(0,2) + mA(2,1) * mB(1,2) + mA(2,2) * mB(2,2), //M(2,2) = M[2][2]
		};
	};
	inline MFmat3 operator+(const MFmat3& mA, const MFmat3& mB)
	{
		return MFmat3
		{
			mA(0,0)+mB(0,0),mA(1,0) + mB(1,0),mA(2,0) + mB(2,0),
			mA(0,1) + mB(0,1),mA(1,1) + mB(1,1),mA(2,1) + mB(2,1),
			mA(0,2) + mB(0,2),mA(1,2) + mB(1,2),mA(2,2) + mB(2,2)
		};
	};

	inline MFvec3 operator*(const MFmat3& mA, const MFvec3& vB)
	{
		return MFvec3
		{													  //N=x/y/z element pos
			mA(0,0) * vB.x + mA(0,1) * vB.y + mA(0,2) * vB.z, //MVx = mAc0rN*vBN
			mA(1,0) * vB.x + mA(1,1) * vB.y + mA(1,2) * vB.z, //MVx = mAc1rN*vBN
			mA(2,0) * vB.x + mA(2,1) * vB.y + mA(2,2) * vB.z  //MVx = mAc2rN*vBN
		};
	};

	inline MFvec3 Transform(const MFmat3& mA, const MFvec3& vB)
	{
		return MFvec3
		{
		 mA[0][0] * vB.x + mA[0][1] * vB.y + mA[0][2] * vB.z ,
		 mA[1][0] * vB.x + mA[1][1] * vB.y + mA[1][2] * vB.z,
		 mA[2][0] * vB.x + mA[2][1] * vB.y + mA[2][2] * vB.z
		};
	};

	//defines three vectors which are orthonormal to one another
	//tangent space can be used as an example where the space TBN
	//represents the forward N, the left T, and the up B
	typedef MFmat3 OrthonormalBase;
}