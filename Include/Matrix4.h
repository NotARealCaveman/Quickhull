#pragma once
#include "Vector4.h"
#include "Matrix3.h"



namespace Manifest_Math
{
	//memory layout
	//rows contain vectors X Y Z
	//columns contain compXYZ
	struct MFmat4
	{
	public:
		//[C][R] column major access
		MFfloat field[4][4];
		/* memory layout
		* f[0][0-2] = vX = vXx vXy vXz
		* f[1][0-2] = vY = vYx vYy vYz
		* f[2][0-2] = vZ = vZx vZy vZz
		*/
	public:
		MFmat4() = default;//identity matrix
		//takes in matrix values as column major
		//stores them in memory layout above
		MFmat4(
			MFfloat f00, MFfloat f01, MFfloat f02, MFfloat f03,
			MFfloat f10, MFfloat f11, MFfloat f12, MFfloat f13,
			MFfloat f20, MFfloat f21, MFfloat f22, MFfloat f23,
			MFfloat f30, MFfloat f31, MFfloat f32, MFfloat f33) :
			/// [00][01][02][03] [10][11][12][13]  [20][21][22][23]  [30][31][32][33]
			field{ 
			f00,f01,f02,f03,
				f10,f11,f12,f13,
				f20,f21,f22,f23,
				f30,f31,f32,f33 }
		{};
		//takes in column major vectors and stores as row major
		//row vUxyzw = x ==Col 0, y == Col 1, z == Col 2, w == Col 3
		//col vKabcd = a ==Row 0, b == Row 1, c == Row 2, d == Row 3
		MFmat4(const MFvec4& a, const MFvec4& b, const MFvec4& c, const MFvec4& d) :
			/// [00][01][02][03]  [10][11][12][13]  [20][21][22][23]   [30][31][32][33]
			field{ a.x,b.x,c.x,d.x,
			a.y,b.y,c.y,d.y,a.z,
			b.z,c.z,d.z ,
			a.w,b.w,c.w,d.w  }
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
		MFvec4& operator[](const int& colVec)
		{
			return *reinterpret_cast<MFvec4*>(field[colVec]);
		}
		const MFvec4& operator[](const int& colVec)const
		{
			return *reinterpret_cast<const MFvec4*>(field[colVec]);
		}			
		explicit operator MFmat3 () const
		{
			return MFmat3
			{
				*reinterpret_cast<const MFvec3*>(&field[0]),
				*reinterpret_cast<const MFvec3*>(&field[1]),
				*reinterpret_cast<const MFvec3*>(&field[2])
			};
		}
		//console shit
		friend std::ostream& operator<<(std::ostream& os, const MFmat4& inMat)
		{
			os << "Printing Matrix@ " << &inMat.field << " \n";
			os << inMat.field[0][0] << " " << inMat.field[1][0] << " " << inMat.field[2][0] << " " << inMat.field[3][0] << " \n";
			os << inMat.field[0][1] << " " << inMat.field[1][1] << " " << inMat.field[2][1] << " " << inMat.field[3][1] << " \n";
			os << inMat.field[0][2] << " " << inMat.field[1][2] << " " << inMat.field[2][2] << " " << inMat.field[3][2] << " \n";
			os << inMat.field[0][3] << " " << inMat.field[1][3] << " " << inMat.field[2][3] << " " << inMat.field[3][3] << " \n";

			return os;
		};

	};	
	
	MFmat4 Transpose(const MFmat4& mA);
	MFmat4 Inverse(const MFmat4& mA);
	MFfloat Determinate(const MFmat4& mA);
	MFmat4 operator*(const MFmat4& mA, const MFmat4& mB);
	MFvec4 operator*(const MFmat4 mA, const MFvec4& vB);	
	MFmat4 XYZtoXZY(const MFmat4& mA);
	MFmat4 XZYtoXYZ(const MFmat4& mA);
}