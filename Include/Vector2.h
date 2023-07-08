#pragma once
#include "Core.h"

namespace Manifest_Math
{
	struct MFvec2
	{
		public:
			MFvec2() = default;
			MFvec2(const MFfloat& u);
			MFvec2(const MFfloat& _x, const MFfloat& _y);

			MFfloat& operator[](const int& i) { return ((&x)[i]); };
			const MFfloat& operator[](const int& i) const { return ((&x)[i]); };
			MFvec2& operator+=(const MFvec2& inVec)
			{
				x += inVec.x;
				y += inVec.y;				
				return (*this);
			};
			MFvec2& operator-=(const MFvec2& inVec)
			{
				x -= inVec.x;
				y -= inVec.y;				
				return (*this);
			
			};
			//default initialzed for neiche case uses as these containers come up
			MFfloat x=0;
			MFfloat y=0;
	};
	inline MFvec2 operator+(const MFvec2& componentVec, const MFvec2& inVec)
	{
		return MFvec2{ componentVec.x + inVec.x, componentVec.y + inVec.y};
	};
	inline MFvec2 operator-(const MFvec2& componentVec, const MFvec2& inVec)
	{
		return MFvec2{ componentVec.x - inVec.x, componentVec.y - inVec.y};
	};
	inline MFvec2 operator/(const MFvec2& dividendVec, const MFfloat& divisor)
	{
		const auto& s = 1.0f / divisor;
		return MFvec2{ dividendVec.x * s, dividendVec.y * s };
	}
	inline MFvec2 operator*(const MFvec2& multiplicandVec, const MFfloat& multiplier)
	{
		return MFvec2{ multiplicandVec.x * multiplier, multiplicandVec.y * multiplier };
	}
	inline MFvec2 operator-(const MFvec2& inverseVec)
	{
		return MFvec2{ -inverseVec.x, -inverseVec.y};
	};
	MFvec2 Normalize(const MFvec2& inVec);
	MFfloat Dot(const MFvec2& a, const MFvec2& b);
	
}
