#pragma once
#include <array>

#include "Point3.h"

namespace Manifest_Math
{
	template<typename T>
	struct Simplex_T
	{
		Simplex_T() = default;
		Simplex_T(std::initializer_list<T> _points)
		{
			*this = _points;
		}
		T& operator[](const MFint32& index) { return points[index]; };
		const T& operator[](const MFint32& index) const { return points[index]; };
		Simplex_T& operator=(std::initializer_list<T> _points) 
		{
			size = _points.size();
			std::transform(_points.begin(), _points.end(), points.begin(), [&](const T& point) {return point; });

			return *this;
		}				
		void PushFront(const T& point)
		{
			points = { point, points[0], points[1], points[2] };
			size = std::min(size + 1, 4);
		}	
		void PushBack(const T& point)
		{
			points[size] = point;
			size = std::min(size + 1, 4);			
		}

		std::array<T, 4> points{ 0,0,0,0 };
		MFu8 size{ 0 };
		
	};
	//TO CHANGE WITH ADAPTIVE EPSILON
	constexpr auto SIMPLEX_EPSILON{ 1e-7 };
	MFbool PositiveHalfSpace(const MFvec3& direction, const MFvec3& AO);
	template<typename T>
	MFbool NextSimplex(Simplex_T<T>& simplex, MFvec3& direction)
	{
		const auto& Line = [&](MFvec3& direction)->MFbool
		{
			const auto& s0{ simplex.points[0] };
			const auto& s1{ simplex.points[1] };
			const auto& a = s0.point;
			const auto& b = s1.point;

			const auto ab = b - a;
			const auto ao = -a;

			if (PositiveHalfSpace(ab, ao))
			{
				direction = Normalize(Cross(Cross(ab, ao), ab));
			}
			else
			{
				simplex = { s0 };
				direction = Normalize(ao);
			}

			return false;
		};
		const auto& Triangle = [&](MFvec3& direction)->MFbool
		{
			const auto& s0{ simplex.points[0] };
			const auto& s1{ simplex.points[1] };
			const auto& s2{ simplex.points[2] };
			const auto& a{ s0.point };
			const auto& b{ s1.point };
			const auto& c{ s2.point };

			const auto ab = b - a;
			const auto ac = c - a;
			const auto ao = -a;

			auto abc = Cross(ab, ac);
			if (PositiveHalfSpace(Cross(abc, ac), ao))
			{
				if (PositiveHalfSpace(ac, ao))
				{
					simplex = { s0,s2 };
					direction = Normalize(Cross(Cross(ac, ao), ac));
				}
				else
				{
					simplex = { s0,s1 };
					return Line(direction);
				}
			}
			else
			{
				if (PositiveHalfSpace(Cross(ab, abc), ao))
				{
					simplex = { s0,s1 };
					return Line(direction);
				}
				else
				{
					if (PositiveHalfSpace(abc, ao))
					{
						direction = Normalize(abc);
					}
					else
					{
						simplex = { s0,s2,s1 };
						direction = -abc;
					}
				}
			}

			//triangle contains origin, use surface normal as new direction
			direction = Normalize(Cross(ab, ac));

			return false;
		};
		const auto& Tetrahedron = [&](MFvec3& direction)->MFbool
		{
			const auto& s0{ simplex.points[0] };
			const auto& s1{ simplex.points[1] };
			const auto& s2{ simplex.points[2] };
			const auto& s3{ simplex.points[3] };
			const auto& a = s0.point;
			const auto& b = s1.point;
			const auto& c = s2.point;
			const auto& d = s3.point;

			const auto ab = b - a;
			const auto ac = c - a;
			const auto ad = d - a;
			const auto ao = -a;

			//calculate surface normals of apex faces
			const auto abc = Cross(ab, ac);
			const auto acd = Cross(ac, ad);
			const auto adb = Cross(ad, ab);
			//check positive half spaces - if origin lies on positive halfspace of the triangle then use this triangle as the new base simplex 
			if (PositiveHalfSpace(abc, ao))
			{
				simplex = { s0,s1,s2 };
				return Triangle(direction);
			}
			if (PositiveHalfSpace(acd, ao))
			{
				simplex = { s0,s2,s3 };
				return Triangle(direction);
			}
			if (PositiveHalfSpace(adb, ao))
			{
				simplex = { s0,s3,s1 };
				return Triangle(direction);
			}

			return true;
		};

		switch (simplex.size)
		{
		case 2:
			return Line(direction);
		case 3:
			return Triangle(direction);
		case 4:
			return Tetrahedron(direction);
		}

		return false;
	}

}