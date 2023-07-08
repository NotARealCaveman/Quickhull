#include "Line.h"

using namespace Manifest_Math;

const MFfloat Manifest_Math::DistancePointLine(const MFpoint3& q, const MFpoint3& l, const MFvec3& v)
{
	const auto& cUV = Cross(q - l, v);

	return std::sqrtf(Dot(cUV, cUV) / Dot(v, v));
};

const MFfloat Manifest_Math::NDistancePointLine(const MFpoint3& q, const MFpoint3& l, const MFvec3& v)
{
	const auto& cUV = Cross(q - l, v);

	return std::sqrtf(Dot(cUV, cUV));// / Dot(v, v)); v has unit length
};

const MFfloat Manifest_Math::DistanceLineLine(const MFpoint3& p1, const MFvec3 v1, const MFpoint3& p2, const MFvec3& v2)
{
	const auto& distP = p2 - p1;

	const auto& v12 = Dot(v1, v1);
	const auto& v22 = Dot(v2, v2);
	const auto& v1v2 = Dot(v1, v2);

	auto&& det = v1v2 * v1v2 - v12 * v22;

	if (std::fabs(det) > FLT_MIN)//a*b = 0 ; lines are parrallel
	{
		det = 1.0f / det;

		const auto& distPv1 = Dot(distP, v1);
		const auto& distPv2 = Dot(distP, v2);
		const auto& t1 = (v1v2 * distPv2 - v22 * distPv1) * det;
		const auto& t2 = (v12 * distPv2 - v1v2 * distPv1) * det;

		return Magnitude(distP + v2 * t2 - v1 * t1);
	}//lines are almost parallel
	const auto& cDistPv1 = Cross(distP, v1);

	return std::sqrt(Dot(cDistPv1, cDistPv1) / v12);
};

//Plücker Coordinate, Implicit lines
MFline::MFline(const MFfloat& vx, const MFfloat& vy, const MFfloat& vz, const MFfloat& mx, const MFfloat& my, const MFfloat& mz)
	:v{ vx,vy,vz }, m{ mx,my,mz }
{};

MFline::MFline(const MFvec3& _v, const MFvec3& _m)
	:v{_v}, m{_m}
{};

//transform l{v|m}->H[Mt] to l{v|m}A 
//v=p1-p2, vA=hAv
//m = p1xp2, pA=Mp+t
//mA = (Mp1+t) x (Mp2 + t)
//   = Mp1 x Mp2 + t x Mp2 - t x Mp1
// Mp1 x Mp2 transforms as the adj(M)
//t x Mp2 - t x Mp1 = t x (Mv)
//transformation of moment m to mA is
//mA = m*adj(M) + t x (Mv)
//due to memory layout of matrices,
//m*adj(M) is adj(M)*m for C-Major M3
MFline Manifest_Math::Transform(const MFline& l, const MFtransform& hA)
{
	const auto& adjugate = MFmat3
	{
		Cross(hA[1],hA[2]),
		Cross(hA[2],hA[0]),
		Cross(hA[0],hA[1])
	};
	const auto& t = hA.GetTranslation();

	const auto& v = hA * l.v;
	const auto& m = adjugate * l.m + Cross(t, v);

	return MFline{ v,m };
}
//if M is orthogonal, such that rows M and columns M are orthonormal
//adj(M)m can be rewritten as M*m
//m is treated as a column vector
//(MFmat3)H->0.037μs, adj(M)->0.199μs
MFline Manifest_Math::TransformOrthonormal(const MFline& l, const MFtransform& hA)
{		
	const auto& t = hA.GetTranslation();

	const auto& v = hA * l.v;
	const auto& m = static_cast<const MFmat3>(hA) * l.m + Cross(t, v);

	return MFline{ v,m };
}