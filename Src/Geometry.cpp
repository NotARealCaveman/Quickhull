#include "Geometry.h"

using namespace Manifest_Math;

//TRIANGLES
MFvec3 Manifest_Math::SurfaceNormal(const MFtriangle& triangle)
{		
	return Normalize(Cross(triangle.vertices[1] - triangle.vertices[0], triangle.vertices[2] - triangle.vertices[0]));
}

MFplane Manifest_Math::CalculateSurfacePlane(const MFtriangle& triangle)
{
	auto normal{ Cross(triangle.vertices[1] - triangle.vertices[0],triangle.vertices[2] - triangle.vertices[0]) * 0.5f };
	auto offset{ Dot(-normal,triangle.vertices[0]) };
	return { normal,offset };
}
