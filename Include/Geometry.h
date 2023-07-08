#pragma once
#include "Point3.h"
#include "Plane.h"


namespace Manifest_Math
{
	struct MFtriangle
	{
		MFpoint3 vertices[3];
	};
	MFvec3 SurfaceNormal(const MFtriangle& triangle);
	MFplane CalculateSurfacePlane(const MFtriangle& triangle);	
}