#include <vector>

#include "Quickhull.h"
#include "RNG.h"

using namespace Manifest_Math;

const std::vector<MFpoint3> pointCloud
{
	{-1, -1, -1},	{-1, 1, -1},
	{-1, 1, -1},	{1, 1, -1},
	{ 1, -1, -1 },	{-1, -1, 1},
	{1, -1, 1},		{1, 1, 1},
	{1, 1, 1},		{-1, 1, 1},
	{-1, -1, -1},	{1, -1, -1},
	{1, -1, 1},		{-1, -1, 1},
	{1, -1, -1},	{1, 1, -1},
	{1, 1, 1},		{1, -1, 1},
	{1, 1, -1},		{-1, 1, -1},
	{-1, 1, 1},		{1, 1, 1}
};

int main()
{
	auto hull{ QuickHull(pointCloud) };	
	LOG(45, "Press q to quit");
	LOG(46, "Press a to add point to hull");	
	const auto lowerBound{ -10 };
	const auto upperBound{ 10 };
	std::vector<MFpoint3>range{ {lowerBound},{upperBound} };
	hull.CONVEXITY_EPSILON = EPSILON_3D(range);
	xoshiro256ss_state ss;
	char key{ ' ' };
	std::cin >> key;	
	while (key != 'q' )
	{
		if (key == 'a')
		{			
			MFpoint3 newPoint{ ss.CrunchRangeFloat(lowerBound,upperBound),ss.CrunchRangeFloat(lowerBound,upperBound),ss.CrunchRangeFloat(lowerBound,upperBound) };
			
			AddPointToHull(hull, newPoint);
		}
		std::cin >> key;
	}
}