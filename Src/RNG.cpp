#include "RNG.h"

using namespace Manifest_Math;

xoshiro256ss_state xoshiro256ss_state::Random{};

void xoshiro256ss_state::Seed(const MFu32& seed)
{
	if (!seed)
		srand(time(NULL)); //set seeds seed
	else
		srand(seed);
	//seed each part of state
	for (MFu32 i = 0; i < 4; ++i)
		states[i] = rand() ^ Crunch();
}

//return random value and modify state object
MFu64 xoshiro256ss_state::Crunch()
{	
	const MFu64 result = rol64(states[1] * 5, 7) * 9;
	const MFu64 t = states[1] << 17;

	states[2] ^= states[0];
	states[3] ^= states[1];
	states[1] ^= states[2];
	states[0] ^= states[3];

	states[3] = rol64(states[3], 45);
	states[2] ^= t;
		 
	currentState = result;

	return result;
}
//returns random float and internally crunches numbers
MFfloat xoshiro256ss_state::CrunchFloat(const MFfloat& scale)
{
	auto crunch = static_cast<MFint32>(Crunch());//capture random, signed 32 int
	auto mask = crunch >> (sizeof(crunch)*8 - 1);//abs mask
	crunch = (crunch + mask) ^ mask;//remove sign if present
	//return random float [0,1) * scale
	return static_cast<MFfloat>(crunch) / static_cast<MFfloat>(INT32_MAX) * scale; 
}

//returns a random number clamped between the upper and lower range
//internally calls the crunch function to shift the RNG state
//ranged clamped to u32s, can be changed to 64s if needed later, i dont feel like testing bandwidth but why convert an extras 4 bytes
//signed type used here as range can be negative
MFint64 xoshiro256ss_state::CrunchRange(const MFint32& lower, const MFint32& upper)
{	
	auto crunch = Crunch();		
	//sign must be removed for range based crunch to handle negative numbers	
	crunch = crunch >> (crunch >> 63);
	auto upper64 = static_cast<MFint64>(upper);
	auto lower64 = static_cast<MFint64>(lower);
	auto range = (upper64 - lower64 + 1);
	auto crunchRange = static_cast<MFint64>(crunch) % range;
	auto shiftedCrunch = crunchRange + lower64;
	return shiftedCrunch;
}

MFfloat xoshiro256ss_state::CrunchRangeFloat(const MFfloat& lower, const MFfloat& upper)
{
	//returns a the lower end + the scaled decimal range
	return lower + (CrunchFloat() * (upper - lower));
}