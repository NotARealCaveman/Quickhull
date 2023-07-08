#pragma once
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <memory>
#include <iostream>

#include "Typenames.h"

using namespace Manifest_Utility;

namespace Manifest_Math
{
	struct xoshiro256ss_state
	{
		private:
			void Seed(const MFu32& seed);
			static MFu64 rol64(MFu64 x, MFint32 k)
			{
				return (x << k) | (x >> (64 - k));
			}
			MFu64 states[4]; //data
			MFu64 currentState;
		public:
			xoshiro256ss_state(const MFu32& seed = 0)
			{			
				for(MFu32 shake = 0; shake < 5;++shake)
				Seed(seed);
			}
			//morph state and return new integer using xoshiro256ss algorithm	
			MFu64 Crunch();
			//returns new float between [0,1)*scale
			MFfloat CrunchFloat(const MFfloat& scale = 1);
			//returns new integer between [lower,upper]
			MFint64 CrunchRange(const MFint32& lower, const MFint32& upper);
			//returns new float between [lower,upper)
			MFfloat CrunchRangeFloat(const MFfloat& lower, const MFfloat& upper);
			
			MFu64 const& GetCurrentState() { return this->currentState; };

			static xoshiro256ss_state Random;
	};
}