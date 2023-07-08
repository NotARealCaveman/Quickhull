#pragma once
#include <stdint.h>
#include <string>

namespace Manifest_Utility
{	
	//tbd
	using MFstring = std::string;
	///single character types
	using MFbool = bool;
	using MFchar = char;	
	//floating point types
	using MFhalf = uint16_t;//not currently implemented
	using MFfloat = float;
	using MFdouble = double;
	//signed integer types	
	using MFint8 = int8_t;
	using MFint16 = int16_t; 
	using MFint32 = int32_t;	
	using MFint64 = int64_t;
	using MFlong = MFint32;
	using MFllong = MFint64;
	//unsigned integer types
	using MFu8 = uint8_t;
	using MFu16 = uint16_t;
	using MFu32 = uint32_t;
	using MFu64 = uint64_t;
	using MFulong = MFu32;
	using MFullong = MFu64;
	using MFsize = size_t;
}