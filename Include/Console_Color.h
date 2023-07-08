#pragma once
#include "DebugLogger.h"

//used for colored codes on windows 
#ifdef _WIN64 
#define	NOMINMAX //turns off windows max function - redef w/ limits::max
#include <Windows.h>
#define STD_OUTPUT_HANDLE_manifest ((DWORD)-11)
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING_manifest 0x0004
void Enable_ASCII_Color()
{
	#ifdef _WIN64 
		HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE_manifest);
		DWORD dwMode = 0;
		GetConsoleMode(hOut, &dwMode);
		SetConsoleMode(hOut, dwMode | ENABLE_VIRTUAL_TERMINAL_PROCESSING_manifest);
	#endif  
}
#endif // _WIN64 

#define WINDOWS_COLOR_CONSOLE Enable_ASCII_Color();