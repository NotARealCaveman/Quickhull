#pragma once
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cassert>

#define CONSOLE_COLOR_TEMPLATE () printf("\x1B[31m'SampleText'033[0m ");
#define LOG(COLOR,x) do { std::cerr <<"\x1B["<<COLOR<<"m"<< x<<"\033[0m" << "\n"; }while (0)

#ifndef _DEBUG //vs built in macro
#define DLOG(COLOR,x) 
#define RLOG(COLOR,x) do { std::cerr <<"\x1B["<<COLOR<<"m"<< x<<"\033[0m" << "\n"; }while (0)
#else
#define DLOG(COLOR,x)  do { std::cerr <<"\x1B["<<COLOR<<"m"<< x<<"\033[0m" << "\n"; }while (0)
#define RLOG(COLOR,x)
#endif


#ifndef DISABLE
#define DISABLE if(0)//disable code path
#endif // !DISABLE

#ifndef DEFAULT_BREAK//used in switch statments to indicate defualt case may be used in the future 
#define DEFAULT_BREAK default:break;
#endif // !DEFAULT_BREAK

 

