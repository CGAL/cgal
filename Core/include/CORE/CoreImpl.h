// The following two lines only for MS Visual C++
#ifdef _MSC_VER
	#pragma warning(disable: 4291)
	#pragma warning(disable: 4800)
#endif

// condition preprocessor for inline function
#ifndef _DEBUG
	#define CORE_ENABLE_INLINES
	#ifdef _MSC_VER
		#define CORE_INLINE __forceinline 
	#else
		#define CORE_INLINE inline
	#endif
#else
	#define CORE_INLINE
#endif

// Macros for defining namespace
#define CORE_BEGIN_NAMESPACE	namespace CORE {
#define CORE_END_NAMESPACE	};

// include some common header files
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <cctype>
#include <climits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
