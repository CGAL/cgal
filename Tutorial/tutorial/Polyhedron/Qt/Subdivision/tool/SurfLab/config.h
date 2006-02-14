// ======================================================================
//
// Copyright (c) 2001 SurfLab of CISE of University of Florida
//
// File          : src/config/config.h
// Description   : configuration functions and macros
// Creation_date : 19 Aug 2001
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

// $Id$

#ifndef _CONFIG_H_08192001
#define _CONFIG_H_08192001

// ======================================================================
// Namespace macro
// ======================================================================
#if defined(_USE_NAMESPACE)

#  define SURFLAB_BEGIN_NAMESPACE namespace SurfLab {
#  define SURFLAB_END_NAMESPACE }
#  define SURFLAB SurfLab
#  define USE_NAMESPACE_SURFLAB using namespace SurfLab;

#  define BBMATH_BEGIN_NAMESPACE namespace BBMath {
#  define BBMATH_END_NAMESPACE }
#  define BBMATH BBMath
#  define USE_NAMESPACE_BBMATH using namespace BBMath;

#  define OGL_BEGIN_NAMESPACE namespace OpenGL {
#  define OGL_END_NAMESPACE }
#  define OGL OpenGL
#  define USE_NAMESPACE_OGL using namespace OpenGL;

#else

#  define SURFLAB_BEGIN_NAMESPACE
#  define SURFLAB_END_NAMESPACE
#  define SURFLAB
#  define USE_NAMESPACE_SURFLAB

#  define BBMATH_BEGIN_NAMESPACE
#  define BBMATH_END_NAMESPACE
#  define BBMATH
#  define USE_NAMESPACE_BBMATH

#  define OGL_BEGIN_NAMESPACE
#  define OGL_END_NAMESPACE
#  define OGL
#  define USE_NAMESPACE_OGL

#endif


// ======================================================================
// Debuging functions
// ======================================================================
#if defined(_NDEBUG)
#  define NOT_IMPLEMENTED() ((void)0)
#else
#  define NOT_IMPLEMENTED() SL_CONF::not_implemented(__FILE__, __LINE__)
#endif

#if defined(_NO_ASSERTIONS) || defined(_NDEBUG)
#  define ASSERTION(EX) ((void)0)
#  define ASSERTION_MSG(EX,MSG) ((void)0)
#  define ASSERTION_CODE(CODE)
#else
#  define ASSERTION(EX) \
   ((EX) ? ((void)0) : SL_CONF::assertion_fail( # EX , __FILE__, __LINE__, 0))
#  define ASSERTION_MSG(EX,MSG) \
   ((EX) ? ((void)0) : SL_CONF::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define ASSERTION_CODE(CODE) CODE
#endif

#if defined(_NO_WARNINGS) || defined(_NDEBUG)
#  define WARNING(EX) ((void)0)
#  define WARNING_MSG(EX,MSG) ((void)0)
#  define WARNING_CODE(CODE)
#else
#  define WARNING(EX) \
   ((EX) ? ((void)0) : SL_CONF::warning_fail( # EX , __FILE__, __LINE__, 0))
#  define WARNING_MSG(EX,MSG) \
   ((EX) ? ((void)0) : SL_CONF::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define WARNING_CODE(CODE) CODE
#endif


// ----------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <cassert>

struct SL_CONF {
  ///
  inline static void not_implemented(const char*,int);
  ///
  inline static void assertion_fail (const char*,const char*,int,const char*);
  ///
  inline static void warning_fail   (const char*,const char*,int,const char*);
};

void SL_CONF::not_implemented(const char* file, int line) {
  std::cerr << "Error: not implemented!" << std::endl
	    << "File(Line): " << file << " ("<< line << ")" << std::endl;
  std::abort();
}

void SL_CONF::assertion_fail (const char* expr, const char* file, 
			   int line, const char* msg) {
  std::cerr << "Error: assertion violation!" << std::endl
	    << "Expr: " << expr << std::endl
	    << "File(Line): " << file << " (" << line << ")" << std::endl;
  if (msg != 0) std::cerr << "Explanation:" << msg << std::endl; 
  std::abort();
}

void SL_CONF::warning_fail (const char* expr, const char* file, 
		   int line, const char* msg) {
  std::cerr << "Warning: check violation!" << std::endl
	    << "Expr: " << expr << std::endl
	    << "File(Line): " << file << " (" << line << ")" << std::endl;
  if (msg != 0) std::cerr << "Explanation:" << msg << std::endl; 
}


// ======================================================================
// Math definitions and functions
// ======================================================================

#ifndef NULL
# define NULL (void *)0
#endif

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
# define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
# define M_PI_4 0.78539816339744830962
#endif

#ifndef M_LOG2E
# define M_LOG2E 1.4426950408889634074
#endif

///
template <class T>
inline T ABS(T a) { return (a >= 0 ? a : -a); }
///
template <class T>
inline T MIN(T a, T b) { return (a < b ? a : b); }
///
template <class T>
inline T MAX(T a, T b) { return (a > b ? a : b); }

///
template <class T>
inline void SWAP(T& a, T& b) {
  T t = a;
  a = b;
  b = t;
}

///
const float FLOAT_EPSILON = (float) 1.0e-7;
///
const double DOUBLE_EPSILON = 1.0e-16;

///
inline bool ALMOST_EQUAL(float a, float b) {
   return ABS(a-b) < FLOAT_EPSILON;
}
///
inline bool ALMOST_EQUAL(double a, double b) {
   return ABS(a-b) < DOUBLE_EPSILON;
}
///
template <class T, class S>
inline bool ALMOST_EQUAL(T a, S b, double eps) {
   return ABS(a-b) < eps;
}



#endif //_CONFIG_H_08192001
