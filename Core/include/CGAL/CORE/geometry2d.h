/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: geometry2d.h
 * Synopsis:
 *      Basic 2-dimensional geometry
 *
 * Written by
 *      Yaping Yuan (yqy0522@cs.nyu.edu), 1999.
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_GEOMETRY2D_H
#define CORE_GEOMETRY2D_H

#ifndef CORE_LEVEL
#  define CORE_LEVEL 3
#endif

#include <CGAL/CORE/geom2d/point2d.h>
#include "CGAL/CORE/geom2d/line2d.h"
#include "CGAL/CORE/geom2d/circle2d.h"
#include "CGAL/CORE/geom2d/segment2d.h"

// automaticall link necessary static library under visual c++
#ifdef _MSC_VER
	#if CORE_LEVEL == 1
		#ifdef _DEBUG
			#pragma comment(lib, "corexDebug_level1.lib")
		#else
			#pragma comment(lib, "corex_level1.lib")
		#endif
	#elif CORE_LEVEL == 2
		#ifdef _DEBUG
			#pragma comment(lib, "corexDebug_level2.lib")
		#else
			#pragma comment(lib, "corex_level2.lib")
		#endif
	#elif CORE_LEVEL == 3
		#ifdef _DEBUG
			#pragma comment(lib, "corexDebug_level3.lib")
		#else
			#pragma comment(lib, "corex_level3.lib")
		#endif
	#endif
#endif

#endif
