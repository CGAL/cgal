/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: geometry3d.h
 *
 * Written by
 *       Shubin Zhao (shubinz@cs.nyu.edu) (2001)
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_GEOMETRY3D_H
#define CORE_GEOMETRY3D_H

#ifndef Level
#  define Level 3
#endif

#include "linearAlgebra.h"

class Point3d;
class Line3d;
class Segment3d;
class Plane3d;
class Triangle3d;
class Polygon3d;

#include "geom3d/point3d.h"
#include "geom3d/line3d.h"
#include "geom3d/segment3d.h"
#include "geom3d/plane3d.h"
#include "geom3d/triangle3d.h"
#include "geom3d/polygon3d.h"

// automaticall link necessary static library under visual c++
#ifdef _MSC_VER
	#if Level == 1
		#ifdef _DEBUG
			#pragma comment(lib, "corexDebug_level1.lib")
		#else
			#pragma comment(lib, "corex_level1.lib")
		#endif
	#elif Level == 2
		#ifdef _DEBUG
			#pragma comment(lib, "corexDebug_level2.lib")
		#else
			#pragma comment(lib, "corex_level2.lib")
		#endif
	#elif Level == 3
		#ifdef _DEBUG
			#pragma comment(lib, "corexDebug_level3.lib")
		#else
			#pragma comment(lib, "corex_level3.lib")
		#endif
	#endif
#endif

#endif
