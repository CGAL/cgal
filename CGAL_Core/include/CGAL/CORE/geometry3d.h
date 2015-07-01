/****************************************************************************
 * Core Library Version 1.7, August 2004                                     
 * Copyright (c) 1995-2004 Exact Computation Project                         
 * All rights reserved.                                                      
 *                                                                           
 * This file is part of CGAL (www.cgal.org).                
 * You can redistribute it and/or modify it under the terms of the GNU       
 * Lesser General Public License as published by the Free Software Foundation,      
 * either version 3 of the License, or (at your option) any later version.   
 *                                                                           
 * Licensees holding a valid commercial license may use this file in         
 * accordance with the commercial license agreement provided with the        
 * software.                                                                 
 *                                                                           
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE   
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. 
 *                                                                           
 *                                                                           
 * $URL$                                                                     
 * $Id$                                                                      
 ***************************************************************************/
/******************************************************************
 * Core Library Version 1.7, August 2004
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

#ifndef CORE_LEVEL
#  define CORE_LEVEL 3
#endif

#include <CGAL/CORE/linearAlgebra.h>

class Point3d;
class Line3d;
class Segment3d;
class Plane3d;
class Triangle3d;
class Polygon3d;

#include <CGAL/CORE/geom3d/point3d.h>
#include <CGAL/CORE/geom3d/line3d.h>
#include <CGAL/CORE/geom3d/segment3d.h>
#include <CGAL/CORE/geom3d/plane3d.h>
#include <CGAL/CORE/geom3d/triangle3d.h>
#include <CGAL/CORE/geom3d/polygon3d.h>

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
