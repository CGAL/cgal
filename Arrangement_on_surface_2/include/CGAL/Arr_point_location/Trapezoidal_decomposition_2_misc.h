// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>
#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_MISC_H
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_MISC_H

#ifndef CGAL_TD_TENTUPLE_H
#include <CGAL/Handle.h>
#include <CGAL/Arr_point_location/Td_ninetuple.h>
#endif

#ifndef CGAL_TD_DAG_H
#include <CGAL/Arr_point_location/Td_dag.h>
#endif

#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_DELETE_SIGNATURE 0xffffffff
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED 0x1
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED 0x2
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED 0x4
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED 0x8
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOTALLY_UNBOUNDED \
	(CGAL_TRAPEZOIDAL_DECOMPOSITION_2_LEFT_UNBOUNDED| \
	CGAL_TRAPEZOIDAL_DECOMPOSITION_2_RIGHT_UNBOUNDED| \
	CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOTTOM_UNBOUNDED| \
	CGAL_TRAPEZOIDAL_DECOMPOSITION_2_TOP_UNBOUNDED) 
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_BOUNDED 0
#define CGAL_TD_DEFAULT_DEPTH_THRESHOLD 2
#define CGAL_TD_DEFAULT_SIZE_THRESHOLD 2

#ifndef _MSC_VER
#if !defined __GNUC__ || __GNUC__>2 || __GNUC__==2 && __GNUC_MINOR__>=95
#define CGAL_PM_FRIEND_CLASS
#endif
#endif

#endif
