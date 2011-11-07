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

#include <CGAL/Handle.h>
#include <CGAL/Arr_point_location/Td_ninetuple.h>
#include <CGAL/Arr_point_location/Td_seventuple.h>
#include <CGAL/Arr_point_location/Td_fivetuple.h>

#ifndef CGAL_TD_DAG_H
#include <CGAL/Arr_point_location/Td_dag_node.h>
#endif

#define CGAL_TD_DELETE_SIGNATURE 0xffffffff

//type
#define CGAL_TD_VERTEX               0
#define CGAL_TD_EDGE                 0x1
#define CGAL_TD_TRAPEZOID            0x2
#define CGAL_TD_TYPE_MASK            0x3

#define CGAL_TD_ON_LEFT_BOUNDARY     0x4
#define CGAL_TD_ON_RIGHT_BOUNDARY    0x8
#define CGAL_TD_ON_BOTTOM_BOUNDARY   0x10
#define CGAL_TD_ON_TOP_BOUNDARY      0x20
#define CGAL_TD_ON_ALL_BOUNDARIES \
   (CGAL_TD_ON_LEFT_BOUNDARY | CGAL_TD_ON_RIGHT_BOUNDARY | \
   CGAL_TD_ON_BOTTOM_BOUNDARY| CGAL_TD_ON_TOP_BOUNDARY)
#define CGAL_TD_INTERIOR             0

#define CGAL_TD_CV_MIN_END           0
#define CGAL_TD_CV_MAX_END           0x1 

#define CGAL_TD_DEFAULT_DEPTH_THRESHOLD 60
#define CGAL_TD_DEFAULT_SIZE_THRESHOLD 12

#ifndef _MSC_VER
#if !defined __GNUC__ || __GNUC__>2 || __GNUC__==2 && __GNUC_MINOR__>=95
#define CGAL_PM_FRIEND_CLASS
#endif
#endif

#endif
