// Copyright (c) 1997-2000  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Thomas Herrmann

#ifndef CGAL_WIDTH_ASSERTIONS_H
#define CGAL_WIDTH_ASSERTIONS_H 1

#include <CGAL/license/Polytope_distance_d.h>


#ifdef SIMPLIFY
#define GCD_COMPUTATION 1
#endif

#ifdef DEBUG

//Turn assertion output on/off 
#define ASSERTION_OUTPUT 0
#define EXPENSIVE_CHECKS_OUTPUT 0

//Turn on/off output in preparation_check
#define PREPARATION_CHECK 0

//Turn on/off output in neighbors_of
#define NEIGHBORS_OF 0

//Turn on/off output of setminus, setunion and setcut
#define SETMINUS 0
#define SETUNION 0
#define SETCUT 0

//Turn on/off output of compute_plane_equation
#define COMPUTE_PLANE_EQUATION 0

//Turn on/off output of solve_3x3
#define SOLVE_3X3 0

//Turn on/off output of solve_4x4
#define SOLVE_4X4 0

//Turn on/off output in check_feasibility
#define CHECK_FEASIBILITY 0

//Turn on/off compilation and output of gcd computation
#define GCD_OUTPUT 0

//Turn on/off output of simplify_solution
#define SIMPLIFY_SOLUTION 0

//Turn on/off output in initial_VF_pair
#define INITIAL_VF_PAIR 0

//Turn on/off output in check_about_VF-pairs
#define CHECK_ABOUT_VF_PAIRS 0
#define VF_PAIR_OUTPUT 0

//Turn on/off output of update_width
#define UPDATE_WIDTH 0

//Turn on/off output of EE_computation and EE_pairs
#define EE_COMPUTATION 0
#define EE_PAIRS 0

//Turn on/off output in origin_inside_CH
#define ORIGIN_INSIDE_CH 0

//Turn on/off output of width_3_convex
#define WIDTH_3_CONVEX 0
#define EDGE_INITIALIZING 0

//Turn on/off output of stack go_on
#define GO_ON_OUTPUT 0

//Turn infos on/off
#define INFO 0

//Turn on/off output of verifications on edges
#define VISITED_CHECK 0
#define IMPASSABLE_CHECK 0

 #include<stream.h>

 #define DEBUGENDL(doit,msg,var)\
 if(doit!=0) std::cout << msg << " " << var << endl;

 #define DEBUGPRINT(doit,msg,var)\
 if(doit!=0) std::cout << msg << " " << var;

 #define DEBUGMSG(doit,msg)\
 if(doit!=0) std::cout << msg << endl;

 #define INFOMSG(doit,msg)\
 if(doit!=0) std::cerr<<msg<<endl;

#else
 #define DEBUGENDL(doit,msg,var)
 #define DEBUGPRINT(doit,msg,var)
 #define DEBUGMSG(doit,msg)
 #define INFOMSG(doit,msg)
#endif

#endif //#WIDTH_DEBUG_H
