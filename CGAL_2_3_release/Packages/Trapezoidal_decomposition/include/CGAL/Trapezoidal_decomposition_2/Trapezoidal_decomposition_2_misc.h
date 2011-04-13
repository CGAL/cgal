// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Trapezoidal_decomposition_2/Trapezoidal_decomposition_2_misc.h
// package       : Trapezoidal decomposition 2
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_TRAPEZOIDAL_DECOMPOSITION_2_MISC_H
#define CGAL_TRAPEZOIDAL_DECOMPOSITION_2_MISC_H

#ifndef CGAL_TD_TENTUPLE_H
#include <CGAL/Trapezoidal_decomposition_2/Td_ninetuple.h>
#endif
#ifndef CGAL_TWOTUPLE_H
#include <CGAL/Twotuple.h>
#endif

#ifndef CGAL_TD_DAG_H
#include <CGAL/Td_dag.h>
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
#ifndef __BORLANDC__
#if !defined __GNUC__ || __GNUC__>2 || __GNUC__==2 && __GNUC_MINOR__>=95
#define CGAL_PM_FRIEND_CLASS
#endif
#endif
#endif

#endif



















