//  -*- Mode: c++ -*-
// ============================================================================
// 
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-1.0 $
// release_date  : $CGAL_Date: 1998/09/12 $
//
// file          : demo/BooleanOperations/include/CGAL/test_bops.h
// source        : demo/BooleanOperations/include/CGAL/test_bops.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#ifndef TEST_BOPS_H
#define TEST_BOPS_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <list>

//#define CGAL__BOPS_DEBUG_ON
//#define CGAL__DCEL_DEBUG_ON
//#define _DCEL__V2E_DEBUG_ON
//#define CGAL__INTERSECTING_POLYGONS_DEBUG_ON

#include <CGAL/boolean_operations_2.h>
#include <iostream>
//#include <CGAL/IO/ostream.h>

#ifndef Bops_test_arithmetic
#define Bops_test_arithmetic 1
#endif

#if Bops_test_arithmetic == 1
#   include <CGAL/leda_rational.h>
    typedef leda_rational TestNum;            /* exact arithmetic */
#endif

#if Bops_test_arithmetic == 2
#   include <CGAL/Quotient.h>
    typedef Quotient<long int> TestNum;  /* exact but (very) finite */
#endif

#if Bops_test_arithmetic == 4
//#   include <CGAL/Double.h>
#   include <CGAL/double.h>
    typedef double TestNum;             /* inexact arithmetic */
#endif

#ifdef BOPS_HOMOGENEOUS
#   include <CGAL/Quotient.h>
  typedef Homogeneous<TestNum>  R_type;
#else
  typedef Cartesian<TestNum>  R_type;
#endif

typedef Point_2<R_type>           Point_2;
typedef Segment_2<R_type>         Segment_2;
typedef Triangle_2<R_type>        Triangle_2;
typedef Iso_rectangle_2<R_type>   Iso_rectangle_2;
typedef list< Point_2 >                Container;
//typedef vector< Point_2 >              Container;
typedef Polygon_2< Polygon_traits_2<R_type>, Container >  Polygon_2;

#include "CGAL/test_bops_io.h"

#endif  // TEST_BOPS_H
