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
// file          : test/BooleanOperations/include/CGAL/test_bops.h
// source        : test/BooleanOperations/include/CGAL/test_bops.h
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

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif
#ifndef CGAL_HOMOGENEOUS_H
#include <CGAL/Homogeneous.h>
#endif
#ifndef CGAL_BOOLEAN_OPERATIONS_2_H
#include <CGAL/boolean_operations_2.h>
#endif

#include <list>

#ifndef CGAL_Bops_test_arithmetic
#define CGAL_Bops_test_arithmetic 2
#endif

#if CGAL_Bops_test_arithmetic == 1
#ifndef CGAL_RATIONAL_H
#include <CGAL/Rational.h>
#endif
    typedef CGAL::Rational TestNum;            /* exact arithmetic */
#endif

#ifndef CGAL_QUOTIENT_H
#include <CGAL/Quotient.h>
#endif
#if CGAL_Bops_test_arithmetic == 2
    typedef CGAL::Quotient<long int> TestNum;  /* exact but (very) finite */
#endif

#if CGAL_Bops_test_arithmetic == 4
#ifndef CGAL_DOUBLE_H
#include <CGAL/Double.h>
#endif
    typedef double TestNum;             /* inexact arithmetic */
#endif

#ifdef CGAL_BOPS_HOMOGENEOUS
  typedef CGAL::Homogeneous<TestNum>  R_type;
#else
  typedef CGAL::Cartesian<TestNum>  R_type;
#endif

typedef CGAL::Point_2<R_type>           Point;
typedef CGAL::Segment_2<R_type>         Segment;
typedef CGAL::Triangle_2<R_type>        Triangle;
typedef CGAL::Iso_rectangle_2<R_type>   Iso_rectangle;
typedef std::list< Point >                 Container;
//typedef vector< Point >              Container;
typedef CGAL::Polygon_2< CGAL::Polygon_traits_2<R_type>, Container > Polygon;

#include "CGAL/test_bops_data.h"

#endif  // TEST_BOPS_H
