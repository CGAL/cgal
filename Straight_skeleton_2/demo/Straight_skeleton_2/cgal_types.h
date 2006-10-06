// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_TYPES
#define CGAL_TYPES

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <boost/shared_ptr.hpp>
#include <vector>

#if USE_BIG_FLOAT

#include <CORE/BigFloat.h>

namespace CGAL
{

template<> struct Number_type_traits<CORE::BigFloat>
{
  typedef Tag_false Has_exact_ring_operations;
  typedef Tag_false Has_exact_division;
  typedef Tag_false Has_exact_sqrt;

  typedef Tag_true Has_gcd;
  typedef Tag_true Has_division;
  typedef Tag_true Has_sqrt;
};

CORE::BigFloat sqrt( CORE::BigFloat const& n ) { return CORE::sqrt(n); }

}

#endif

namespace demo
{

//typedef CGAL::Simple_cartesian<double> K ;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt K;
//typedef CGAL::Simple_cartesian<CORE::BigFloat> K ;

typedef CGAL::Polygon_2<K>            CGAL_Polygon ;

typedef K::FT                         FT ;
typedef K::Point_2                    Point;
typedef std::vector<Point>            Polygon;
typedef boost::shared_ptr<Polygon>    PolygonPtr;
typedef CGAL::Segment_2<K>            Segment;
typedef std::vector<PolygonPtr>       Region ;
typedef boost::shared_ptr<Region>     RegionPtr ;
typedef std::vector<RegionPtr>        Regions ;
typedef std::set<double>              Doubles ;

}
#endif

