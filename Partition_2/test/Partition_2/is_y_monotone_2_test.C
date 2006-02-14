// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/is_y_monotone_2_test.C
// package       : $CGAL_Package: Partition_2 1.0 (27 Jul 2000) $
// chapter       : Planar Polygon Partitioning
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: testing of is_y_monotone_2 function
// ============================================================================

#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/is_y_monotone_2.h>
#include <list>
#include <vector>
#include <cassert>

typedef CGAL::Cartesian<double>          CR;
typedef CR::Point_2                      CPoint_2;
typedef std::list<CPoint_2>              CContainer;
typedef CGAL::Polygon_2<CR, CContainer>  CPolygon_2;

typedef CGAL::Homogeneous<double>        HR;
typedef HR::Point_2                      HPoint_2;
typedef std::vector<HPoint_2>            HContainer;
typedef CGAL::Polygon_2<HR, HContainer>  HPolygon_2;

template <class Polygon_2>
void make_nonmonotone_polygon(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(380, 485));
   polygon.push_back(Point_2(304, 335));
   polygon.push_back(Point_2(182, 405));
   polygon.push_back(Point_2(102, 319));
   polygon.push_back(Point_2( 72, 176));
   polygon.push_back(Point_2(204, 202));
   polygon.push_back(Point_2(286,  96));
   polygon.push_back(Point_2(404, 190));
   polygon.push_back(Point_2(428, 263));
   polygon.push_back(Point_2(213, 224));
   polygon.push_back(Point_2(400, 320));
}

template <class Polygon_2>
void make_monotone_polygon(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(266, 251));
   polygon.push_back(Point_2(131, 170));
   polygon.push_back(Point_2(186, 146));
   polygon.push_back(Point_2(266, 59));
   polygon.push_back(Point_2(277, 136));
   polygon.push_back(Point_2(260, 186));
   polygon.push_back(Point_2(370, 233));
   polygon.push_back(Point_2(399, 321));
   polygon.push_back(Point_2(318, 291));
}

int main(void)
{
   CPolygon_2 c_polygon;

   make_monotone_polygon(c_polygon);
   assert(CGAL::is_y_monotone_2(c_polygon.vertices_begin(), 
                                c_polygon.vertices_end()));

   c_polygon.erase(c_polygon.vertices_begin(), c_polygon.vertices_end());
   make_nonmonotone_polygon(c_polygon);
   assert(!CGAL::is_y_monotone_2(c_polygon.vertices_begin(), 
                                 c_polygon.vertices_end()));

   HPolygon_2 h_polygon;

   make_monotone_polygon(h_polygon);
   assert(CGAL::is_y_monotone_2(h_polygon.vertices_begin(), 
                                h_polygon.vertices_end()));

   h_polygon.erase(h_polygon.vertices_begin(), h_polygon.vertices_end());
   make_nonmonotone_polygon(h_polygon);
   assert(!CGAL::is_y_monotone_2(h_polygon.vertices_begin(), 
                                 h_polygon.vertices_end()));
   return 0;
}
