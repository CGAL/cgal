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
// file          : include/CGAL/convex_partiton_2_test.C
// package       : $CGAL_Package: Partition_2 1.0 (27 Jul 2000) $
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Testing of convex partitioning functions
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_traits_2.h>
#include <CGAL/partition_2.h>
#include <list>
#include <cassert>

typedef CGAL::Cartesian<double>               CR;
typedef CR::Point_2                           CPoint_2;
typedef CGAL::Polygon_traits_2<CR>            CTraits;
typedef std::list<CPoint_2>                   CContainer;
typedef CGAL::Polygon_2<CTraits, CContainer>  CPolygon_2;

typedef CGAL::Homogeneous<double>             HR;
typedef HR::Point_2                           HPoint_2;
typedef CGAL::Polygon_traits_2<HR>            HTraits;
typedef std::list<HPoint_2>                   HContainer;
typedef CGAL::Polygon_2<HTraits, HContainer>  HPolygon_2;

template <class Polygon_2>
void make_convex_w_collinear_points(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(140, 226));
   polygon.push_back(Point_2(293, 226));
   polygon.push_back(Point_2(335, 226));
   polygon.push_back(Point_2(358, 226));
   polygon.push_back(Point_2(358, 264));
   polygon.push_back(Point_2(358, 286));
   polygon.push_back(Point_2(358, 315));
   polygon.push_back(Point_2(358, 334));
   polygon.push_back(Point_2(297, 334));
   polygon.push_back(Point_2(246, 334));
   polygon.push_back(Point_2(211, 334));
   polygon.push_back(Point_2(163, 334));
   polygon.push_back(Point_2(146, 307));
   polygon.push_back(Point_2(132, 278));
   polygon.push_back(Point_2(132, 259));
   polygon.push_back(Point_2(132, 249));
   polygon.push_back(Point_2(132, 236));

}

template <class Polygon_2>
void make_monotone_convex(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(239, 108));
   polygon.push_back(Point_2(358, 170));
   polygon.push_back(Point_2(387, 305));
   polygon.push_back(Point_2(340, 416));
   polygon.push_back(Point_2(215, 432));
   polygon.push_back(Point_2(106, 358));
   polygon.push_back(Point_2(106, 236));
   polygon.push_back(Point_2(138, 162));
}

template <class Polygon_2>
void make_nonconvex_w_collinear_points(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(129, 147));
   polygon.push_back(Point_2(232, 147));
   polygon.push_back(Point_2(344, 147));
   polygon.push_back(Point_2(365, 234));
   polygon.push_back(Point_2(239, 251));
   polygon.push_back(Point_2(249, 210));
   polygon.push_back(Point_2(190, 212));
   polygon.push_back(Point_2(205, 323));
   polygon.push_back(Point_2(104, 323));
   polygon.push_back(Point_2(104, 279));
   polygon.push_back(Point_2(104, 243));
}

template <class Polygon_2>
void make_nonconvex(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(391, 374));
   polygon.push_back(Point_2(240, 431));
   polygon.push_back(Point_2(252, 340));
   polygon.push_back(Point_2(374, 320));
   polygon.push_back(Point_2(289, 214));
   polygon.push_back(Point_2(134, 390));
   polygon.push_back(Point_2( 68, 186));
   polygon.push_back(Point_2(154, 259));
   polygon.push_back(Point_2(161, 107));
   polygon.push_back(Point_2(435, 108));
   polygon.push_back(Point_2(208, 148));
   polygon.push_back(Point_2(295, 160));
   polygon.push_back(Point_2(421, 212));
   polygon.push_back(Point_2(441, 303));
}

template <class Polygon_2>
void test_optimal_convex(Polygon_2& polygon)
{
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_monotone_convex(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));

   assert (partition_polys.size() == 1 && 
           partition_polys.front().size() == polygon.size());
   assert(CGAL::is_convex_2(partition_polys.front().vertices_begin(), 
                            partition_polys.front().vertices_end()));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_convex_w_collinear_points(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex_w_collinear_points(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                   polygon.vertices_end(),
                                   std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
}

template <class Polygon_2>
void test_one_optimal_convex(Polygon_2& polygon)
{
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
}

template <class Polygon_2>
void test_approx_convex(Polygon_2& polygon)
{
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_monotone_convex(polygon);
   CGAL::approx_convex_partition_2(polygon.vertices_begin(), 
                                   polygon.vertices_end(),
                                   std::back_inserter(partition_polys));

   assert (partition_polys.size() == 1 &&
           partition_polys.front().size() == polygon.size());
   assert(CGAL::is_convex_2(partition_polys.front().vertices_begin(),
                            partition_polys.front().vertices_end()));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_convex_w_collinear_points(polygon);
   CGAL::approx_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex_w_collinear_points(polygon);
   CGAL::approx_convex_partition_2(polygon.vertices_begin(),
                                   polygon.vertices_end(),
                                   std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex(polygon);
   CGAL::approx_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
}


template <class Polygon_2>
void test_one_approx_convex(Polygon_2& polygon)
{
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex(polygon);
   CGAL::approx_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
}

template <class Polygon_2>
void test_greene_approx_convex(Polygon_2& polygon)
{
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_monotone_convex(polygon);
   CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
                                          polygon.vertices_end(),
                                          std::back_inserter(partition_polys));

   assert (partition_polys.size() == 1 &&
           partition_polys.front().size() == polygon.size());
   assert(CGAL::is_convex_2(partition_polys.front().vertices_begin(),
                            partition_polys.front().vertices_end()));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_convex_w_collinear_points(polygon);
   CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
                                          polygon.vertices_end(),
                                          std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex_w_collinear_points(polygon);
   CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
                                          polygon.vertices_end(),
                                          std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex(polygon);
   CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
                                          polygon.vertices_end(),
                                          std::back_inserter(partition_polys));
}

template <class Polygon_2>
void test_one_greene_approx_convex(Polygon_2& polygon)
{
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex(polygon);
   CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
                                          polygon.vertices_end(),
                                          std::back_inserter(partition_polys));
}

int main(void)
{
   CPolygon_2 c_polygon;
   HPolygon_2 h_polygon;
   test_optimal_convex(c_polygon);
   test_one_optimal_convex(h_polygon);
   test_approx_convex(c_polygon);
   test_one_approx_convex(h_polygon);
   test_greene_approx_convex(c_polygon);
   test_one_greene_approx_convex(h_polygon);

   return 0;
}
