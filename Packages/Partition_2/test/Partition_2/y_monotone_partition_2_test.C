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
// file          : include/CGAL/y_monotone_partition_2_test.C
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
// implementation: Testing of y-monotone partitioning function
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
void make_monotone_w_collinear_points(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

/*
   polygon.push_back(Point_2(186, 101));
   polygon.push_back(Point_2(304, 101));
   polygon.push_back(Point_2(351, 101));
   polygon.push_back(Point_2(367, 122));
   polygon.push_back(Point_2(388, 151));
   polygon.push_back(Point_2(407, 179));
   polygon.push_back(Point_2(339, 208));
   polygon.push_back(Point_2(280, 220));
   polygon.push_back(Point_2(214, 204));
   polygon.push_back(Point_2(214, 178));
   polygon.push_back(Point_2(214, 158));
   polygon.push_back(Point_2(164, 156));
*/   
/*
   polygon.push_back(Point_2(149,  91));
   polygon.push_back(Point_2(309,  91));
   polygon.push_back(Point_2(351,  91));
   polygon.push_back(Point_2(384,  91));
   polygon.push_back(Point_2(384, 127));
   polygon.push_back(Point_2(384, 158));
   polygon.push_back(Point_2(384, 186));
   polygon.push_back(Point_2(220, 186));
   polygon.push_back(Point_2(152, 186));
*/
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
void make_nonmonotone(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(336, 218));
   polygon.push_back(Point_2(444, 290));
   polygon.push_back(Point_2(416, 402));
   polygon.push_back(Point_2(282, 331));
   polygon.push_back(Point_2(366, 290));
   polygon.push_back(Point_2(314, 245));
   polygon.push_back(Point_2(209, 269));
   polygon.push_back(Point_2(202, 371));
   polygon.push_back(Point_2( 75, 261));
   polygon.push_back(Point_2(148, 267));
   polygon.push_back(Point_2(166, 231));
   polygon.push_back(Point_2(176, 256));
   polygon.push_back(Point_2(215, 152));
   polygon.push_back(Point_2(145, 134));
   polygon.push_back(Point_2(286, 123));
   polygon.push_back(Point_2(226, 197));
   polygon.push_back(Point_2(274, 208));
   polygon.push_back(Point_2(375, 125));
}

template <class Polygon_2>
void test_y_monotone(Polygon_2& polygon)
{
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_monotone_convex(polygon);
   CGAL::y_monotone_partition_2(polygon.vertices_begin(), 
                                polygon.vertices_end(),
                                std::back_inserter(partition_polys));

   assert (partition_polys.size() == 1 && 
           partition_polys.front().size() == polygon.size());
   assert(CGAL::is_y_monotone_2(partition_polys.front().vertices_begin(), 
                                partition_polys.front().vertices_end()));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_monotone_w_collinear_points(polygon);
   CGAL::y_monotone_partition_2(polygon.vertices_begin(), 
                                polygon.vertices_end(),
                                std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonmonotone(polygon);
   CGAL::y_monotone_partition_2(polygon.vertices_begin(), 
                                polygon.vertices_end(),
                                std::back_inserter(partition_polys));
}

int main(void)
{
   CPolygon_2 c_polygon;
   HPolygon_2 h_polygon;
   test_y_monotone(c_polygon);
   test_y_monotone(h_polygon);

   return 0;
}
