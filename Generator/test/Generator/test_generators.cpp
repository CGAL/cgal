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
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : test_generators.C
// chapter       : $CGAL_Chapter: Geometric Object Generators $
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : generators.fw
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//                 Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Test Random Sources and Geometric Object Generators
// ============================================================================


#include <cstddef>
#include <vector>
#include <algorithm>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/random_selection.h>

using namespace CGAL;

typedef Cartesian<double>                R;
typedef CGAL::Cartesian_d<double>        Kd;




void test_point_generators_2() {
    typedef Point_2<R>                       Point_2;
    typedef Creator_uniform_2<double,Point_2>  Creator;

    // Create test point set.
    std::vector<Point_2> points;
    points.reserve(1000);
    Random_points_in_disc_2<Point_2,Creator>     g1( 100.0);
    CGAL::cpp11::copy_n( g1, 100, std::back_inserter(points));
    Random_points_on_circle_2<Point_2,Creator>   g2( 100.0);
    Random_points_in_square_2<Point_2,Creator>   g3( 100.0);
    Random_points_on_square_2<Point_2,Creator>   g4( 100.0);
    Random_points_on_segment_2<Point_2,Creator>  g5( Point_2(-100,-100),
                                                    Point_2( 100, 100));
    Points_on_segment_2<Point_2>                g5a( Point_2( 50,-50),
                                                    Point_2(-50, 50),
                                                   50);
    CGAL::cpp11::copy_n( g2, 100, std::back_inserter(points));
    CGAL::cpp11::copy_n( g3, 100, std::back_inserter(points));
    CGAL::cpp11::copy_n( g4, 100, std::back_inserter(points));
    CGAL::cpp11::copy_n( g5,  50, std::back_inserter(points));
    CGAL::cpp11::copy_n( g5a, 50, std::back_inserter(points));
    points_on_square_grid_2( 50.0, (std::size_t)1,
                             std::back_inserter(points), Creator());
    points_on_square_grid_2( 50.0, (std::size_t)2,
                             std::back_inserter(points), Creator());
    points_on_square_grid_2( 50.0, (std::size_t)3,
                             std::back_inserter(points), Creator());
    points_on_square_grid_2( 50.0, (std::size_t)94,
                             std::back_inserter(points), Creator());
    points_on_segment_2( Point_2(-100, 100), Point_2( 100,-100),
                         (std::size_t)100, std::back_inserter(points));
    random_selection( points.begin(), points.end(), (std::size_t)100,
                      std::back_inserter(points));
    random_collinear_points_2( points.begin(), points.end(),
                               (std::size_t)100,
                               std::back_inserter(points));

    // Check perturbation. Make sure that the result stays within
    // the 100 x 100 square. 10 pixel perturbation allowed.
    Random_points_in_square_2<Point_2,Creator>   g6( 90.0);
    int count = 100 ;
    CGAL::cpp11::copy_n( g6, count, std::back_inserter(points));
    std::vector<Point_2>::iterator i2 = points.end();
    std::vector<Point_2>::iterator i1 = i2 ;
    std::advance(i1,-count);
    perturb_points_2( i1, i2, 10.0);

    // Create a random permutation.
    std::random_shuffle( points.begin(), points.end(), get_default_random());

    assert( points.size() == 1000);
    for ( std::vector<Point_2>::iterator i = points.begin();
          i != points.end(); i++){
        assert( i->x() <=  100);
        assert( i->x() >= -100);
        assert( i->y() <=  100);
        assert( i->y() >= -100);
    }
}

void test_point_generators_3() {
    typedef Point_3<R>                       Point_3;
    typedef Creator_uniform_3<double,Point_3>  Creator;

    /* Create test point set. */
    std::vector<Point_3> points;
    points.reserve(500);
    Random_points_in_sphere_3<Point_3,Creator>      g1( 100.0);
    CGAL::cpp11::copy_n( g1, 100, std::back_inserter(points));
    Random_points_on_sphere_3<Point_3,Creator>      g2( 100.0);
    Random_points_in_cube_3<Point_3,Creator>        g3( 100.0);
    CGAL::cpp11::copy_n( g2, 100, std::back_inserter(points));
    CGAL::cpp11::copy_n( g3, 100, std::back_inserter(points));
    points_on_cube_grid_3( 50.0, (std::size_t)1,
                           std::back_inserter(points), Creator());
    points_on_cube_grid_3( 50.0, (std::size_t)2,
                           std::back_inserter(points), Creator());
    points_on_cube_grid_3( 50.0, (std::size_t)3,
                           std::back_inserter(points), Creator());
    points_on_cube_grid_3( 50.0, (std::size_t)94,
                           std::back_inserter(points), Creator());
    random_selection( points.begin(), points.end(), 100,
                      std::back_inserter(points));

    assert( points.size() == 500);
    for ( std::vector<Point_3>::iterator i = points.begin();
          i != points.end(); i++){
        assert( i->x() <=  100);
        assert( i->x() >= -100);
        assert( i->y() <=  100);
        assert( i->y() >= -100);
        assert( i->z() <=  100);
        assert( i->z() >= -100);
    }
}

void test_point_generators_d()
{
 typedef Kd::Point_d Point;
 typedef CGAL::Creator_uniform_d<std::vector<double>::iterator, Point>Creator_d;

    std::vector<Point> points;
    int nb_g=10000;
    points.reserve (nb_g+403);
    int i=0,ii;
    
    {
    // 100 random points in dim 36
      std::cout<<"    cube dim 36"<<std::flush;
      CGAL::Random_points_in_cube_d<Point> gen (36, 1.0);
      CGAL::cpp11::copy_n( gen, 100, std::back_inserter(points));
      i+=100;
      std::cout<<" done"<<std::endl;
    }
    double o[26]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    {
    // 100 random points in dim 4
      std::cout<<"    in ball 4D"<<std::flush;
      Point o4(4,o,o+4);
      CGAL::Random_points_in_ball_d<Point> gen (4, 100.0);
      CGAL::cpp11::copy_n( gen, 100, std::back_inserter(points));
      std::cout<<" done"<<std::flush;
      for(ii=i+100; i<ii; ++i)
	assert( CGAL::squared_distance(o4,points[i]) <= 10000.0);
      std::cout<<" checked"<<std::endl;
    }
    {
      // nb_g random points in dim 3
      std::cout<<"    in ball 3D"<<std::flush;
      Point o3(3,o,o+3);
      Point g=o3;
      CGAL::Random_points_in_ball_d<Point> gen (3, 1.0);
      CGAL::cpp11::copy_n( gen, nb_g, std::back_inserter(points));
      std::cout<<" done"<<std::flush;
      for(ii=i+nb_g; i<ii; ++i){
	assert( CGAL::squared_distance(o3,points[i]) <= 1.0);
	if (points[i][0] >0) 
	  g = g + (points[i] - o3);
      }
      assert( std::fabs( g[0]/nb_g - 3.0/16.0) < 0.01 );
      std::cout<<" center of mass 3/16~="<<g[0]/nb_g<<" checked"<<std::endl;
    }
    {
      // 100 random points in dim 26
      std::cout<<"    on sphere 26D"<<std::flush;
      Point o26(26,o,o+26);
      CGAL::Random_points_on_sphere_d<Point> gen (26, 1.0);
      CGAL::cpp11::copy_n( gen, 100, std::back_inserter(points));
      std::cout<<" done"<<std::flush;
      for(ii=i+100; i<ii; ++i) {
	assert( CGAL::squared_distance(o26,points[i]) - 1.0 <= 0.1);
      }
      std::cout<<" checked"<<std::endl;
    }
    
    {
      std::cout<<"    on grid "<<std::flush;
      int dim =4;
      // 100 grid points in dim 4
      CGAL::points_on_cube_grid_d (dim, 1.0, (std::size_t) 100, 
				   std::back_inserter(points), Creator_d(dim) );
      // 1 grid points in dim 4
      CGAL::points_on_cube_grid_d (dim, 1.0, (std::size_t) 1, 
				   std::back_inserter(points), Creator_d(dim) );
      // 2 grid points in dim 4
      CGAL::points_on_cube_grid_d (dim, 1.0, (std::size_t) 2, 
				   std::back_inserter(points), Creator_d(dim) );
      i=i+103;
      std::cout<<" done"<<std::endl;
    }
}


int main(){
    std::cout<<"testing 2D"<<std::flush;
    test_point_generators_2();
    std::cout<<" done"<<std::endl;

    std::cout<<"testing 3D"<<std::flush;
    test_point_generators_3();
    std::cout<<" done"<<std::endl;

    std::cout<<"testing high D"<<std::endl;
    test_point_generators_d();
    std::cout<<" done"<<std::endl;
    return 0;
}
// EOF //
