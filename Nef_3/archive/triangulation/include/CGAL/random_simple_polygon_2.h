#ifndef RANDOM_SIMPLE_POLYGON_2_H
#define RANDOM_SIMPLE_POLYGON_2_H

#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random_polygon_2_sweep.h>
#include <CGAL/function_objects.h>

namespace CGAL {

template <class Integer, class Real, class Polygon_2>
void random_simple_polygon_2( Integer num_verts, Real radius, Polygon_2& P) {

  typedef typename Polygon_2::Traits Traits;
  typedef typename Polygon_2::Point_2 Point_2;
  typedef Creator_uniform_2< double, Point_2>  Creator;

  Random_points_in_square_2< Point_2, Creator> pg(radius);
  random_polygon_2( num_verts, std::back_inserter(P), pg);
  make_simple_polygon( P.vertices_begin(), P.vertices_end(), Traits());
  CGAL_assertion( P.is_simple());
  CGAL_assertion( P.orientation() == COUNTERCLOCKWISE);
}

}

#endif // RANDOM_SIMPLE_POLYGON_2_H
