#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/point_generators_3.h>

#include <vector>
#include <CGAL/IO/Polyhedron_iostream.h>

template <class K>
void test()
{
  typedef typename K::Plane_3                                            Plane;
  typedef typename K::Point_3                                            Point;
  typedef CGAL::Polyhedron_3<K>                                 Polyhedron_3;

  // generates supporting planes of the facets of a cube
  std::vector<Plane> planes;

  planes.push_back( Plane( 1, 0, 0,-1) ); // x= 1
  planes.push_back( Plane( 0, 1, 0,-1) ); // y= 1
  planes.push_back( Plane( 0,-1, 0,-1) ); // y=-1
  planes.push_back( Plane( 0, 0, 1,-1) ); // z= 1
  planes.push_back( Plane( 0, 0,-1,-1) ); // z=-1
  planes.push_back( Plane(-1, 0, 0,-1) ); // x=-1

  // define polyhedron to hold the intersection
  Polyhedron_3 P1, P2, P3, P4, P5;

  // test halfspace_intersection_3 with a point inside
  CGAL::halfspace_intersection_3(planes.begin(),
                                 planes.end(),
                                 P1,
                                 boost::make_optional(Point(0, 0, 0)) );

  // test halfspace_intersection_3 with non point inside
  CGAL::halfspace_intersection_3(planes.begin(),
                                 planes.end(),
                                 P2);

  // test halfspace_intersection_with_constructions_3 with a point inside
  CGAL::halfspace_intersection_with_constructions_3( planes.begin(),
                                                     planes.end(),
                                                     P3,
                                                     boost::make_optional(Point(0, 0, 0)) );

  // test halfspace_intersection_with_constructions_3 with non point inside
  CGAL::halfspace_intersection_with_constructions_3( planes.begin(),
                                                     planes.end(),
                                                     P4);

  // test halfspace_intersection_with_constructions_3 with non point inside but with the kernel
  CGAL::halfspace_intersection_with_constructions_3( planes.begin(),
                                                     planes.end(),
                                                     P5,
                                                     boost::optional<Point>(),
                                                     K());

  assert(P1.size_of_vertices()==8 && P1.size_of_facets()==6);
  assert(P2.size_of_vertices()==8 && P2.size_of_facets()==6);
  assert(P3.size_of_vertices()==8 && P3.size_of_facets()==6);
  assert(P4.size_of_vertices()==8 && P4.size_of_facets()==6);
  assert(P5.size_of_vertices()==8 && P5.size_of_facets()==6);

  using CGAL::Convex_hull_3::internal::point_inside_convex_polyhedron;
  assert(point_inside_convex_polyhedron(P1, Point(0,0,0)));
  assert(point_inside_convex_polyhedron(P2, Point(0,0,0)));
  assert(point_inside_convex_polyhedron(P3, Point(0,0,0)));
  assert(point_inside_convex_polyhedron(P4, Point(0,0,0)));
  assert(point_inside_convex_polyhedron(P5, Point(0,0,0)));
}


int main()
{
  test<CGAL::Exact_predicates_inexact_constructions_kernel>();
  test<CGAL::Exact_predicates_exact_constructions_kernel>();
}
