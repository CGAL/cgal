#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/point_generators_3.h>

#include <vector>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>


template <class K, class Polyhedron_3>
void test()
{
  typedef typename K::Plane_3                                            Plane;
  typedef typename K::Point_3                                            Point;

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

  assert(num_vertices(P1)==8 && num_faces(P1)==6);
  assert(num_vertices(P2)==8 && num_faces(P2)==6);
  assert(num_vertices(P3)==8 && num_faces(P3)==6);
  assert(num_vertices(P4)==8 && num_faces(P4)==6);
  assert(num_vertices(P5)==8 && num_faces(P5)==6);

  using CGAL::Convex_hull_3::internal::point_inside_convex_polyhedron;
  assert(point_inside_convex_polyhedron(P1, Point(0,0,0)));
  assert(point_inside_convex_polyhedron(P2, Point(0,0,0)));
  assert(point_inside_convex_polyhedron(P3, Point(0,0,0)));
  assert(point_inside_convex_polyhedron(P4, Point(0,0,0)));
  assert(point_inside_convex_polyhedron(P5, Point(0,0,0)));
}


int main()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
  typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;

  test<Epic,CGAL::Polyhedron_3<Epic> >();
  test<Epec,CGAL::Polyhedron_3<Epec> >();
  test<Epic,CGAL::Surface_mesh<Epic::Point_3> >();
  test<Epec,CGAL::Surface_mesh<Epec::Point_3> >();
}
