#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;

typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef std::vector<Point_3>::iterator point_iterator;

int main() {

  std::vector<Point_3> points;

  points.push_back(Point_3(1,2,3));
  points.push_back(Point_3(1,2,5));

  Nef_polyhedron N1(points.begin(), points.end(), Nef_polyhedron::Points_tag());
  Nef_polyhedron N2(points.begin(), points.end(), Nef_polyhedron::Points_tag());

  Nef_polyhedron N3(Segment_3(Point_3(1,2,3),Point_3(1,2,5)));
  Nef_polyhedron N4(Segment_3(Point_3(1,2,3),Point_3(1,2,-5)));

  Nef_polyhedron N5(Point_3(1,2,3));

  // test non emptyness
  assert(!N1.is_empty());
  assert(!N2.is_empty());
  assert(!N3.is_empty());
  assert(!N4.is_empty());
  assert(!N5.is_empty());
  // test validitiy
  assert(N1.is_valid());
  assert(N2.is_valid());
  assert(N3.is_valid());
  assert(N4.is_valid());
  assert(N5.is_valid());
  // test points-points intersection
  N1*=N2;
  assert(N1.is_valid());
  assert(N1==N2);
  // test segment-segment intersection
  assert( N3*N4 == N5 );
  // test point-segment intersection
  assert( (N1*N3).is_valid() );
  assert( N1*N3 == N1 );
  // test points-points difference
  assert( (N1-N2).is_valid() );
  assert( (N1-N2).is_empty() );
  // test opening and closing segment
  assert( ((N3-N4)+N5).is_valid() );
  assert( ((N3-N4)+N5)==N3 );
  assert( ((N3-N4)*N5).is_empty() );

  return 0;
}
