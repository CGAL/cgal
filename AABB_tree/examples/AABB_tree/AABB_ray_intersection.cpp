#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/point_generators_3.h>

#include <algorithm>
#include <iterator>
#include <boost/functional/value_factory.hpp>
#include <boost/timer/timer.hpp>


typedef CGAL::Epeck K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;

int main()
{
  Point p(1.0, 0.0, 0.0);
  Point q(0.0, 1.0, 0.0);
  Point r(0.0, 0.0, 1.0);
  Point s(0.0, 0.0, 0.0);
  Polyhedron polyhedron;
  polyhedron.make_tetrahedron(p, q, r, s);

  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);


  const int NB_RAYS = 1000;
  std::vector<Point> v1, v2;
  v1.reserve(NB_RAYS); v2.reserve(NB_RAYS);

  const float rad = 2.0;
  // Generate NB_RAYS*2 points that lie on a sphere of radius r.
  CGAL::Random rand = CGAL::Random(23); // fix the seed to yield the same results each run
  CGAL::cpp11::copy_n(CGAL::Random_points_on_sphere_3<Point>(rad, rand), NB_RAYS, std::back_inserter(v1));
  CGAL::cpp11::copy_n(CGAL::Random_points_on_sphere_3<Point>(rad, rand), NB_RAYS, std::back_inserter(v2));

  // Generate NB_RAYS using v1 as source and v2 as target.
  std::vector<Ray> rays;
  rays.reserve(NB_RAYS);
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 std::back_inserter(rays), boost::value_factory<Ray>());

  // Calculate intersections between rays and primitives.
  std::vector< boost::optional<Tree::template Intersection_and_primitive_id<Ray>::Type > >
    first_intersections;
  for(std::vector<Ray>::iterator it = rays.begin(); it != rays.end(); ++it)
    first_intersections.push_back(tree.first_intersection(*it));

  //alternative if we are only intersected by the primitives
  std::vector< boost::optional<Tree::Primitive_id> > intersected_primitives;
  for(std::vector<Ray>::iterator it = rays.begin(); it != rays.end(); ++it) {
    intersected_primitives.push_back(tree.first_intersected_primitive(*it));
  }

  return 0;
}
