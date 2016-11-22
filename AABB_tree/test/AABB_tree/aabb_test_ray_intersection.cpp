#include <fstream>
#include <algorithm>
#include <iterator>

#include <boost/functional/value_factory.hpp>
#include <boost/array.hpp>

#include <CGAL/assertions.h>
#include <CGAL/algorithm.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Timer.h>

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
typedef CGAL::Timer Timer;

FT point_on_ray_dist(const Ray& ray, const Point& point) {
  Vector i_ray(point, ray.source());
  return i_ray.squared_length();
}

std::size_t accum = 0;

boost::optional<
  Tree::Intersection_and_primitive_id<Ray>::Type
  >
min_intersection(const Tree& tree, const Ray& ray) {
  typedef std::vector< Tree::Intersection_and_primitive_id<Ray>::Type > IntersectionVector;
  IntersectionVector all_intersections;

  tree.all_intersections(ray, std::back_inserter(all_intersections));
  accum += all_intersections.size();
  Tree::FT min_distance = DBL_MAX;
  boost::optional<
    Tree::Intersection_and_primitive_id<Ray>::Type
    > mini = boost::none;

  for(IntersectionVector::iterator it2 = all_intersections.begin(); it2 != all_intersections.end(); ++it2) {
    if(Point* point = boost::get<Point>(&(it2->first))) {
      Vector i_ray(*point, ray.source());
      Tree::FT new_distance = i_ray.squared_length();
      if(new_distance < min_distance) {
        mini = *it2;
        min_distance = new_distance;
      }
    } else {
      std::cout << "ERROR ignored a segment" << std::endl;
    }
  }
  return mini;
}

int main()
{
  Polyhedron polyhedron;
  {
    // Point p(1.0, 0.0, 0.0);
    // Point q(0.0, 1.0, 0.0);
    // Point r(0.0, 0.0, 1.0);
    // Point s(0.0, 0.0, 0.0);
    // polyhedron.make_tetrahedron(p, q, r, s);
  }

  std::ifstream in("data/bunny00.off");
  if(in)
    in >> polyhedron;
  else{
    std::cout << "error reading bunny" << std::endl;
    return 1;
  }

  Timer t;
  t.start();

  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
  Tree::Bounding_box bbox = tree.bbox();
  Vector bbox_center((bbox.xmin() + bbox.xmax()) / 2, 
                     (bbox.ymin() + bbox.ymax()) / 2, 
                     (bbox.zmin() + bbox.zmax()) / 2);
  boost::array<double, 3> extents;
  extents[0] = bbox.xmax() - bbox.xmin();
  extents[1] = bbox.ymax() - bbox.ymin();
  extents[2] = bbox.zmax() - bbox.zmin();
  double max_extent = *std::max_element(extents.begin(), extents.end());

  std::cout << bbox << std::endl;
  std::cout << bbox_center << std::endl;
  std::cout << max_extent << std::endl;
  
  const int NB_RAYS = 1000;
  std::vector<Point> v1, v2;
  v1.reserve(NB_RAYS); v2.reserve(NB_RAYS);

  const double r = max_extent / 2;
  // Generate NB_RAYS*2 points that lie on a sphere of radius r, centered around bbox_center
  CGAL::Random rand = CGAL::Random(23); // fix the seed to yield the same results each run
  CGAL::cpp11::copy_n(CGAL::Random_points_on_sphere_3<Point>(r, rand), NB_RAYS, std::back_inserter(v1));
  CGAL::cpp11::copy_n(CGAL::Random_points_on_sphere_3<Point>(r, rand), NB_RAYS, std::back_inserter(v2));

  for(std::vector<Point>::iterator it = v1.begin(); it != v1.end(); ++it) {
    *it = *it + bbox_center;
  }

  for(std::vector<Point>::iterator it = v2.begin(); it != v2.end(); ++it) {
    *it = *it + bbox_center;
  }
  
  // Generate NB_RAYS using v1 as source and v2 as target.
  std::vector<Ray> rays;
  rays.reserve(NB_RAYS);
  std::transform(v1.begin(), v1.end(), v2.begin(),
                 std::back_inserter(rays), boost::value_factory<Ray>());
  std::vector< boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type > > primitives1, primitives2;
  primitives1.reserve(NB_RAYS); primitives2.reserve(NB_RAYS);


  {
    for(std::vector<Ray>::iterator it = rays.begin(); it != rays.end(); ++it) {
      primitives1.push_back(min_intersection(tree, *it));
    }
  }

  for(std::vector<Ray>::iterator it = rays.begin(); it != rays.end(); ++it) {
    primitives2.push_back(tree.first_intersection(*it));
  }
  CGAL_assertion_msg(primitives1.size() == primitives2.size(), "Different amount of primitives intersected.");
  CGAL_assertion_msg(std::equal(primitives1.begin(), primitives1.end(), primitives2.begin()),
                     "Primitives mismatch.");
  std::size_t c = primitives1.size() - std::count(primitives1.begin(), primitives1.end(), boost::none);
  std::cout << "Intersected " << c << " primitives with " << NB_RAYS << " rays" << std::endl;
  std::cout << "Primitive method had to sort " << accum/NB_RAYS 
            << " intersections on average." << std::endl;
  t.stop();
  std::cout << t.time() << std::endl;
  return 0;
}
