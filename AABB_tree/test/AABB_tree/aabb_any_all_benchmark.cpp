#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/function_objects.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polyhedron_3.h>


#include <CGAL/Timer.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>


const std::size_t elements = 50000;
const int runs = 10;

template<typename Tree>
struct FilterP {
  const Tree* t;

  template<typename T>
  bool operator()(const T& tt) { return !t->do_intersect(tt); }
};


template<typename ForwardIterator, typename Tree>
std::size_t intersect(ForwardIterator b, ForwardIterator e, const Tree& tree, long& counter) {
      typedef
        typename Tree::AABB_traits::template Intersection_and_primitive_id<typename ForwardIterator::value_type>::Type
        Obj_type;

  std::vector<Obj_type> v;
  // bad educated guess
  v.reserve(elements);
  for(; b != e; ++b) {
    tree.all_intersections(*b, std::back_inserter(v));
    boost::optional<Obj_type> o = tree.any_intersection(*b);
    if(o)
      ++counter;
  }

  return v.size();
}

template<typename K>
boost::tuple<std::size_t, std::size_t, std::size_t, long> test(const char* name) {
  typedef typename K::FT FT;
  typedef typename K::Ray_3 Ray;
  typedef typename K::Line_3 Line;
  typedef typename K::Point_3 Point;
  typedef typename K::Segment_3 Segment;
  typedef CGAL::Polyhedron_3<K> Polyhedron;

  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  std::ifstream ifs(name);

  Polyhedron polyhedron;
  ifs >> polyhedron;

  // Random seeded to 23, cube size equal to the magic number 2
  CGAL::Random r(23);
  CGAL::Random_points_in_cube_3<Point, CGAL::Creator_uniform_3<FT, Point> > g( 2., r);

  std::vector<Point> points;
  points.reserve(elements * 2);
  std::copy_n(g, elements * 2, std::back_inserter(points));

  // generate a bunch of happy random primitives
  std::vector<Line> lines;
  lines.reserve(elements);

  // forward
  for(std::size_t i = 0; i < points.size(); i += 2)
  {
    lines.push_back(Line(points[i], points[i + 1]));
  }

  std::vector<Ray> rays;
  rays.reserve(elements);

  // backwards
  for(std::size_t i = points.size(); i != 0; i -= 2)
  {
    rays.push_back(Ray(points[i - 1], points[i - 2]));
  }

  std::vector<Segment> segments;
  segments.reserve(elements);
  // from both sides
  for(std::size_t i = 0, j = points.size() - 1; i < j; ++i, --j)
  {
    segments.push_back(Segment(points[i], points[j]));
  }

  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);

  // filter all primitives that do not intersect

  FilterP<Tree> p = { &tree };

  lines.erase(std::remove_if(lines.begin(), lines.end(), p), lines.end());

  rays.erase(std::remove_if(rays.begin(), rays.end(), p), rays.end());

  segments.erase(std::remove_if(segments.begin(), segments.end(), p), segments.end());

  boost::tuple<std::size_t, std::size_t, std::size_t, long> tu;

    {
      CGAL::Timer t;
      t.start();

      for(int i = 0; i < runs; ++i)
      {
        long counter = 0L;
        tu = boost::make_tuple(intersect(lines.begin(), lines.end(), tree, counter),
                               intersect(rays.begin(), rays.end(), tree, counter),
                               intersect(segments.begin(), segments.end(), tree, counter),
                               // cant use counter here
                               0);
        boost::get<3>(tu) = counter;
      }
      std::cout << t.time();
    }

  return tu;
}

int main()
{
  const char* filename = "./data/finger.off";

  std::cout << "| Simple cartesian float kernel | ";
  boost::tuple<std::size_t, std::size_t, std::size_t, long> t1 = test<CGAL::Simple_cartesian<float> >(filename);
  std::cout << " | " << std::endl;

  std::cout << "| Cartesian float kernel | ";
  boost::tuple<std::size_t, std::size_t, std::size_t, long> t2 = test<CGAL::Cartesian<float> >(filename);
  std::cout << " | " << std::endl;

  std::cout << "| Simple cartesian double kernel |";
  boost::tuple<std::size_t, std::size_t, std::size_t, long> t3 = test<CGAL::Simple_cartesian<double> >(filename);
  std::cout << " | " << std::endl;

  std::cout << "| Cartesian double kernel |";
  boost::tuple<std::size_t, std::size_t, std::size_t, long> t4 = test<CGAL::Cartesian<double> >(filename);
  std::cout << " | " << std::endl;

  std::cout << "| Epic kernel |";
  boost::tuple<std::size_t, std::size_t, std::size_t, long> t5 = test<CGAL::Exact_predicates_inexact_constructions_kernel>(filename);
  std::cout << " | " << std::endl;

  std::size_t a, b, c;
  long d;

  boost::tie(a, b, c, d) = t5;
  std::cout << a << " " << b << " " << c << " " << d << std::endl;

  boost::tie(a, b, c, d) = t4;
  std::cout << a << " " << b << " " << c << " " << d << std::endl;

  boost::tie(a, b, c, d) = t3;
  std::cout << a << " " << b << " " << c << " " << d << std::endl;

  boost::tie(a, b, c, d) = t2;
  std::cout << a << " " << b << " " << c << " " << d << std::endl;

  boost::tie(a, b, c, d) = t1;
  std::cout << a << " " << b << " " << c << " " << d << std::endl;
  return 0;
}
