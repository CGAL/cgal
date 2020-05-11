#include <vector>
#include <cassert>
#include <string>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/Point_set_3/point_set_io.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef std::pair<Point_3, Vector_3> PointVectorPair;

template <typename PointRange>
void fill_points(PointRange& points)
{
  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(1,0,0));
  points.push_back(Point_3(0,1,0));
  points.push_back(Point_3(0,0,1));
  points.push_back(Point_3(1,1,1));
}

bool write_ps(std::string s)
{
  std::vector<Point_3> points;
  fill_points(points);

  return CGAL::write_points(s, points);
}

//todo
bool write_ps_with_np(std::string s)
{
  return false;
}
bool read_point_set(std::string s)
{
  CGAL::Point_set_3<Point_3, Vector_3> ps;
  return CGAL::read_point_set(s, ps);
}

bool read_ps(std::string s)
{
  std::vector<Point_3> points;
  return CGAL::read_points(s, std::back_inserter(points));
}

bool read_ps_with_np(std::string s)
{
  std::vector<PointVectorPair> pv_pairs;
  return CGAL::read_points(s,
                           std::back_inserter(pv_pairs),
                           CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                           normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
}


int main()
{
  assert(write_ps("test.xyz"));
  assert(write_ps("test.off"));
  assert(write_ps("test.ply"));
#ifdef CGAL_LINKED_WITH_LASLIB
  assert(write_ps("test.las"));
#endif

  assert(read_point_set("data/read_test/ok_1.xyz"));
  assert(read_point_set("data/read_test/ok_1.off"));
  assert(read_point_set("data/read_test/simple.ply"));
#ifdef CGAL_LINKED_WITH_LASLIB
  assert(read_point_set("data/read_test/pig_points.las"));
#endif

  assert(read_ps("data/read_test/ok_1.xyz"));
  assert(read_ps("data/read_test/ok_1.off"));
  assert(read_ps("data/read_test/simple.ply"));
#ifdef CGAL_LINKED_WITH_LASLIB
  assert(read_ps("data/read_test/pig_points.las"));
#endif

  assert(read_ps_with_np("data/read_test/ok_1.xyz"));
  assert(read_ps_with_np("data/read_test/ok_1.off"));
  assert(read_ps_with_np("data/read_test/simple.ply"));
#ifdef CGAL_LINKED_WITH_LASLIB
  assert(read_ps_with_np("data/read_test/pig_points.las"));
#endif



  return 0;
}
