#include <CGAL/Simple_cartesian.h>

#include <CGAL/config.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/property_map.h>
#include <boost/iterator/function_output_iterator.hpp>

#include <cassert>
#include <fstream>
#include <string>
#include <vector>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef std::pair<Point_3, Vector_3> PointVectorPair;
typedef CGAL::First_of_pair_property_map<PointVectorPair> Point_map;
typedef CGAL::Second_of_pair_property_map<PointVectorPair> Normal_map;

int main()
{
  std::vector<PointVectorPair> pv_pairs;
  const std::function<void(const PointVectorPair& p)> lambda =
    [&](const PointVectorPair& p) {
    FT len = p.second.squared_length();
    if (len > 0 || len != 1.0) {
      Vector_3 n = p.second * (1.0 / CGAL::sqrt(len));
      pv_pairs.push_back(std::make_pair(p.first, n));
    }
    else pv_pairs.push_back(p);
    };

  pv_pairs.clear();
  std::ifstream file("data/simple_ascii.ply");
  CGAL::IO::read_PLY_with_properties<PointVectorPair>(file, boost::function_output_iterator(lambda),
    CGAL::make_ply_point_reader(Point_map()),
    CGAL::make_ply_normal_reader(Normal_map()));

  assert(pv_pairs[0].first == Point_3(1, 1, 1));
  assert(pv_pairs[1].first == Point_3(3, 3, 3));
  assert(pv_pairs[2].first == Point_3(5, 5, 5));

  for (std::size_t i = 0; i < pv_pairs.size(); i++) {
    FT dev = CGAL::abs(1.0 - pv_pairs[i].second.squared_length());
    assert(dev < 0.01);
  }

  return 0;
}
