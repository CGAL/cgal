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


bool write_ps(std::string s)
{
  std::vector<Point_3> points;
  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(1,0,0));
  points.push_back(Point_3(0,1,0));
  points.push_back(Point_3(0,0,1));
  points.push_back(Point_3(1,1,1));

  return CGAL::write_points(s, points);
}

bool write_ps_with_np(std::string s)
{
  std::vector<PointVectorPair> points;
  points.push_back(std::make_pair(Point_3(0,0,0), Vector_3(0,0,0)));
  points.push_back(std::make_pair(Point_3(1,0,0), Vector_3(1,0,0)));
  points.push_back(std::make_pair(Point_3(0,1,0), Vector_3(0,1,0)));
  points.push_back(std::make_pair(Point_3(0,0,1), Vector_3(0,0,1)));
  points.push_back(std::make_pair(Point_3(1,1,1), Vector_3(1,1,1)));

  return CGAL::write_points(s, points,
                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                            normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
}

bool write_point_set(std::string s)
{
  CGAL::Point_set_3<Point_3, Vector_3> point_set;
  point_set.insert (Point_3(0,0,0));
  point_set.insert (Point_3(1,0,0));
  point_set.insert (Point_3(0,1,0));
  point_set.insert (Point_3(0,0,1));
  point_set.insert (Point_3(1,1,1));

  point_set.add_normal_map();
  return CGAL::write_point_set(s, point_set);
}

bool read_point_set(std::string s)
{
  CGAL::Point_set_3<Point_3, Vector_3> ps;
  if(!CGAL::read_point_set(s, ps))
    return false;
  if(ps.size() != 5)
  {
    std::cerr <<"error: wrong number of points read."<<std::endl;
    return false;
  }
  return true;
}

bool read_ps(std::string s)
{
  std::vector<Point_3> points;
  if(!CGAL::read_points(s, std::back_inserter(points)))
    return false;
  if(points.size() != 5)
  {
    std::cerr <<"error: wrong number of points read."<<std::endl;
    return false;
  }
  return true;
}

bool read_ps_with_np(std::string s)
{
  std::vector<PointVectorPair> pv_pairs;
  if(!CGAL::read_points(s,
                        std::back_inserter(pv_pairs),
                        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                        normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())))
    return false;
  if(pv_pairs.size() != 5)
  {
    std::cerr <<"error: wrong number of points read."<<std::endl;
    return false;
  }
  return true;
}


int main()
{
  std::vector<std::string> filenames = {"test.xyz" ,
                                        "test.off" ,
                                        "test.ply"};
#ifdef CGAL_LINKED_WITH_LASLIB
  filenames.push_back("test.las");
#endif
  for(const std::string filename : filenames)
  {
    assert(write_ps(filename));
    assert(read_ps(filename));
    assert(write_ps_with_np(filename));
    assert(read_ps_with_np(filename));
    assert(write_point_set(filename));
    assert(read_point_set(filename));
  }
  return 0;
}
