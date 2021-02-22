#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>

// Just to try and create ambiguities
#include <CGAL/boost/graph/io.h>
#include <CGAL/IO/io.h>

#include <CGAL/property_map.h>

#include <vector>
#include <cassert>
#include <string>
#include <fstream>

typedef CGAL::Simple_cartesian<double>         Kernel;
typedef Kernel::Point_3                        Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef std::array<unsigned short, 4> Color;
typedef std::pair<Point_3, Color> PointWithColor;
typedef CGAL::Nth_of_tuple_property_map<1, PointWithColor> Color_map;

struct GetRedMap{
  typedef PointWithColor key_type;
  typedef unsigned short value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;
};
unsigned short get(const GetRedMap&, const PointWithColor& p)
{
  return p.second[0];
}

struct GetGreenMap{
  typedef PointWithColor key_type;
  typedef unsigned short value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;
};
unsigned short get(const GetGreenMap&, const PointWithColor& p)
{
  return p.second[1];
}

struct GetBlueMap{
  typedef PointWithColor key_type;
  typedef unsigned short value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;
};
unsigned short get(const GetBlueMap&, const PointWithColor& p)
{
  return p.second[2];
}

struct GetAlphaMap{
  typedef PointWithColor key_type;
  typedef unsigned short value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;
};
unsigned short get(const GetAlphaMap&, const PointWithColor& p)
{
  return p.second[3];
}

int main()
{

std::vector<PointWithColor> points(3);
points[0] = std::make_pair(Point_3(1,0,0), Color{255,0,0,255});
points[1] = std::make_pair(Point_3(0,1,0), Color{0,255,0,255});
points[2] = std::make_pair(Point_3(0,0,1), Color{0,0,255,255});


std::ofstream os;
std::ifstream is;
bool ok;
std::vector<Point_3> ps;
ps.push_back(Point_3(1,0,0));
ps.push_back(Point_3(0,1,0));
ps.push_back(Point_3(0,0,1));
//LAS
#ifdef CGAL_LINKED_WITH_LASLIB
os.open("tmp.las", std::ios::binary);
ok = CGAL::write_las_points_with_properties(os, points,
                                                 CGAL::make_las_point_writer(CGAL::First_of_pair_property_map<PointWithColor>()),
                                                 std::make_pair(GetRedMap(),CGAL::LAS_property::R()),
                                                 std::make_pair(GetGreenMap(), CGAL::LAS_property::G()),
                                                 std::make_pair(GetBlueMap(), CGAL::LAS_property::B()),
                                                 std::make_pair(GetAlphaMap(), CGAL::LAS_property::I())
                                                 );
os.close();
assert(ok);
points.clear();
is.open("tmp.las", std::ios::binary);
ok = CGAL::read_las_points_with_properties(is, std::back_inserter (points),
                                    CGAL::make_las_point_reader(CGAL::First_of_pair_property_map<PointWithColor>()),
                                    std::make_tuple(CGAL::Second_of_pair_property_map<PointWithColor>(),
                                                    CGAL::Construct_array(),
                                                    CGAL::LAS_property::R(),
                                                    CGAL::LAS_property::G(),
                                                    CGAL::LAS_property::B(),
                                                    CGAL::LAS_property::I()));
is.close();
assert(ok);
assert(points.size() == 3);
assert(points[1].second[1] == 255);

os.open("tmp.las", std::ios_base::binary);
CGAL::write_las_points(os, ps, CGAL::parameters::all_default());
os.close();
assert(ok);
ps.clear();
is.open("tmp.las", std::ios::binary);
ok = CGAL::read_las_points(is, std::back_inserter (ps),CGAL::parameters::all_default());
is.close();
assert(ok);
assert(ps.size() == 3);
#endif
//PLY
os.open("tmp.ply");
ok = CGAL::write_ply_points_with_properties(os, points,
                                            CGAL::make_ply_point_writer (CGAL::First_of_pair_property_map<PointWithColor>()),
                                            std::make_pair(GetRedMap(),CGAL::PLY_property<unsigned short>("red")),
                                            std::make_pair(GetGreenMap(), CGAL::PLY_property<unsigned short>("green")),
                                            std::make_pair(GetBlueMap(), CGAL::PLY_property<unsigned short>("blue")),
                                            std::make_pair(GetAlphaMap(), CGAL::PLY_property<unsigned short>("alpha"))
                                            );
os.close();
assert(ok);

is.open("tmp.ply");
points.clear();
ok = CGAL::read_ply_points_with_properties(is, std::back_inserter (points),
                                    CGAL::make_ply_point_reader(CGAL::First_of_pair_property_map<PointWithColor>()),
                                    std::make_tuple(CGAL::Second_of_pair_property_map<PointWithColor>(),
                                                    CGAL::Construct_array(),
                                                    CGAL::PLY_property<unsigned short>("red"),
                                                    CGAL::PLY_property<unsigned short>("green"),
                                                    CGAL::PLY_property<unsigned short>("blue"),
                                                    CGAL::PLY_property<unsigned short>("alpha")));
is.close();
assert(ok);
assert(points.size() == 3);
assert(points[1].second[1] == 255);

os.open("tmp.ply");
ok = CGAL::write_ply_points(os, ps, CGAL::parameters::all_default());
os.close();
assert(ok);

is.open("tmp.ply");
ps.clear();
ok = CGAL::read_ply_points(is, std::back_inserter (ps),
                           CGAL::parameters::all_default());
is.close();
assert(ok);
assert(ps.size() == 3);

//OFF
os.open("tmp.off");
ok = CGAL::write_off_points(os, ps, CGAL::parameters::all_default());
os.close();
assert(ok);

is.open("tmp.off");
ps.clear();
ok = CGAL::read_off_points(is, std::back_inserter (ps),
                           CGAL::parameters::all_default());
is.close();
assert(ok);
assert(ps.size() == 3);

//XYZ
os.open("tmp.xyz");
ok = CGAL::write_xyz_points(os, ps, CGAL::parameters::all_default());
os.close();
assert(ok);

is.open("tmp.xyz");
ps.clear();
ok = CGAL::read_xyz_points(is, std::back_inserter (ps),
                           CGAL::parameters::all_default());
is.close();
assert(ok);
assert(ps.size() == 3);
}
