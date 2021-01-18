//#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

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


int main()
{

std::vector<PointWithColor> points(3);
points[0] = std::make_pair(Point_3(1,0,0), Color(255,0,0,255));
points[1] = std::make_pair(Point_3(0,1,0), Color(0,255,0,255));
points[2] = std::make_pair(Point_3(0,0,1), Color(0,0,255,255));


std::ofstream os;
std::ifstream is;
#ifdef CGAL_LINKED_WITH_LASLIB
os.open("tmp.las", std::ios::binary);
bool ok = CGAL::write_las_points_with_properties(os, points, std::make_tuple(
                                                   ),
                                                 );
is.open("data/read_test/pig_points.las", std::ios::binary);


assert(ok);
#endif

}
