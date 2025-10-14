#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/IO/VTK.h>

#include <vector>


template <typename Kernel>
void test_VTK()
{
  typedef Kernel::Point_3                               Point;
  typedef std::vector<std::size_t>                      Face;

  const std::vector<Point> points = { Point(0,0,0), Point(1,0,0), Point(0,1,0), Point(0,0,1) };
  const std::vector<Face> polygons = { Face{0,1,2}, Face{0,1,3}, Face{0,2,3}, Face{1,2,3} };

  bool ok = CGAL::IO::write_VTP("tmp.vtp", points, polygons, CGAL::parameters::use_binary_mode(true));
  assert(ok);

  std::vector<Point> rpoints;
  std::vector<Face> rpolygons;
  ok = CGAL::IO::read_VTP("tmp.vtp", rpoints, rpolygons, CGAL::parameters::use_binary_mode(true));
  assert(points == rpoints);
  assert(polygons == rpolygons);

  ok = CGAL::IO::write_VTP("tmp2.vtp", points, polygons, CGAL::parameters::use_binary_mode(false));
  assert(ok);

  rpoints.clear();
  rpolygons.clear();
  ok = CGAL::IO::read_VTP("tmp2.vtp", rpoints, rpolygons, CGAL::parameters::use_binary_mode(false));
  assert(points == rpoints);
  assert(polygons == rpolygons);

}

int main(int argc, char** argv)
{
  test_VTK<CGAL::Simple_cartesian<double>>();
  test_VTK<CGAL::Simple_cartesian<CGAL::Exact_rational>>();
  test_VTK<CGAL::Exact_predicates_exact_constructions_kernel>();
  return 0;
}