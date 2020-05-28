#include <iostream>
#include <fstream>

#include <CGAL/IO/OFF.h>

#include <CGAL/Simple_cartesian.h>
#include <vector>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef std::vector<std::size_t> Face;
int main()
{
  std::ifstream in("data/cube.off");
  std::vector<Point> points;
  std::vector<Face> faces;
  CGAL::read_OFF(in, points, faces);
  in.close();
  assert(points.size() == 8);
  assert(faces.size() == 12);

  assert(CGAL::write_OFF(std::cout, points, faces));

  return 0;
}
