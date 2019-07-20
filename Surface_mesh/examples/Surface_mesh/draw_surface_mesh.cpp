#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

int main(int argc, char* argv[])
{
  Mesh sm1;
  std::ifstream in1((argc>1)?argv[1]:"data/elephant.off");
  in1 >> sm1;

  CGAL::draw(sm1);

  return EXIT_SUCCESS;
}
