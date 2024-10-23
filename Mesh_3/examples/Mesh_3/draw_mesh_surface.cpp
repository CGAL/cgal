#include "draw_c3t3_surface.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>

using Kernel= CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;

using Polyhedron = Mesh;

#include "ui_draw_c3t3_surface.h"

int main(int argc, char* argv[])
{
  Mesh mesh;
  std::ifstream in1((argc>1)?argv[1]:"data/elephant.off");
  std::string comments;
  CGAL::IO::read_PLY(in1, mesh, comments);
  std::cerr << "comments: " << comments << '\n';

  return draw_c3t3_surface_from_surface_mesh(mesh);
}
