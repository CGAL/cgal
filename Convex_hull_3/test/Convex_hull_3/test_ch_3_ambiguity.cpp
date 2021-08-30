#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/convex_hull_3.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_3                                           Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;

int main(int argc, char* argv[])
{
  const char* filename = (argc>1)? argv[1] : "data/cross.off";

  Surface_mesh poly;
  if(!CGAL::IO::read_polygon_mesh(filename, poly))
  {
    std::cerr<<"Could not find a correct input file."<<std::endl;
    return 1;
  }

  Surface_mesh chull;
  // compute convex hull
  auto np = CGAL::parameters::all_default();
  CGAL::convex_hull_3(poly, chull, np);
  std::cout << "The convex hull contains " << chull.number_of_vertices() << " vertices" << std::endl;
  return 0;
}
