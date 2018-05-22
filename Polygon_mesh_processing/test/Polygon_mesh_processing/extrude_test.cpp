#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SMesh;



int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/quad.off";
  SMesh in, out; 
  std::ifstream input(filename);
  
  if (!input || !(input >> in))
  {
    std::cerr << "Error: cannot read Surface Mesh : " << filename << "\n";
    assert(!CGAL::is_empty(in));
    assert(false);
    return 1;
  }
  CGAL::Polygon_mesh_processing::extrude_mesh(in, out, Kernel::Vector_3(0.0, 0.0, -1.0), 1.0);
  std::ofstream extruded_off("extruded.off");
  extruded_off << out;
  extruded_off.close();
  std::cerr << "All done." << std::endl;
  
  return 0;
}
