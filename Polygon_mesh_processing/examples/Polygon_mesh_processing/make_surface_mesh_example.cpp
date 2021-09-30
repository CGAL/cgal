#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/make_surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/anchor_dense.off";

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  double target_edge_length = (argc > 2) ? std::stod(std::string(argv[2])) : 0.02;

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  Mesh outmesh;
  PMP::make_surface_mesh(mesh, outmesh,
    PMP::parameters::protect_constraints(true)
    .mesh_edge_size(0.025));

  std::cout << "Remeshing done." << std::endl;

  std::ofstream ofs("anchor_remeshed.off");
  ofs << outmesh;
  ofs.close();

  return 0;
}
