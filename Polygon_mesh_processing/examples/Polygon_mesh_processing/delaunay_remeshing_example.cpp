#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/surface_Delaunay_remeshing.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Mesh_constant_domain_field_3.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

typedef CGAL::Mesh_constant_domain_field_3<K, int> Sizing_field;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  std::string filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/anchor_dense.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  double target_edge_length = (argc > 2) ? std::stod(std::string(argv[2])) : 0.02;
  Sizing_field size(target_edge_length);

  double fdist = (argc > 3) ? std::stod(std::string(argv[3])) : 0.01;

  std::cout << "Detect features..." << std::endl;

  using EIFMap = boost::property_map<Mesh, CGAL::edge_is_feature_t>::type;
  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  PMP::detect_sharp_edges(mesh, 45, eif);

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  Mesh outmesh = PMP::surface_Delaunay_remeshing(mesh,
                                  CGAL::parameters::protect_constraints(true)
                                 .mesh_edge_size(size)
                                 .mesh_facet_distance(fdist)
                                 .edge_is_constrained_map(eif));

  std::cout << "Remeshing done." << std::endl;

  std::ofstream ofs("anchor_remeshed.off");
  ofs << outmesh;
  ofs.close();

  return 0;
}
