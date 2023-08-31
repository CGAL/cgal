#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/acvd/acvd.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> Surface_Mesh;
typedef boost::graph_traits<Surface_Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::property_map<Surface_Mesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type VertexColorMap;


int main(int argc, char* argv[])
{
  Surface_Mesh smesh;
  const std::string filename = (argc > 1) ?
    CGAL::data_file_path(argv[1]) :
    CGAL::data_file_path("meshes/better_dragon.obj");

  const int nb_clusters = (argc > 2) ? atoi(argv[2]) : 98;

  if (!CGAL::IO::read_polygon_mesh(filename, smesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  PMP::acvd_isotropic_simplification(smesh, nb_clusters);
  // std::cout << "kak3" << std::endl;
  return 0;

  // Output the simplified mesh, use write_OFF()
  //CGAL::IO::write_OFF("sphere966_clustered_0.off", smesh, CGAL::parameters::stream_precision(17).vertex_color_map(vcm));

}

