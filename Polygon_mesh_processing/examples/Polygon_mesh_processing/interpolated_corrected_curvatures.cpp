#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_corrected_curvature_measures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <unordered_map>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel EpicKernel;
typedef CGAL::Surface_mesh<EpicKernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  Mesh g1;
  const std::string filename = (argc>1) ?
      argv[1] :
      CGAL::data_file_path("meshes/small_bunny.obj");

  if(!CGAL::IO::read_polygon_mesh(filename, g1))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  bool created = false;

  Mesh::Property_map<vertex_descriptor, EpicKernel::FT> mean_curvature_map, gaussian_curvature_map;
  boost::tie(mean_curvature_map, created) = g1.add_property_map<vertex_descriptor, EpicKernel::FT>("v:mean_curvature_map", 0);
  assert(created);

  boost::tie(gaussian_curvature_map, created) = g1.add_property_map<vertex_descriptor, EpicKernel::FT>("v:gaussian_curvature_map", 0);
  assert(created);

  Mesh::Property_map<vertex_descriptor, std::tuple<
      EpicKernel::FT,
      EpicKernel::FT,
      EpicKernel::Vector_3,
      EpicKernel::Vector_3
      >> principle_curvature_map;

  boost::tie(principle_curvature_map, created) = g1.add_property_map<vertex_descriptor, std::tuple<
      EpicKernel::FT,
      EpicKernel::FT,
      EpicKernel::Vector_3,
      EpicKernel::Vector_3
      >>("v:principle_curvature_map", { 0, 0,
          EpicKernel::Vector_3 (0,0,0),
          EpicKernel::Vector_3 (0,0,0)});
  assert(created);

  PMP::interpolated_corrected_mean_curvature(
      g1,
      mean_curvature_map
  );
  PMP::interpolated_corrected_gaussian_curvature(
      g1,
      gaussian_curvature_map
  );
  PMP::interpolated_corrected_principal_curvatures(
      g1,
      principle_curvature_map
  );

  for (vertex_descriptor v : vertices(g1))
  {
    auto PC = principle_curvature_map[v];
      std::cout << v.idx() << ": HC = " << mean_curvature_map[v]
                           << ", GC = " << gaussian_curvature_map[v] << "\n"
                           << ", PC = [ " << std::get<0>(PC) << " , " << std::get<1>(PC) << " ]\n";
  }
}
