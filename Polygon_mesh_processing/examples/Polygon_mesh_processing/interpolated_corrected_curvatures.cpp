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
typedef std::unordered_map<face_descriptor, EpicKernel::FT> FaceMeasureMap_tag;
typedef std::unordered_map<vertex_descriptor, EpicKernel::FT> vertexMeasureMap_tag;


int main(int argc, char* argv[])
{
  Mesh g1;
  const std::string filename = (argc>1) ?
      argv[1] :
      CGAL::data_file_path("meshes/cylinder.off");

  if(!CGAL::IO::read_polygon_mesh(filename, g1))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }


  vertexMeasureMap_tag mean_curvature_init, gaussian_curvature_init;
  boost::associative_property_map<vertexMeasureMap_tag>
      mean_curvature_map(mean_curvature_init), gaussian_curvature_map(gaussian_curvature_init);

  PMP::interpolated_corrected_mean_curvature(
      g1,
      mean_curvature_map
      );
  PMP::interpolated_corrected_gaussian_curvature(
      g1,
      gaussian_curvature_map
  );

  for (vertex_descriptor v : vertices(g1))
      std::cout << v.idx() << ": HC = " << get(mean_curvature_map, v)
                           << ", GC = " << get(gaussian_curvature_map, v) << "\n";
}
