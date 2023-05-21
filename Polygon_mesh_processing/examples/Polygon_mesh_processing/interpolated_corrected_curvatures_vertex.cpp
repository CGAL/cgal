#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> Surface_Mesh;
typedef boost::graph_traits<Surface_Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  // instantiating and reading mesh
  Surface_Mesh smesh;
  const std::string filename = (argc > 1) ?
    argv[1] :
    CGAL::data_file_path("meshes/sphere.off");

  if (!CGAL::IO::read_polygon_mesh(filename, smesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // loop over vertices and use vertex_descriptor to compute a curvature on one vertex
  for (vertex_descriptor v : vertices(smesh))
  {
    double h = PMP::interpolated_corrected_mean_curvature_one_vertex(smesh, v);
    double g = PMP::interpolated_corrected_Gaussian_curvature_one_vertex(smesh, v);
    PMP::Principal_curvatures_and_directions<Epic_kernel> p =
      PMP::interpolated_corrected_principal_curvatures_and_directions_one_vertex(smesh, v);

    // we can also specify a ball radius for expansion and a user defined vertex normals map using
    // named parameters. Refer to interpolated_corrected_curvatures_SM.cpp to see example usage.

    // Can also use interpolated_corrected_curvatures_one_vertex() to compute multiple curvatures
    // on the vertex at the same time. This is more efficient than computing each one separately.
    // The following commented lines show this (all mentioned named parameters work on it as well)
    // we specify which curvatures we want to compute by passing pointers as named parameters
    // as shown. These pointers are used for storing the result as well. in this example we
    // selected mean and Gaussian curvatures
    // PMP::interpolated_corrected_curvatures_one_vertex(
    //   smesh,
    //   v,
    //   CGAL::parameters::vertex_mean_curvature(&h)
    //   .vertex_Gaussian_curvature(&g)
    // );

    std::cout << v.idx() << ": HC = " << h
      << ", GC = " << g << "\n"
      << ", PC = [ " << p.min_curvature << " , " << p.max_curvature << " ]\n";
  }
}
