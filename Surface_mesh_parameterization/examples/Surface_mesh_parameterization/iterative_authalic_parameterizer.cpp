#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Unique_hash_map.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                           Kernel;
typedef Kernel::Point_2                                          Point_2;
typedef Kernel::Point_3                                          Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor   halfedge_descriptor;

typedef CGAL::Unique_hash_map<vertex_descriptor, Point_2>        UV_uhm;
typedef boost::associative_property_map<UV_uhm>                  UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;

int main(int argc, char** argv)
{
  std::ifstream in((argc>1) ? argv[1] : "data/nefertiti.off");
  if(!in)
  {
    std::cerr << "Error: problem loading the input data" << std::endl;
    return EXIT_FAILURE;
  }

  Surface_mesh sm;
  in >> sm;

  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  // The 2D points of the uv parametrisation will be written into this map
  UV_uhm uv_uhm;
  UV_pmap uv_map(uv_uhm);

  typedef SMP::Circular_border_arc_length_parameterizer_3<Surface_mesh> Border_parameterizer;
  Border_parameterizer border_parameterizer; // the border parameterizer will automatically compute the corner vertices

  typedef SMP::Iterative_authalic_parameterizer_3<Surface_mesh, Border_parameterizer> Parameterizer;
  Parameterizer parameterizer(border_parameterizer);

  const unsigned int iterations = (argc > 2) ? std::atoi(argv[2]) : 15;
  SMP::Error_code err = parameterizer.parameterize(sm, bhd, uv_map, iterations);

  if(err != SMP::OK)
  {
    std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
  }

  std::ofstream out("iterative_result.off");
  SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);

  return EXIT_SUCCESS;
}
