#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Unique_hash_map.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                           Kernel;
typedef Kernel::Point_2                                          Point_2;
typedef Kernel::Point_3                                          Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor       halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor         vertex_descriptor;

typedef CGAL::Unique_hash_map<vertex_descriptor, Point_2>        UV_uhm;
typedef boost::associative_property_map<UV_uhm>                  UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;

int main(int argc, char** argv)
{
  std::ifstream in((argc>1) ? argv[1] : "data/airplane_0.off");
  if(!in){
    std::cerr << "Error: problem loading the input data" << std::endl;
    return 1;
  }

  Surface_mesh sm;
  in >> sm;

  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  // The 2D points of the uv parametrisation will be written into this map
  UV_uhm uv_uhm;
  UV_pmap uv_map(uv_uhm);

  typedef SMP::Square_border_arc_length_parameterizer_3<Surface_mesh> Border_parameterizer;
  typedef SMP::Discrete_authalic_parameterizer_3<Surface_mesh, Border_parameterizer> Parameterizer;

  Border_parameterizer border_param; // the border parameterizer will compute the corner vertices

  SMP::Error_code err = SMP::parameterize(sm, Parameterizer(border_param), bhd, uv_map);

  if(err != SMP::OK) {
    std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
    return 1;
  }

  std::ofstream out("discrete_result.off");
  SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);

  return EXIT_SUCCESS;
}
