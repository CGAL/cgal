#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>

#include <CGAL/Polygon_mesh_processing/measure.h>


#include <cstdlib>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>       Kernel;
typedef Kernel::Point_2                      Point_2;
typedef Kernel::Point_3                      Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>  SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor      face_descriptor;

namespace SMP = CGAL::Surface_mesh_parameterization;

int main(int argc, char** argv)
{
  std::ifstream in((argc>1) ? argv[1] : "data/three_peaks.off");
  if(!in) {
    std::cerr << "Problem loading the input data" << std::endl;
    return EXIT_FAILURE;
  }

  SurfaceMesh sm;
  in >> sm;

  // A halfedge on the border
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  // The 2D points of the uv parametrisation will be written into this map
  typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2>  UV_pmap;
  UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("v:uv").first;

  typedef SMP::Circular_border_arc_length_parameterizer_3<SurfaceMesh>  Border_parameterizer;
  typedef SMP::Discrete_authalic_parameterizer_3<SurfaceMesh, Border_parameterizer> Parameterizer;

  SMP::Error_code err = SMP::parameterize(sm, Parameterizer(), bhd, uv_map);

  if(err != SMP::OK) {
    std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
  }

  std::ofstream out("result.off");
  SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);

  return EXIT_SUCCESS;
}
