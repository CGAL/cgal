#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>
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
#ifdef CGAL_SMP_USE_SUITESPARSE_SOLVERS
  std::cout << "Using SuiteSparse" << std::endl;
#else
  std::cout << "Using Eigen3" << std::endl;
#endif
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/head.off");

  SurfaceMesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // A halfedge on the border
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  // The 2D points of the uv parametrisation will be written into this map
  typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2>  UV_pmap;
  UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("v:uv").first;

  SMP::ARAP_parameterizer_3<SurfaceMesh> parameterizer;
  SMP::Error_code err = SMP::parameterize(sm, parameterizer, bhd, uv_map);

  if(err != SMP::OK) {
    std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
  }

  std::ofstream out("result.off");
  SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);

  return EXIT_SUCCESS;
}
