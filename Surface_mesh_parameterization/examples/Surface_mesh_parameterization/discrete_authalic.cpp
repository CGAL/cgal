#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/Circular_border_parameterizer_3.h>
#include <CGAL/Square_border_parameterizer_3.h>
#include <CGAL/IO/Surface_mesh_parameterization/File_off.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/parameterize.h>

#include <CGAL/Polygon_mesh_processing/measure.h>

#include <boost/foreach.hpp>

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

int main(int argc, char * argv[])
{
  SurfaceMesh sm;

  std::ifstream in((argc>1) ? argv[1] : "data/nefertiti.off");

  in >> sm;

  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  // The 2D points of the uv parametrisation will be written into this map
  SurfaceMesh::Property_map<vertex_descriptor, Point_2> uvpm;
  bool created;
  boost::tie(uvpm, created) = sm.add_property_map<vertex_descriptor, Point_2>("v:uv");

  // Four different border parameterizers (pick one)
  typedef CGAL::Circular_border_uniform_parameterizer_3<SurfaceMesh>     Border_parameterizer;
//  typedef CGAL::Circular_border_arc_length_parameterizer_3<SurfaceMesh>  Border_parameterizer;
//  typedef CGAL::Square_border_uniform_parameterizer_3<SurfaceMesh>       Border_parameterizer;
//  typedef CGAL::Square_border_arc_length_parameterizer_3<SurfaceMesh>    Border_parameterizer;

  typedef CGAL::Discrete_authalic_parameterizer_3<SurfaceMesh, Border_parameterizer> Parameterizer;
  Parameterizer::Error_code err = CGAL::parameterize(sm, Parameterizer(), bhd, uvpm);

  if(err != Parameterizer::OK) {
    std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
    return 1;
  }

  std::ofstream out("result.off");
  CGAL::Parameterization::output_uvmap_to_off(sm, bhd, uvpm, out);

  return 0;
}
