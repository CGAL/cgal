#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/IO/Surface_mesh_parameterization/File_off.h>
#include <CGAL/parameterize.h>

#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>          Kernel;
typedef Kernel::Point_2                         Point_2;
typedef Kernel::Point_3                         Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>     SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor       face_descriptor;

int main(int argc, char * argv[])
{
  std::ifstream in((argc>1) ? argv[1] : "data/three_peaks.off");
  if(!in) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }

  SurfaceMesh sm;
  in >> sm;

  // An edge on the border
  halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  // The UV property map that holds the parameterized values
  typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2>  UV_pmap;
  UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("h:uv").first;

  CGAL::parameterize(sm, hd, uv_map);

  std::ofstream out("result.off");
  CGAL::Parameterization::output_uvmap_to_off(sm, hd, uv_map, out);

  return 0;
}
