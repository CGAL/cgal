#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>
#include <CGAL/Heat_method_3/Intrinsic_Delaunay_triangulation_3.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Point_2                                      Point_2;
typedef CGAL::Surface_mesh<Point>                            Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<Surface_mesh,Kernel, Vertex_distance_map> Idt;

typedef CGAL::Heat_method_3::Heat_method_3<Idt,Kernel, Idt::Vertex_distance_map> Heat_method_idt;


int main(int argc, char* argv[])
{
  //read in mesh
  Surface_mesh sm;
  const char* filename = (argc > 1) ? argv[1] : "./data/bunny.off";
  std::ifstream in(filename);
  in >> sm;
  //the vertex distance map will hold the distance values from the source set
  Vertex_distance_map vdm_idt = sm.add_property_map<vertex_descriptor,double>("v:idt",0).first;
  //first call the idt remeshing algorithm
  Idt idt(sm, vdm_idt);

  //pass in the idt object and its vertex_distance_map
  Heat_method_idt hm_idt(idt, idt.vertex_distance_map());

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(sm).first);
  hm_idt.add_source(source);
  hm_idt.update();

  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    std::cout << vd << "  is at distance " << get(vdm_idt, vd) << " from " << source << std::endl;
  }

  return 0;
}
