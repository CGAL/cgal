#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Heat_method_3/Intrinsic_Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Point_2                                      Point_2;
typedef CGAL::Surface_mesh<Point>                            Mesh;
//typedef CGAL::Polyhedron_3<Kernel> Mesh;

typedef CGAL::dynamic_vertex_property_t<double> Vertex_distance_tag;
typedef boost::property_map<Mesh, Vertex_distance_tag >::type Vertex_distance_map;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<Mesh> IDT;


void bglstyle(const IDT& idt)
{
  typedef boost::graph_traits<IDT>::vertex_descriptor IDT_vertex_descriptor;
  typedef boost::graph_traits<IDT>::halfedge_descriptor IDT_halfedge_descriptor;

  std::cout << num_vertices(idt) << std::endl;
  vertices(idt);
  IDT_halfedge_descriptor hd = *(halfedges(idt).first);
  IDT_vertex_descriptor vd = vertex(hd,idt);
  CGAL_USE(vd);

  //get(idt, CGAL::vertex_point);
}


int main()
{
  Mesh sm;

  std::ifstream in("data/brain100k.off");
  in >> sm;
  if(!in || num_vertices(sm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  num_faces(sm);
  IDT im(sm);

  std::cout<<"success \n";
  return 0;

  //tests to add: for all edges, angles opposite edge sum up to less than pi
}
