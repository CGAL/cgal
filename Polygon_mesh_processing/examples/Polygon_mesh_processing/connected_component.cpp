#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/Connected_components.h>
#include <iostream>
#include <fstream>


namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                           Mesh;


typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;


template <typename G>
struct Constraint {
  Constraint() { }

  Constraint(G & g) : g(&g) { }

  bool operator[](const edge_descriptor&) const {
    return false; // no constraint
  }

  G* g;
};

 
int main(int, char* argv[]) 
{
  Mesh sm;
  std::ifstream in(argv[1]);
  in >> sm;

  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(sm).first;
  CGAL::Polygon_mesh_processing::connected_component(fd,
                                                     sm,
                                                     std::back_inserter(cc));


  std::cerr << cc.size() << " faces in the CC of " << fd << std::endl;

  Mesh::Property_map<face_descriptor,std::size_t> fccmap;
  fccmap = sm.add_property_map<face_descriptor,std::size_t>("f:CC").first; 
  std::size_t num = PMP::connected_components(sm,
                                              fccmap,
                                              CGAL::parameters::edge_is_constrained_map(Constraint<Mesh>())
                                              );
  
 std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
 BOOST_FOREACH(face_descriptor f , faces(sm)){
   std::cout  << f << " in connected component " << fccmap[f] << std::endl;
  }
 
 CGAL::Polygon_mesh_processing::keep_largest_connected_components(sm,2);

 std::cout << "mesh:\n" << sm << std::endl;
  return 0;
}
