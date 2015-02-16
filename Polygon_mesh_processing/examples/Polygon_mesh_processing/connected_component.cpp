#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/Connected_components.h>
#include <iostream>
#include <fstream>


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

  bool operator[](const edge_descriptor& e) const {
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
                                                     Constraint<Mesh>(sm),
                                                     std::back_inserter(cc));


  std::cerr << cc.size() << " faces in the CC of " << fd << std::endl;

  Mesh::Property_map<face_descriptor,int> fccmap;
  fccmap = sm.add_property_map<face_descriptor,int>("f:CC").first; 
  int num = CGAL::Polygon_mesh_processing::connected_components(sm,
                                                                //Constraint<Mesh>(sm),
                                                                fccmap);
  
 std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
 BOOST_FOREACH(face_descriptor f , faces(sm)){
   std::cout  << f << " in connected component " << fccmap[f] << std::endl;
  }
 
 CGAL::Polygon_mesh_processing::keep_largest_connected_components(sm,2);

 std::cout << "mesh:\n" << sm << std::endl;
  return 0;
}
