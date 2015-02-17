#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/Connected_components.h>
#include <iostream>
#include <fstream>


namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>                Mesh;


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

  int i=0;
  BOOST_FOREACH(face_descriptor f, faces(sm)){
    f->id() = i++;
  } 
  i=0;
  BOOST_FOREACH(vertex_descriptor v, vertices(sm)){
    v->id() = i++;
  }

  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(sm).first;
  CGAL::Polygon_mesh_processing::connected_component(fd,
                                                     sm,
                                                     std::back_inserter(cc));


  std::cerr << cc.size() << " faces in the CC of " << &*fd << std::endl;

  boost::vector_property_map<int, typename boost::property_map<Mesh, boost::face_index_t>::type> fccmap(get(boost::face_index,sm));

  int num = PMP::connected_components(sm,
                                      fccmap// ,
                                      //      CGAL::parameters::edge_is_constrained_map(Constraint<Mesh>())
                                      //                  .face_index_map(get(CGAL::face_index, sm))
                                                                );
  
 std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
 BOOST_FOREACH(face_descriptor f , faces(sm)){
   std::cout  << &*f << " in connected component " << fccmap[f] << std::endl;
  }
 
 CGAL::Polygon_mesh_processing::keep_largest_connected_components(sm,2);

 std::cout << "mesh:\n" << sm << std::endl;
  return 0;
}
