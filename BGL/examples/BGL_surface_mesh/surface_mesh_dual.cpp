#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>

#include <iostream>
#include <fstream>

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>

typedef CGAL::Simple_cartesian<double>             Kernel;
typedef Kernel::Point_3                            Point;
typedef CGAL::Surface_mesh<Point>                  Mesh;
typedef CGAL::Dual<Mesh>                           Dual;
typedef boost::graph_traits<Dual>::edge_descriptor edge_descriptor;

template <typename G>
struct noborder {
  noborder() : g(NULL) {} // default-constructor required by filtered_graph
  noborder(G& g) : g(&g) {}

  bool operator()(const edge_descriptor& e) const
  { return !is_border(e,*g); }

  G* g;
};


// A dual border edge has a null_face as the source or target "vertex"
// BGL algorithms won't like that, so we remove border edges through a
// boost::filtered_graph.
typedef boost::filtered_graph<Dual, noborder<Mesh> >   FiniteDual;
typedef boost::graph_traits<Mesh>::vertex_descriptor   vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor     face_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;

int main(int argc, char* argv[])
{
  Mesh primal;
  const char* filename = (argc > 1) ? argv[1] : "data/prim.off";
  std::ifstream in(filename);
  if(!(in >> primal)) {
    std::cerr << "Error reading polyhedron from file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  Dual dual(primal);
  FiniteDual finite_dual(dual,noborder<Mesh>(primal));

  std::cout << "dual has " << num_vertices(dual) << " vertices" << std::endl;

  std::cout << "The vertices of dual are faces in primal"<< std::endl;
  for(boost::graph_traits<Dual>::vertex_descriptor dvd : vertices(dual)) {
    std::cout << dvd << std::endl;
  }

  std::cout << "The edges in primal and dual with source and target" << std::endl;
  for(edge_descriptor e : edges(dual)) {
   std::cout << e << " in primal:  " << source(e,primal)      << " -- " << target(e,primal)       << "   "
             <<      " in dual  :  " << source(e,finite_dual) << " -- " << target(e,finite_dual)  << std::endl;
  }


 std::cout << "edges of the finite dual graph" << std::endl;
 for(boost::graph_traits<FiniteDual>::edge_descriptor e : CGAL::make_range(edges(finite_dual))) {
   std::cout << e << "  " << source(e,primal) << " " << source(e,finite_dual)  << std::endl;
 }

 // the storage of a property map is in primal
 Mesh::Property_map<face_descriptor,int> fccmap;
 fccmap = primal.add_property_map<face_descriptor,int>("f:CC").first;
 int num = connected_components(finite_dual, fccmap);

 std::cout << "The graph has " << num << " connected components (face connectivity)" << std::endl;
 for(face_descriptor f : faces(primal)) {
   std::cout << f << " in connected component " << fccmap[f] << std::endl;
 }

 Mesh::Property_map<vertex_descriptor,int> vccmap;
 vccmap = primal.add_property_map<vertex_descriptor,int>("v:CC").first;
 num = connected_components(primal, vccmap);

 std::cout << "The graph has " << num << " connected components (edge connectvity)" << std::endl;
 for(vertex_descriptor v : vertices(primal)) {
   std::cout << v << " in connected component " << vccmap[v] << std::endl;
 }
  return 0;
}
