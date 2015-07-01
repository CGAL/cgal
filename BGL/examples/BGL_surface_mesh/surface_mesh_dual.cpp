#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>
#include <iostream>
#include <fstream>

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                           Mesh;
typedef CGAL::Dual<Mesh> Dual;
typedef boost::graph_traits<Dual>::edge_descriptor edge_descriptor;

template <typename G>
struct noborder {
  noborder() { }

  noborder(G & g) : g(&g) { }

  bool operator()(const edge_descriptor& e) const {
    return ! is_border(e,*g);
  }

  G* g;
};


// A dual border edge has a null_face as source or target "vertex"
// BGL algorithms won't like that
typedef boost::filtered_graph<Dual, noborder<Mesh> > FiniteDual;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

int main(int, char* argv[]) 
{
  Mesh primal;
  std::ifstream in(argv[1]);
  in >> primal;

  Dual dual(primal);
  FiniteDual finite_dual(dual,noborder<Mesh>(primal));

  std::cout << "dual has " << num_vertices(dual) << " vertices" << std::endl;
  
  std::cout << "The vertices of dual are faces in primal"<< std::endl;
  BOOST_FOREACH(boost::graph_traits<Dual>::vertex_descriptor dvd , vertices(dual)){
    std::cout  << dvd << std::endl;
  }
  std::cerr << "The halfedges in primal and dual with source and target"<< std::endl;
 BOOST_FOREACH(halfedge_descriptor h , halfedges(dual)){
   std::cout  << h << " in primal:  " << source(h,primal) << " -- " << target(h,primal) << "   "
              <<      " in dual  :  " << source(h,finite_dual)<< " -- " << target(h,finite_dual)  << std::endl;
  }
  

 std::cout << "edges of the finite dual graph" << std::endl;
 BOOST_FOREACH(boost::graph_traits<FiniteDual>::edge_descriptor e , edges(finite_dual)){
   std::cout  << e << "  " <<  source(e,primal) << " " << source(e,finite_dual)  << std::endl;
  }

 // the storage of a property map is in primal
 Mesh::Property_map<face_descriptor,int> fccmap;
 fccmap = primal.add_property_map<face_descriptor,int>("f:CC").first; 
 int num = connected_components(finite_dual, fccmap);
 
 std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
 BOOST_FOREACH(face_descriptor f , faces(primal)){
   std::cout  << f << " in connected component " << fccmap[f] << std::endl;
  }

 Mesh::Property_map<vertex_descriptor,int> vccmap;
 vccmap = primal.add_property_map<vertex_descriptor,int>("v:CC").first; 
 num = connected_components(primal, vccmap);
 
 std::cerr << "The graph has " << num << " connected components (edge connectvity)" << std::endl;
 BOOST_FOREACH(vertex_descriptor v , vertices(primal)){
   std::cout  << v << " in connected component " << vccmap[v] << std::endl;
  }
  return 0;
}
