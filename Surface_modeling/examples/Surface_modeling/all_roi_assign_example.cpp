#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// HalfedgeGraph adapters for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <CGAL/Deform_mesh.h>

#include <fstream>
#include <map>
#include <boost/property_map/property_map.hpp>
#include <boost/utility.hpp>

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator      vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor      edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator        edge_iterator;

typedef std::map<vertex_descriptor, std::size_t>   Internal_vertex_map;
typedef std::map<edge_descriptor, std::size_t>     Internal_edge_map;

typedef boost::associative_property_map<Internal_vertex_map>   Vertex_index_map;
typedef boost::associative_property_map<Internal_edge_map>     Edge_index_map;

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map> Deform_mesh;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< "Cannot open  data/plane.off" << std::endl;
    return 1;
  }
  // index maps must contain an index unique per vertex starting from 0
  // to the total number of vertices
  Internal_vertex_map internal_vertex_index_map;
  Vertex_index_map vertex_index_map(internal_vertex_index_map);
  vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(boost::tie(vb, ve) = boost::vertices(mesh); vb != ve; ++vb, ++counter) {
    put(vertex_index_map, *vb, counter);
  }

  Internal_edge_map internal_edge_index_map;
  Edge_index_map edge_index_map(internal_edge_index_map);
  counter = 0;
  edge_iterator eb, ee;
  for(boost::tie(eb, ee) = boost::edges(mesh); eb != ee; ++eb, ++counter) {
    put(edge_index_map, *eb, counter);
  }
//// PREPROCESS SECTION ////
  Deform_mesh deform_mesh(mesh, vertex_index_map, edge_index_map);

  // insert region of interest
  boost::tie(vb, ve) = boost::vertices(mesh);

  deform_mesh.insert_roi_vertices(vb, ve); // insert whole mesh as ROI

  // insert controls
  vertex_descriptor control_1 = *boost::next(vb, 213);
  vertex_descriptor control_2 = *boost::next(vb, 157);

  deform_mesh.insert_control_vertex(control_1); // insert controls
  deform_mesh.insert_control_vertex(control_2);

  // insertion of ROI and controls completed, call preprocess
  bool is_matrix_factorization_OK = deform_mesh.preprocess();
  if(!is_matrix_factorization_OK){ 
    std::cerr << "Check documentation of preprocess()" << std::endl; 
    return 1;
  }

//// DEFORM SECTION ////
  // now use set_target_position() to provide constained positions of controls
  Deform_mesh::Point constrained_pos_1(-0.35, 0.40, 0.60); // target position of control_1
  deform_mesh.set_target_position(control_1, constrained_pos_1);
  // note that we only assign a constraint for control_1, other controls will be constrained to last assigned positions

  // deform the mesh, now positions of vertices of 'mesh' will be changed
  deform_mesh.deform();
  deform_mesh.deform(); // you can call deform multiple times if you like

  Deform_mesh::Point constrained_pos_2(0.55, -0.30, 0.70);
  deform_mesh.set_target_position(control_2, constrained_pos_2);
  // note that control_1 will be still constrained to constrained_pos_1,

  deform_mesh.deform(10, 0.0); // deform(unsigned int iterations, double tolerance) can be called with instant parameters
  // this time iterate 10 times, and do not use energy based termination

  std::ofstream output("deform_1.off");
  output << mesh; // save deformed mesh
  output.close();

  // want to add another control
//// PREPROCESS SECTION AGAIN////
  vertex_descriptor control_3 = *boost::next(vb, 92);
  deform_mesh.insert_control_vertex(control_3); // now I need to prepocess again

  if(!deform_mesh.preprocess()){ 
    std::cerr << "Check documentation of preprocess()" << std::endl;
    return 1;
  }

//// DEFORM SECTION AGAIN////
  Deform_mesh::Point constrained_pos_3(0.55, 0.30, -0.70);
  deform_mesh.set_target_position(control_3, constrained_pos_3);

  deform_mesh.deform(15, 0.0);

  output.open("deform_2.off");
  output << mesh;
}