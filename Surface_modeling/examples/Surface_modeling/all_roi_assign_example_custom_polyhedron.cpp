#include <fstream>
#include <map>
#include <cmath>
#include <boost/property_map/property_map.hpp>

struct Custom_point_3{
  //needed by File_scanner_OFF
  struct R{
    typedef double RT;
  };

  double coords[3];
  Custom_point_3(){}
  Custom_point_3(double x, double y, double z)
  { coords[0]=x; coords[1]=y; coords[2]=z; }
  Custom_point_3(double x, double y, double z, double w)
  { coords[0]=x/w; coords[1]=y/w; coords[2]=z/w; }

  double x() const {return coords[0];}
  double y() const {return coords[1];}
  double z() const {return coords[2];}

  double& operator[](int i)       { return coords[i]; }
  double  operator[](int i) const { return coords[i]; }

  friend std::ostream& operator<<(std::ostream& out, const Custom_point_3& p)
  {
    out << p.x() << " " << p.y() << " " << p.z();
    return out;
  }

  friend std::istream& operator<<(std::istream& in, Custom_point_3& p)
  {
    in >> p.coords[0] >> p.coords[1] >> p.coords[2];
    return in;
  }
};

//needed by CGAL/boost/graph/properties_Polyhedron_3.h
namespace CGAL{
  double squared_distance(const Custom_point_3& p1, const Custom_point_3& p2){
    return std::pow(p1.x()-p2.x(),2) + std::pow(p1.y()-p2.y(), 2) + std::pow(p1.z()-p2.z(), 2);
  }
}

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// Halfedge adapters for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <CGAL/Deform_mesh.h>

struct Custom_traits{
  typedef Custom_point_3 Point_3;
  struct Plane_3{};
};

typedef CGAL::Polyhedron_3<Custom_traits>       Polyhedron;

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

  // index maps should contain unique indices with 0 offset
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

  deform_mesh.insert_roi(vb, ve); // insert whole mesh as roi

  vertex_descriptor handle_1 = *boost::next(vb, 213);
  vertex_descriptor handle_2 = *boost::next(vb, 157);

  deform_mesh.insert_handle(handle_1); // insert handles
  deform_mesh.insert_handle(handle_2);

  // insertion of roi and handles completed, call preprocess
  bool is_matrix_factorization_OK = deform_mesh.preprocess();
  if(!is_matrix_factorization_OK){ 
    std::cerr << "Check documentation of preprocess()" << std::endl; 
    return 1;
  }

//// DEFORM SECTION ////
  // now use assign() to provide constained positions of handles
  Deform_mesh::Point constrained_pos_1(-0.35, 0.40, 0.60); // target position of handle_1
  deform_mesh.assign(handle_1, constrained_pos_1);
  // note that we only assign a constraint for handle_1, other handles will be constrained to last assigned positions

  // deform the mesh, now positions of vertices of 'mesh' will be changed
  deform_mesh.deform();
  deform_mesh.deform(); // you can call deform multiple times if you like

  Deform_mesh::Point constrained_pos_2(0.55, -0.30, 0.70);
  deform_mesh.assign(handle_2, constrained_pos_2);
  // note that handle_1 will be still constrained to constrained_pos_1,

  deform_mesh.deform(10, 0.0); // deform(unsigned int iterations, double tolerance) can be called with instant parameters
  // this time iterate 10 times, and do not use energy based termination

  std::ofstream output("data/deform_1.off");
  output << mesh; // save deformed mesh
  output.close();

  // want to add another handle
//// PREPROCESS SECTION AGAIN////
  vertex_descriptor handle_3 = *boost::next(vb, 92);
  deform_mesh.insert_handle(handle_3); // now I need to prepocess again

  if(!deform_mesh.preprocess()) {
    std::cerr << "Check documentation of preprocess()" << std::endl; 
    return 1;
  }

//// DEFORM SECTION AGAIN////
  Deform_mesh::Point constrained_pos_3(0.55, 0.30, -0.70);
  deform_mesh.assign(handle_3, constrained_pos_3);

  deform_mesh.deform(15, 0.0);

  output.open("data/deform_2.off");
  output << mesh;
}