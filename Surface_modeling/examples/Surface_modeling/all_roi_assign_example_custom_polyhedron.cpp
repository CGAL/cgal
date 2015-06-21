#include <fstream>
#include <map>
#include <cmath>
#include <CGAL/property_map.h>

struct Custom_point_3{
  // Required by File_scanner_OFF
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

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// Halfedge adapters for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <CGAL/Surface_mesh_deformation.h>

struct Custom_traits{
  typedef Custom_point_3 Point_3;
  struct Plane_3{};
};

typedef CGAL::Polyhedron_3<Custom_traits>       Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator      vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_iterator    halfedge_iterator;

typedef std::map<vertex_descriptor, std::size_t>   Internal_vertex_map;
typedef std::map<halfedge_descriptor, std::size_t>     Internal_hedge_map;

typedef boost::associative_property_map<Internal_vertex_map>   Vertex_index_map;
typedef boost::associative_property_map<Internal_hedge_map>     Hedge_index_map;

typedef CGAL::Surface_mesh_deformation<Polyhedron, Vertex_index_map, Hedge_index_map> Surface_mesh_deformation;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< "Cannot open  data/plane.off" << std::endl;
    return 1;
  }

  // Index maps must contain an index unique per vertex starting from 0
  // to the total number of vertices
  Internal_vertex_map internal_vertex_index_map;
  Vertex_index_map vertex_index_map(internal_vertex_index_map);
  vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(boost::tie(vb, ve) = vertices(mesh); vb != ve; ++vb, ++counter) {
    put(vertex_index_map, *vb, counter);
  }

  Internal_hedge_map internal_hedge_index_map;
  Hedge_index_map hedge_index_map(internal_hedge_index_map);
  counter = 0;
  halfedge_iterator eb, ee;
  for(boost::tie(eb, ee) = halfedges(mesh); eb != ee; ++eb, ++counter) {
    put(hedge_index_map, *eb, counter);
  }

  Surface_mesh_deformation deform_mesh(mesh, vertex_index_map, hedge_index_map);

  // Insert the whole mesh as region of interest
  boost::tie(vb, ve) = vertices(mesh);
  deform_mesh.insert_roi_vertices(vb, ve);

  // Insert two control vertices
  vertex_descriptor control_1 = *CGAL::cpp11::next(vb, 213);
  vertex_descriptor control_2 = *CGAL::cpp11::next(vb, 157);
  deform_mesh.insert_control_vertex(control_1);
  deform_mesh.insert_control_vertex(control_2);

  // The definition of the ROI and the control vertices is done, call preprocess
  bool is_matrix_factorization_OK = deform_mesh.preprocess();
  if(!is_matrix_factorization_OK){
    std::cerr << "Check documentation of preprocess()" << std::endl;
    return 1;
  }

  // Use set_target_position() to set the constained position
  // of control_1. control_2 remains at the last assigned positions
  Surface_mesh_deformation::Point constrained_pos_1(-0.35, 0.40, 0.60);
  deform_mesh.set_target_position(control_1, constrained_pos_1);

  // Deform the mesh, the positions of vertices of 'mesh' are updated
  deform_mesh.deform();
  // The function deform() can be called several times if the convergence has not been reached yet
  deform_mesh.deform();

  // Set the constained position of control_2
  Surface_mesh_deformation::Point constrained_pos_2(0.55, -0.30, 0.70);
  deform_mesh.set_target_position(control_2, constrained_pos_2);


  // Call the function deform() with one-time parameters:
  // iterate 10 times and do not use energy based termination criterion
  deform_mesh.deform(10, 0.0);

  std::ofstream output("deform_1.off");
  output << mesh; // save deformed mesh
  output.close();

  // Add another control vertex
  vertex_descriptor control_3 = *CGAL::cpp11::next(vb, 92);
  deform_mesh.insert_control_vertex(control_3);

  // The prepocessing step is again needed
  if(!deform_mesh.preprocess()) {
    std::cerr << "Check documentation of preprocess()" << std::endl;
    return 1;
  }

  Surface_mesh_deformation::Point constrained_pos_3(0.55, 0.30, -0.70);
  deform_mesh.set_target_position(control_3, constrained_pos_3);

  deform_mesh.deform(15, 0.0);

  output.open("deform_2.off");
  output << mesh;
}
