#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// HalfedgeGraph adapters for Polyhedron_3
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <CGAL/Deform_mesh.h>

#include <fstream>


typedef CGAL::Simple_cartesian<double>                                   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator        vertex_iterator;

typedef CGAL::Deform_mesh<Polyhedron> Deform_mesh;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< "Cannot open  data/plane.off" << std::endl;
    return 1;
  }

  //Init the indices of the halfedges and the vertices.
  set_halfedgeds_items_id(mesh);

  //Create deformation object
  Deform_mesh deform_mesh(mesh);

  // Definition of the region of interest (use the whole mesh)
  vertex_iterator vb,ve;
  boost::tie(vb, ve) = boost::vertices(mesh);
  deform_mesh.insert_roi_vertices(vb, ve);

  // Select two control vertices ...
  vertex_descriptor control_1 = *boost::next(vb, 213);
  vertex_descriptor control_2 = *boost::next(vb, 157);

  // ... and insert them
  deform_mesh.insert_control_vertex(control_1);
  deform_mesh.insert_control_vertex(control_2);

  // The definition of the      ROI and the control vertices is done, call preprocess
  bool is_matrix_factorization_OK = deform_mesh.preprocess();
  if(!is_matrix_factorization_OK){
    std::cerr << "Error in preprocessing, check documentation of preprocess()" << std::endl;
    return 1;
  }

  // now use set_target_position() to provide constained positions of control vertices
  Deform_mesh::Point constrained_pos_1(-0.35, 0.40, 0.60); // target position of control_1
  deform_mesh.set_target_position(control_1, constrained_pos_1);
  // note that we only assign a constraint for control_1, other control vertices will be constrained to last assigned positions

  // deform the mesh, now the positions of vertices of 'mesh' will be changed
  deform_mesh.deform();
  // deform can be called several times if the convergence has not been reached yet
  deform_mesh.deform();

  Deform_mesh::Point constrained_pos_2(0.55, -0.30, 0.70);
  deform_mesh.set_target_position(control_2, constrained_pos_2);
  // note that control_1 will be still constrained to constrained_pos_1,

  deform_mesh.deform(10, 0.0); // deform(unsigned int iterations, double tolerance) can be called with instant parameters
  // this time iterate 10 times, and do not use energy based termination

  std::ofstream output("deform_1.off");
  output << mesh; // save deformed mesh
  output.close();

  // We add another control vertex which require another call to preprocess
  vertex_descriptor control_3 = *boost::next(vb, 92);
  deform_mesh.insert_control_vertex(control_3); // now I need to prepocess again

  if(!deform_mesh.preprocess()){
    std::cerr << "Error in preprocessing, check documentation of preprocess()" << std::endl;
    return 1;
  }

  // Deform the mesh
  Deform_mesh::Point constrained_pos_3(0.55, 0.30, -0.70);
  deform_mesh.set_target_position(control_3, constrained_pos_3);

  deform_mesh.deform(15, 0.0);

  output.open("deform_2.off");
  output << mesh;
}
