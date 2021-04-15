#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Surface_mesh_deformation.h>

#include <fstream>


typedef CGAL::Simple_cartesian<double>                                   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator        vertex_iterator;

typedef CGAL::Surface_mesh_deformation<Polyhedron> Surface_mesh_deformation;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< "Cannot open  data/plane.off" << std::endl;
    return 1;
  }

  // Init the indices of the halfedges and the vertices.
  set_halfedgeds_items_id(mesh);

  // Create a deformation object
  Surface_mesh_deformation deform_mesh(mesh);

  // Definition of the region of interest (use the whole mesh)
  vertex_iterator vb,ve;
  boost::tie(vb, ve) = vertices(mesh);
  deform_mesh.insert_roi_vertices(vb, ve);

  // Select two control vertices ...
  vertex_descriptor control_1 = *std::next(vb, 213);
  vertex_descriptor control_2 = *std::next(vb, 157);

  // ... and insert them
  deform_mesh.insert_control_vertex(control_1);
  deform_mesh.insert_control_vertex(control_2);

  // The definition of the ROI and the control vertices is done, call preprocess
  bool is_matrix_factorization_OK = deform_mesh.preprocess();
  if(!is_matrix_factorization_OK){
    std::cerr << "Error in preprocessing, check documentation of preprocess()" << std::endl;
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

  // Save the deformed mesh into a file
  std::ofstream output("deform_1.off");
  output << mesh;
  output.close();

  // Add another control vertex which requires another call to preprocess
  vertex_descriptor control_3 = *std::next(vb, 92);
  deform_mesh.insert_control_vertex(control_3);

  // The prepocessing step is again needed
  if(!deform_mesh.preprocess()){
    std::cerr << "Error in preprocessing, check documentation of preprocess()" << std::endl;
    return 1;
  }

  // Deform the mesh
  Surface_mesh_deformation::Point constrained_pos_3(0.55, 0.30, -0.70);
  deform_mesh.set_target_position(control_3, constrained_pos_3);

  deform_mesh.deform(15, 0.0);

  output.open("deform_2.off");
  output << mesh;
}
