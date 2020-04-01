#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

// HalfedgeGraph adapters
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <CGAL/Surface_mesh_deformation.h>

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/>                     Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Mesh>::vertex_iterator        vertex_iterator;
typedef boost::graph_traits<Mesh>::halfedge_iterator    halfedge_iterator;

typedef CGAL::Surface_mesh_deformation<Mesh> Surface_mesh_deformation;

int main()
{
  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, "data/plane.off");

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
  OpenMesh::IO::write_mesh(mesh,"deform_1.off");

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

  OpenMesh::IO::write_mesh(mesh,"deform_2.off");
}
