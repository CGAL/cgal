#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// HalfedgeGraph adapters for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <boost/foreach.hpp>
// #define CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SR_ARAP
#include <CGAL/Surface_mesh_deformation.h>

#include <fstream>


typedef CGAL::Simple_cartesian<double>                                   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator        vertex_iterator;

typedef CGAL::Surface_mesh_deformation<Polyhedron,CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP> Surface_mesh_deformation;

int main(int argc,char** argv)
{
  if  ( argc!=4){
    std::cerr <<"Usage " << argv[0] << " input.off input.sel input.def\n";
    return 1;
  }
  Polyhedron mesh;
  std::ifstream input(argv[1]);

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< argv[1] << " is not a valid off file" << std::endl;
    return 1;
  }
  input.close();

  // Init the indices of the halfedges and the vertices.
  set_halfedgeds_items_id(mesh);

  // Create a deformation object
  Surface_mesh_deformation deform_mesh(mesh);

  // Definition of the region of interest (use the whole mesh)
  vertex_iterator vb,ve;
  boost::tie(vb, ve) = boost::vertices(mesh);

  //the selection is set by a file
  input.open(argv[2]);
  std::string line;
  std::vector<vertex_descriptor> control_vertices;
  while(getline(input, line))
  {
    if (line[0]=='#') continue;
    if (line[0]=='1') deform_mesh.insert_roi_vertex(*vb);
    if (line[0]=='2') {
      deform_mesh.insert_control_vertex(*vb);
      control_vertices.push_back(*vb);
    }
    ++vb;
    if (vb==ve) break;
  }
  input.close();

  std::cout << "Using " << control_vertices.size() << " control vertices\n";
  // The definition of the ROI and the control vertices is done, call preprocess
  bool is_matrix_factorization_OK = deform_mesh.preprocess();
  if(!is_matrix_factorization_OK){
    std::cerr << "Error in preprocessing, check documentation of preprocess()" << std::endl;
    return 1;
  }

  //define the transformation
  input.open(argv[3]);
  double m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, hw, sink;
  getline(input, line); // skip first comment line
  input >> m00 >> m01 >> m02 >> m03;
  input >> m10 >> m11 >> m12 >> m13;
  input >> m20 >> m21 >> m22 >> m23;
  input >> sink >> sink >> sink >> hw;

  Kernel::Aff_transformation_3 aff(m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23);
  BOOST_FOREACH(vertex_descriptor vd, control_vertices)
  {
    Surface_mesh_deformation::Point pos = vd->point().transform(aff);
    deform_mesh.set_target_position(vd, pos);
  }

  // Call the function deform() with one-time parameters:
  deform_mesh.deform(1000, 1e-4);

  // Save the deformed mesh into a file
  std::ofstream output("deform_res.off");
  output << mesh;
  output.close();
}
