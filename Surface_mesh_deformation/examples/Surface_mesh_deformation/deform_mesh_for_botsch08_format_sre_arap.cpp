#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh_deformation.h>

#include <fstream>


typedef CGAL::Simple_cartesian<double>                                   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator        vertex_iterator;

typedef CGAL::Surface_mesh_deformation<Polyhedron,CGAL::Default, CGAL::Default, CGAL::SRE_ARAP> Surface_mesh_deformation;

int main(int argc,char** argv)
{
  std::string off_name=CGAL::data_file_path("meshes/cactus.off"),
              sel_name="data/cactus.sel",
              def_name="data/cactus.def";
  if  ( argc!=4){
    std::cerr <<"Usage " << argv[0] << " input.off input.sel input.def\n";
    std::cerr <<"Using default values: " << off_name << " " << sel_name << " " << def_name << "\n";
  }
  else
  {
    off_name=argv[1];
    sel_name=argv[2];
    def_name=argv[3];
  }

  Polyhedron mesh;
  std::ifstream input(off_name);

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< off_name << " is not a valid off file" << std::endl;
    return 1;
  }
  input.close();

  // Init the indices of the halfedges and the vertices.
  set_halfedgeds_items_id(mesh);

  // Create a deformation object
  Surface_mesh_deformation deform_mesh(mesh);

  // Changing alpha value
  deform_mesh.set_sre_arap_alpha(0.02);


  // Definition of the region of interest (use the whole mesh)
  vertex_iterator vb,ve;
  boost::tie(vb, ve) = vertices(mesh);

  //the selection is set by a file
  input.open(sel_name);
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
  input.open(def_name);
  double m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, hw, sink;
  getline(input, line); // skip first comment line
  input >> m00 >> m01 >> m02 >> m03;
  input >> m10 >> m11 >> m12 >> m13;
  input >> m20 >> m21 >> m22 >> m23;
  input >> sink >> sink >> sink >> hw;

  std::cout << "Setting target positions\n";
  Kernel::Aff_transformation_3 aff(m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23);
  for(vertex_descriptor vd : control_vertices)
  {
    Surface_mesh_deformation::Point pos = vd->point().transform(aff);
    deform_mesh.set_target_position(vd, pos);
  }

  // Call the function deform() with one-time parameters:
  std::cout << "Deforming the mesh\n";
  deform_mesh.deform(1000, 1e-4);

  // Save the deformed mesh into a file
  std::ofstream output("deform_res.off");
  output << mesh;
  output.close();
}
