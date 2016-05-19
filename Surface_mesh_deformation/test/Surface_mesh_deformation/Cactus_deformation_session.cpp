#include "Surface_mesh_deformation_test_commons.h"
#include <algorithm>
#include <sstream>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_deformation.h>

#include <CGAL/Timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;

typedef CGAL::Surface_mesh_deformation<Polyhedron, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP>  Deform_mesh_arap;
typedef CGAL::Surface_mesh_deformation<Polyhedron, CGAL::Default, CGAL::Default, CGAL::SPOKES_AND_RIMS> Deform_mesh_spoke;

const double squared_threshold = 0.001; // alert if average difs between precomputed and deformed mesh models is above threshold

void compare_mesh(const Polyhedron& mesh_1, const Polyhedron& mesh_2)
{
  Polyhedron::Vertex_const_iterator it_1 = mesh_1.vertices_begin();
  Polyhedron::Vertex_const_iterator it_2 = mesh_2.vertices_begin();
  Kernel::Vector_3 total_dif(0,0,0);
  for( ; it_1 != mesh_1.vertices_end(); ++it_1 , ++it_2)
  {
    total_dif = total_dif + (it_1->point() - it_2->point());
  }
  double n = static_cast<double>(mesh_1.size_of_vertices());
  double average_mesh_dif = total_dif.squared_length() / n / n;



  std::cerr << "Average mesh difference: " << average_mesh_dif << std::endl;

  assert( average_mesh_dif < squared_threshold);
}

// read deformation session saved as a handle differences
template<class DeformMesh, class InputIterator>
void read_handle_difs_and_deform(DeformMesh& deform_mesh, InputIterator begin, InputIterator end)
{
  typedef CGAL::Simple_cartesian<double>::Vector_3 Vector;

  if(!deform_mesh.preprocess()) {
    std::cerr << "Error: preprocess() failed!" << std::endl;
    assert(false);
  }

  std::ifstream dif_stream("data/cactus_handle_differences.txt");
  std::vector<Vector> dif_vector;
  double x, y, z;
  while(dif_stream >> x >> y >> z)
  { dif_vector.push_back(Vector(x, y, z)); }

  CGAL::Timer timer;

  //the original behavior of translate was to overwrite the previous
  //translation. Now that it is cumulative, we need to substract the
  //previous translation vector to mimic the overwrite
  Vector previous(0,0,0);
  for(std::size_t i = 0; i < dif_vector.size(); ++i)
  {
    timer.start();
    deform_mesh.translate(begin, end, dif_vector[i]-previous);
    deform_mesh.deform();
    timer.stop();
    previous=dif_vector[i];

    // read pre-deformed cactus
    std::stringstream predeformed_cactus_file;
    predeformed_cactus_file << "data/cactus_deformed/cactus_deformed_" << i << ".off";
    Polyhedron predeformed_cactus;

    read_to_polyhedron(predeformed_cactus_file.str().c_str(), predeformed_cactus);
    compare_mesh(predeformed_cactus, deform_mesh.halfedge_graph());

    // for saving deformation
    //std::ofstream(predeformed_cactus_file) << deform_mesh.halfedge_graph();
    //std::cerr << predeformed_cactus_file << std::endl;
  }
  std::cerr << "Deformation performance (with default number_of_iteration and tolerance) " << std::endl
    << dif_vector.size() << " translation: " << timer.time() << std::endl;
}

int main()
{
  Polyhedron mesh_1;
  read_to_polyhedron("data/cactus.off", mesh_1);
  Polyhedron mesh_2 = mesh_1;

  init_indices(mesh_1);
  init_indices(mesh_2);

  Deform_mesh_arap deform_mesh_arap(mesh_1);
  Deform_mesh_spoke deform_mesh_spoke(mesh_2);
  // For original arap
  std::vector<vertex_descriptor> hg_1 =
    read_rois(deform_mesh_arap, "data/cactus_roi.txt", "data/cactus_handle.txt");

  read_handle_difs_and_deform(deform_mesh_arap, hg_1.begin(), hg_1.end());
  std::cerr << "ORIGINAL ARAP Success!" << std::endl;
  // For spokes rims
  std::vector<vertex_descriptor> hg_2 =
    read_rois(deform_mesh_spoke, "data/cactus_roi.txt", "data/cactus_handle.txt");

  read_handle_difs_and_deform(deform_mesh_spoke, hg_2.begin(), hg_2.end());
  std::cerr << "SPOKES AND RIMS ARAP Success!" << std::endl;
  std::cerr << "All done!" << std::endl;
}

