#include "Surface_mesh_deformation_test_commons.h"
#include <algorithm>
#include <sstream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

// HalfedgeGraph adapters
#define CGAL_USE_OM_POINTS // use OpenMesh point type
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <CGAL/Surface_mesh_deformation.h>

#include <CGAL/Timer.h>

struct DoubleTraits : public OpenMesh::DefaultTraits
{
  typedef OpenMesh::Vec3d Point;
  typedef OpenMesh::Vec3d Normal;
};


typedef OpenMesh::PolyMesh_ArrayKernelT<DoubleTraits>               Mesh;
typedef Mesh::Point                                                 Point;
typedef boost::graph_traits<Mesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Mesh>::vertex_iterator        vertex_iterator;

typedef CGAL::Surface_mesh_deformation<Mesh, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP>  Deform_mesh_arap;
typedef CGAL::Surface_mesh_deformation<Mesh, CGAL::Default, CGAL::Default, CGAL::SPOKES_AND_RIMS> Deform_mesh_spoke;
const double squared_threshold = 0.001; // alert if average difs between precomputed and deformed mesh models is above threshold

double squared_length(const Point& p)
{
  std::cerr << typeid(p[0]).name() << std::endl;
  return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
}

void compare_mesh(const Mesh& mesh_1, const Mesh& mesh_2)
{
  vertex_iterator it_1,end_1, it_2, end_2;
  CGAL::cpp11::tie(it_1, end_1) = vertices(mesh_1);
  CGAL::cpp11::tie(it_2, end_2) = vertices(mesh_2);
  boost::property_map<Mesh, boost::vertex_point_t>::type
    ppmap_1 = get(boost::vertex_point, mesh_1), ppmap_2 = get(boost::vertex_point, mesh_2);
  Point total_dif(0,0,0);
  for( ; it_1 != end_1; ++it_1 , ++it_2)
  {
    total_dif = total_dif + (get(ppmap_1, *it_1) - get(ppmap_2, *it_2));
  }
  const double n = static_cast<double>(num_vertices(mesh_1));
  double average_mesh_dif = squared_length(total_dif) / n / n;

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
    Mesh predeformed_cactus;

    OpenMesh::IO::read_mesh(predeformed_cactus, predeformed_cactus_file.str().c_str());
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
  Mesh mesh_1;
  OpenMesh::IO::read_mesh(mesh_1, "data/cactus.off");
  Mesh mesh_2 = mesh_1;

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

