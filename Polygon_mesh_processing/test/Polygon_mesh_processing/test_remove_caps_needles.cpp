#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>

#include <CGAL/Polyhedral_envelope.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <fstream>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor     face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Polyhedral_envelope<K> Envelope;

void general_test(std::string filename)
{
  std::cout << "Removing caps/needles no extra parameters\n";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "  Input mesh has " << edges(mesh).size() << " edges\n";
  if (PMP::does_self_intersect(mesh))
    std::cout << "  Input mesh has self-intersections\n";

  PMP::remove_almost_degenerate_faces(mesh,
                                      params::cap_threshold(std::cos(160. / 180 * CGAL_PI))
                                             .needle_threshold(4)
                                             .collapse_length_threshold(0.14));


  CGAL::IO::write_polygon_mesh("cleaned_mesh.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "  Output mesh has " << edges(mesh).size() << " edges\n";
  if (PMP::does_self_intersect(mesh))
    std::cout << "  Output mesh has self-intersections\n";
}

void test_with_envelope(std::string filename, double eps)
{
  std::cout << "Removing caps/needles with envelope, epsilon = " << eps << "\n";
  std::ifstream input(filename);

  Mesh mesh, bk;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    exit(EXIT_FAILURE);
  }

  if (PMP::does_self_intersect(mesh))
    std::cout << "  Input mesh has self-intersections\n";

  // first a test where nothing should be done
  struct No_modification_allowed
  {
    bool operator()(K::Point_3, K::Point_3, K::Point_3) const { return false; }
  };
  No_modification_allowed no_modif;
  const std::size_t nbv = vertices(mesh).size();
  PMP::remove_almost_degenerate_faces(mesh,
                                      params::cap_threshold(std::cos(160. / 180 * CGAL_PI))
                                             .needle_threshold(4)
                                             .collapse_length_threshold(0.14)
                                             .filter(no_modif));
  assert(nbv == vertices(mesh).size());

  // now the real test with a fixed envelope
  std::cout << "Using fixed envelope\n";
  std::cout << "  Input mesh has " << edges(mesh).size() << " edges\n";
  bk=mesh;
  Envelope envelope(mesh, eps);
  PMP::remove_almost_degenerate_faces(mesh,
                                      params::cap_threshold(std::cos(160. / 180 * CGAL_PI))
                                             .needle_threshold(4)
                                             .collapse_length_threshold(0.14)
                                             .filter(std::ref(envelope)));

  CGAL::IO::write_polygon_mesh("cleaned_mesh_with_envelope.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "  Output mesh has " << edges(mesh).size() << " edges\n";
  if (PMP::does_self_intersect(mesh))
    std::cout << "  Output mesh has self-intersections\n";


  // now the real test with an iterative envelope
  mesh=bk;
  std::cout << "Using iteratively created fixed envelope\n";
  std::cout << "  Input mesh has " << edges(mesh).size() << " edges\n";
  auto create_envelope = [&](const std::vector<Mesh::Face_index>& frange) -> Envelope
                          {
                            return Envelope(frange, mesh, eps);
                          };
  std::function<Envelope(const std::vector<Mesh::Face_index>&)> filter(create_envelope);
  PMP::remove_almost_degenerate_faces(mesh,
                                      params::cap_threshold(std::cos(160. / 180 * CGAL_PI))
                                             .needle_threshold(4)
                                             .collapse_length_threshold(0.14)
                                             .filter(filter));

  CGAL::IO::write_polygon_mesh("cleaned_mesh_with_iterative_envelope.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "  Output mesh has " << edges(mesh).size() << " edges\n";
  if (PMP::does_self_intersect(mesh))
    std::cout << "  Output mesh has self-intersections\n";
}

bool same_meshes(const Mesh& m1, const Mesh& m2)
{
  std::size_t c=0, m1_only=0, m2_only=0;
  PMP::match_faces(m1, m2, CGAL::Counting_output_iterator(&c)
                         , CGAL::Counting_output_iterator(&m1_only)
                         , CGAL::Counting_output_iterator(&m2_only));
  return m1_only==0 && m2_only==0;
}
void test_parameters_on_pig(std::string filename)
{
  std::ifstream input(filename);

  Mesh mesh, bk;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    exit(EXIT_FAILURE);
  }

  bk=mesh;

  PMP::remove_almost_degenerate_faces(mesh,
                                      params::cap_threshold(std::cos(160. / 180 * CGAL_PI))
                                             .needle_threshold(4));
  assert(vertices(mesh).size()!=vertices(bk).size());

  mesh=bk;
  PMP::remove_almost_degenerate_faces(mesh,
                                      params::cap_threshold(std::cos(160. / 180 * CGAL_PI))
                                             .needle_threshold(4)
                                             .collapse_length_threshold(0.000000000000001)); // no-collapse but flips
  assert(vertices(mesh).size()==vertices(bk).size());
  assert(!same_meshes(mesh,bk));

  mesh=bk;
  PMP::remove_almost_degenerate_faces(mesh,
                                      params::cap_threshold(std::cos(160. / 180 * CGAL_PI))
                                             .needle_threshold(4)
                                             .collapse_length_threshold(0.000000000000001)
                                             .flip_triangle_height_threshold(0.000000000000001)); // no-collapse and no flip
  assert(vertices(mesh).size()==vertices(bk).size());
  assert(same_meshes(mesh,bk));
}

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/pig.off");

  general_test(filename);
  if (argc==2)
    test_with_envelope(filename, 0.01);
  else
    if (argc==3)
      test_with_envelope(filename, atof(argv[2]));

  // only run that test with pig.off
  if (argc==1)
    test_parameters_on_pig(filename);

  return 0;
}
