#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polyhedral_envelope.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor     face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

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

  PMP::experimental::remove_almost_degenerate_faces(mesh,
                                                    std::cos(160. / 180 * CGAL_PI),
                                                    4,
                                                    0.14);


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
  PMP::experimental::remove_almost_degenerate_faces(mesh,
                                                    std::cos(160. / 180 * CGAL_PI),
                                                    4,
                                                    0.14,
                                                    CGAL::parameters::filter(no_modif));
  assert(nbv == vertices(mesh).size());

  // now the real test with a fixed envelope
  std::cout << "Using fixed envelope\n";
  std::cout << "  Input mesh has " << edges(mesh).size() << " edges\n";
  bk=mesh;
  Envelope envelope(mesh, eps);
  PMP::experimental::remove_almost_degenerate_faces(mesh,
                                                    std::cos(160. / 180 * CGAL_PI),
                                                    4,
                                                    0.14,
                                                    CGAL::parameters::filter(std::ref(envelope)));

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
  PMP::experimental::remove_almost_degenerate_faces(mesh,
                                                    std::cos(160. / 180 * CGAL_PI),
                                                    4,
                                                    0.14,
                                                    CGAL::parameters::filter(filter));

  CGAL::IO::write_polygon_mesh("cleaned_mesh_with_iterative_envelope.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "  Output mesh has " << edges(mesh).size() << " edges\n";
  if (PMP::does_self_intersect(mesh))
    std::cout << "  Output mesh has self-intersections\n";
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

  return 0;
}
