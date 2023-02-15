#define CGAL_PMP_DEBUG_SMALL_CC_REMOVAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <iostream>
#include <fstream>

//note : when
//CGAL::get_default_random()::get_seed() = 1473902576
//the last test (on trihole.off) does not terminate
//

namespace PMP = CGAL::Polygon_mesh_processing;
namespace CP = CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel             EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel               EPECK;

template <typename K, typename EdgeRange, typename FaceRange, typename Mesh>
void detect_degeneracies(const EdgeRange& edge_range,
                         const FaceRange& face_range,
                         const Mesh& mesh,
                         const std::size_t expected_dedges_n,
                         const std::size_t expected_dfaces_n)
{
  typedef typename boost::graph_traits<Mesh>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor   face_descriptor;

  // API tests
  std::set<edge_descriptor> dedges;
  PMP::degenerate_edges(mesh, std::inserter(dedges, dedges.end()));
  PMP::degenerate_edges(edge_range, mesh, std::inserter(dedges, dedges.begin()));
  PMP::degenerate_edges(mesh, std::inserter(dedges, dedges.end()), CP::default_values());

  dedges.clear();
  PMP::degenerate_edges(edge_range, mesh, std::inserter(dedges, dedges.begin()), CP::default_values());
  std::cout << "\t" << dedges.size() << " degenerate edges vs " <<  expected_dedges_n << std::endl;
  assert(dedges.size() == expected_dedges_n);

  // API tests
  std::vector<face_descriptor> dfaces;
  PMP::degenerate_faces(mesh, std::back_inserter(dfaces));
  PMP::degenerate_faces(face_range, mesh, std::back_inserter(dfaces));
  PMP::degenerate_faces(mesh, std::back_inserter(dfaces), CP::default_values());

  dfaces.clear();
  PMP::degenerate_faces(face_range, mesh, std::back_inserter(dfaces), CP::default_values());
  std::cout << "\t" << dfaces.size() << " degenerate faces vs " << expected_dfaces_n << std::endl;
  assert(dfaces.size() == expected_dfaces_n);
}

template <typename K, typename Mesh>
void detect_degeneracies(const std::vector<std::size_t>& edges_selection_ids,
                         const std::vector<std::size_t>& faces_selection_ids,
                         const Mesh& mesh,
                         const std::size_t expected_dedges_n,
                         const std::size_t expected_dfaces_n)
{
  typedef typename boost::graph_traits<Mesh>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor   face_descriptor;

  // Convert IDs to descriptors
  std::vector<edge_descriptor> all_edges(edges(mesh).begin(), edges(mesh).end());
  std::vector<face_descriptor> all_faces(faces(mesh).begin(), faces(mesh).end());

  std::vector<edge_descriptor> edge_range;
  for(std::size_t edge_id : edges_selection_ids)
    edge_range.push_back(all_edges[edge_id]);

  std::set<face_descriptor> face_range;
  for(std::size_t face_id : faces_selection_ids)
    face_range.insert(all_faces[face_id]);

  return detect_degeneracies<K>(edge_range, face_range, mesh, expected_dedges_n, expected_dfaces_n);
}

template <typename K, typename Mesh>
void detect_degeneracies(const Mesh& mesh,
                         const std::size_t expected_dedges_n,
                         const std::size_t expected_dfaces_n)
{
  return detect_degeneracies<K>(edges(mesh), faces(mesh), mesh, expected_dedges_n, expected_dfaces_n);
}

template <typename Mesh>
bool remove_dedges(const std::vector<std::size_t>& edges_selection_ids,
                   Mesh& mesh)
{
  typedef typename boost::graph_traits<Mesh>::edge_descriptor   edge_descriptor;

  // Convert IDs to descriptors
  std::vector<edge_descriptor> all_edges(edges(mesh).begin(), edges(mesh).end());

  std::vector<edge_descriptor> edge_range;
  for(std::size_t edge_id : edges_selection_ids)
    edge_range.push_back(all_edges[edge_id]);

  return CGAL::Polygon_mesh_processing::remove_degenerate_edges(edge_range, mesh, CP::default_values());
}

template <typename Mesh>
bool remove_dfaces(const std::vector<std::size_t>& faces_selection_ids,
                   Mesh& mesh)
{
  typedef typename boost::graph_traits<Mesh>::face_descriptor   face_descriptor;

  // Convert IDs to descriptors
  std::vector<face_descriptor> all_faces(faces(mesh).begin(), faces(mesh).end());

  std::vector<face_descriptor> face_range;
  for(std::size_t face_id : faces_selection_ids)
    face_range.push_back(all_faces[face_id]);

  return CGAL::Polygon_mesh_processing::remove_degenerate_faces(face_range, mesh, CP::default_values());
}

template <typename K, typename Mesh>
void remove_degeneracies(const std::string filename,
                         const std::vector<std::size_t>& edges_selection_ids,
                         const std::vector<std::size_t>& faces_selection_ids,
                         const std::size_t expected_all_degen_edges_n,
                         const std::size_t expected_all_degen_faces_n,
                         const std::size_t expected_partial_degen_edges_n,
                         const std::size_t expected_partial_degen_faces_n,
                         const std::size_t expected_post_removal_degen_edges_n,
                         const std::size_t expected_post_removal_degen_faces_n)
{
  std::cout << "  remove_degeneracies, file: " << filename << std::endl;

  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << filename << " is not a valid off file." << std::endl;
    exit(1);
  }

  Mesh mesh_cpy = mesh;

  // Count
  std::cout << "    Count..." << std::endl;
  detect_degeneracies<K>(mesh, expected_all_degen_edges_n, expected_all_degen_faces_n);
  detect_degeneracies<K>(edges_selection_ids, faces_selection_ids, mesh, expected_partial_degen_edges_n, expected_partial_degen_faces_n);

  // Complete remove
  std::cout << "    Remove all..." << std::endl;
  mesh = mesh_cpy;
  /* bool all_removed = */ CGAL::Polygon_mesh_processing::remove_degenerate_edges(mesh, CP::default_values());
  //assert(all_removed);
  assert(CGAL::is_valid_polygon_mesh(mesh));

  mesh = mesh_cpy;
  /* all_removed = */ CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh, CP::default_values());
  // assert(all_removed);
  assert(CGAL::is_valid_polygon_mesh(mesh));

  // Check that everything is gone
  std::cout << "    Count (Post All-Remove)..." << std::endl;
  detect_degeneracies<K>(mesh, 0, 0);

  // Partial remove
  std::cout << "    Partial remove..." << std::endl;
  mesh = mesh_cpy;
  remove_dedges(edges_selection_ids, mesh);
  assert(CGAL::is_valid_polygon_mesh(mesh));

  mesh = mesh_cpy;
  remove_dfaces(faces_selection_ids, mesh);
  assert(CGAL::is_valid_polygon_mesh(mesh));

  // Count how much is left after partial removal
  std::cout << "    Count (Post Partial-Remove)..." << std::endl;
  detect_degeneracies<K>(mesh, expected_post_removal_degen_edges_n, expected_post_removal_degen_faces_n);

  std::cout << "  Done" << std::endl;
}

template <typename Mesh>
void initialize_IDs(const Mesh&) { }

template <typename Kernel>
void initialize_IDs(const CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>& mesh)
{
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Mesh;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor           face_descriptor;

  std::size_t i=0;
  for(vertex_descriptor v : vertices(mesh))
    v->id() = i++;

  i=0;
  for(face_descriptor f : faces(mesh))
    f->id() = i++;
}

template <typename K, typename Mesh>
void remove_negligible_connected_components(const std::string filename)
{
  typedef typename boost::graph_traits<Mesh>::face_descriptor           face_descriptor;

  std::cout << "  remove negligible CCs, file: " << filename << std::endl;

  std::ifstream input(filename);
  Mesh mesh, mesh_cpy;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << filename << " is not a valid off file." << std::endl;
    exit(1);
  }

  mesh_cpy = mesh;

  initialize_IDs(mesh);
  initialize_IDs(mesh_cpy);

  std::size_t ini_nv = num_vertices(mesh);
  std::size_t ini_nf = num_faces(mesh);

  std::cout << "before: " << ini_nv << " nv & " << ini_nf << " nf" << std::endl;

  // negative thresholds --> doesn't remove anything
  std::cout << "---------\nnull or negative threshold, nothing happens..." << std::endl;
  std::size_t res = PMP::remove_connected_components_of_negligible_size(mesh, CP::area_threshold(-1e15)
                                                                                 .volume_threshold(0));
  assert(PMP::internal::number_of_connected_components(mesh) == 4);
  assert(num_vertices(mesh) == ini_nv);
  assert(num_faces(mesh) == ini_nf);

  // threshold too small, doesn't remove anything
  std::cout << "---------\ntiny threshold, nothing happens..." << std::endl;
  PMP::remove_connected_components_of_negligible_size(mesh, CP::area_threshold(1e-15)
                                                               .volume_threshold(1e-15));
  assert(PMP::internal::number_of_connected_components(mesh) == 4);
  assert(num_vertices(mesh) == ini_nv);
  assert(num_faces(mesh) == ini_nf);

  // that removes the CCs with small volumes
  std::cout << "---------\nremove small volumes..." << std::endl;
  res = PMP::remove_connected_components_of_negligible_size(mesh, CP::area_threshold(1e-15)
                                                                     .volume_threshold(0.1));
  std::cout << "res: " << res << std::endl;
  initialize_IDs(mesh);
  assert(PMP::internal::number_of_connected_components(mesh) == 2);

  // that removes the open CC with a small area
  std::cout << "---------\nremove small areas..." << std::endl;
  PMP::remove_connected_components_of_negligible_size(mesh, CP::area_threshold(20));
  initialize_IDs(mesh);
  assert(PMP::internal::number_of_connected_components(mesh) == 1);

  // Remove everything with a too-large value
  std::cout << "---------\nremove everything..." << std::endl;
  PMP::remove_connected_components_of_negligible_size(mesh, CP::area_threshold(1e15));
  assert(is_empty(mesh));

  // Could also have used default parameters, which does the job by itself
  std::cout << "---------\ndefault values..." << std::endl;

  std::vector<face_descriptor> faces_to_be_removed;
  std::size_t nb_to_be_rm = PMP::remove_connected_components_of_negligible_size(
                              mesh_cpy, CP::dry_run(true)
                                           .output_iterator(std::back_inserter(faces_to_be_removed)));
  assert(nb_to_be_rm == 3);
  assert(PMP::internal::number_of_connected_components(mesh_cpy) == 4); // a dry run does not remove anything
  assert(faces_to_be_removed.size() == 216); // sum of #faces of the small CCs

  assert(nb_to_be_rm == PMP::remove_connected_components_of_negligible_size(mesh_cpy));
  assert(PMP::internal::number_of_connected_components(mesh_cpy) == 1);
}

template <typename K, typename Mesh>
void test()
{
  remove_degeneracies<K, Mesh>("data_degeneracies/degtri_2dt_1edge_split_twice.off",
                               std::initializer_list<std::size_t>({0, 1, 4, 3}), // edge selection
                               std::initializer_list<std::size_t>({0}), // face selection
                               0, 2, // expected number of degenerate edges/faces in the complete mesh
                               0, 1, // expected number of degenerate edges/faces in the selection
                               0, 0); // expected number of degenerate edges/faces in the mesh after partial removal

  remove_degeneracies<K, Mesh>("data_degeneracies/degtri_four.off",
                               std::initializer_list<std::size_t>({1}),
                               std::initializer_list<std::size_t>({3}),
                               0, 1, 0, 0, 0, 1);

  remove_degeneracies<K, Mesh>("data_degeneracies/degtri_four-2.off",
                               std::initializer_list<std::size_t>({2}),
                               std::initializer_list<std::size_t>({3}),
                               0, 1, 0, 0, 0, 1);

  remove_degeneracies<K, Mesh>("data_degeneracies/degtri_on_border.off",
                               std::initializer_list<std::size_t>({2}),
                               std::initializer_list<std::size_t>({0}),
                               0, 1, 0, 1, 0, 0);

  remove_degeneracies<K, Mesh>("data_degeneracies/degtri_three.off",
                               std::initializer_list<std::size_t>({2}),
                               std::initializer_list<std::size_t>({1}),
                               0, 1, 0, 0, 0, 1);

  remove_degeneracies<K, Mesh>("data_degeneracies/degtri_single.off",
                               std::initializer_list<std::size_t>({0, 1, 2}),
                               std::initializer_list<std::size_t>({0}),
                               0, 1, 0, 1, 0, 0);

  remove_degeneracies<K, Mesh>("data_degeneracies/degtri_nullface.off",
                               std::initializer_list<std::size_t>({3, 6, 7}),
                               std::initializer_list<std::size_t>({0, 1, 2}),
                               3, 4, 1, 2, 0, 0);

  remove_degeneracies<K, Mesh>("data_degeneracies/trihole.off",
                               std::initializer_list<std::size_t>({12}),
                               std::initializer_list<std::size_t>({4, 5}),
                               1, 3, 1, 2, 0, 0);

  remove_degeneracies<K, Mesh>(CGAL::data_file_path("meshes/degtri_sliding.off"),
                               std::initializer_list<std::size_t>({2}),
                               std::initializer_list<std::size_t>({2, 4}),
                               0, 4, 0, 2, 0, 0);

  remove_degeneracies<K, Mesh>("data_degeneracies/fused_vertices.off",
                               std::initializer_list<std::size_t>({5, 10, 13, 15, 27, 45}),
                               std::initializer_list<std::size_t>({1, 3, 5, 10, 19}),
                               6, 7, 2, 4, 3, 3);

  remove_negligible_connected_components<K, Mesh>("data_degeneracies/small_ccs.off");
}

template <typename Kernel>
void test()
{
  typedef typename Kernel::Point_3                                      Point_3;
  typedef CGAL::Surface_mesh<Point_3>                                   Surface_mesh;
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Polyhedron_with_ID;

  std::cout << "SM TESTS" << std::endl;
  test<Kernel, Surface_mesh>();

  std::cout << "POLYHEDRON TESTS" << std::endl;
  test<Kernel, Polyhedron_with_ID>();
}

int main(int /*argc*/, char** /*argv*/)
{
  std::cout << "EPICK TESTS" << std::endl;
  test<EPICK>();

  std::cout << "EPECK TESTS" << std::endl;
  test<EPECK>();

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
