#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/repair.h>

#include <iostream>
#include <fstream>

//note : when
//CGAL::get_default_random()::get_seed() = 1473902576
//the last test (on trihole.off) does not terminate
//

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel     EPICK;

typedef CGAL::Surface_mesh<EPICK::Point_3>                      Surface_mesh;
typedef CGAL::Polyhedron_3<EPICK>                               Polyhedron;

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
  PMP::degenerate_edges(mesh, std::inserter(dedges, dedges.end()), CGAL::parameters::all_default());

  dedges.clear();
  PMP::degenerate_edges(edge_range, mesh, std::inserter(dedges, dedges.begin()), CGAL::parameters::all_default());
  std::cout << "\t" << dedges.size() << " degenerate edges vs " <<  expected_dedges_n << std::endl;
  assert(dedges.size() == expected_dedges_n);

  // API tests
  std::vector<face_descriptor> dfaces;
  PMP::degenerate_faces(mesh, std::back_inserter(dfaces));
  PMP::degenerate_faces(face_range, mesh, std::back_inserter(dfaces));
  PMP::degenerate_faces(mesh, std::back_inserter(dfaces), CGAL::parameters::all_default());

  dfaces.clear();
  PMP::degenerate_faces(face_range, mesh, std::back_inserter(dfaces), CGAL::parameters::all_default());
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

  return CGAL::Polygon_mesh_processing::remove_degenerate_edges(edge_range, mesh, CGAL::parameters::all_default());
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

  return CGAL::Polygon_mesh_processing::remove_degenerate_faces(face_range, mesh, CGAL::parameters::all_default());
}

template <typename K, typename Mesh>
void test(const char* filename,
          const std::vector<std::size_t>& edges_selection_ids,
          const std::vector<std::size_t>& faces_selection_ids,
          const std::size_t expected_all_degen_edges_n,
          const std::size_t expected_all_degen_faces_n,
          const std::size_t expected_partial_degen_edges_n,
          const std::size_t expected_partial_degen_faces_n,
          const std::size_t expected_post_removal_degen_edges_n,
          const std::size_t expected_post_removal_degen_faces_n)
{
  std::cout << "  test file: " << filename << std::endl;

  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << filename << " is not a valid off file.\n";
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
  /* bool all_removed = */ CGAL::Polygon_mesh_processing::remove_degenerate_edges(mesh, CGAL::parameters::all_default());
  //assert(all_removed);
  assert(CGAL::is_valid_polygon_mesh(mesh));

  mesh = mesh_cpy;
  /* all_removed = */ CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh, CGAL::parameters::all_default());
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

template <typename K, typename Mesh>
void test()
{
  test<K, Mesh>("data_degeneracies/degtri_2dt_1edge_split_twice.off",
                std::initializer_list<std::size_t>({0, 1, 4, 3}), // edge selection
                std::initializer_list<std::size_t>({0}), // face selection
                0, 2, // expected number of degenerate edges/faces in the complete mesh
                0, 1, // expected number of degenerate edges/faces in the selection
                0, 0); // expected number of degenerate edges/faces in the mesh after partial removal

  test<K, Mesh>("data_degeneracies/degtri_four.off",
                std::initializer_list<std::size_t>({1}),
                std::initializer_list<std::size_t>({3}),
                0, 1, 0, 0, 0, 1);

  test<K, Mesh>("data_degeneracies/degtri_four-2.off",
                std::initializer_list<std::size_t>({2}),
                std::initializer_list<std::size_t>({3}),
                0, 1, 0, 0, 0, 1);

  test<K, Mesh>("data_degeneracies/degtri_on_border.off",
                std::initializer_list<std::size_t>({2}),
                std::initializer_list<std::size_t>({0}),
                0, 1, 0, 1, 0, 0);

  test<K, Mesh>("data_degeneracies/degtri_three.off",
                std::initializer_list<std::size_t>({2}),
                std::initializer_list<std::size_t>({1}),
                0, 1, 0, 0, 0, 1);

  test<K, Mesh>("data_degeneracies/degtri_single.off",
                std::initializer_list<std::size_t>({0, 1, 2}),
                std::initializer_list<std::size_t>({0}),
                0, 1, 0, 1, 0, 0);

  test<K, Mesh>("data_degeneracies/degtri_nullface.off",
                std::initializer_list<std::size_t>({3, 6, 7}),
                std::initializer_list<std::size_t>({0, 1, 2}),
                3, 4, 1, 2, 0, 0);

  test<K, Mesh>("data_degeneracies/trihole.off",
                std::initializer_list<std::size_t>({12}),
                std::initializer_list<std::size_t>({4, 5}),
                1, 3, 1, 2, 0, 0);

  test<K, Mesh>("data_degeneracies/degtri_sliding.off",
                std::initializer_list<std::size_t>({2}),
                std::initializer_list<std::size_t>({2, 4}),
                0, 4, 0, 2, 0, 0);

  test<K, Mesh>("data_degeneracies/fused_vertices.off",
                std::initializer_list<std::size_t>({5, 10, 13, 15, 27, 45}),
                std::initializer_list<std::size_t>({1, 3, 5, 10, 19}),
                6, 7, 2, 4, 3, 3);
}

int main()
{
  std::cout << "EPICK SM TESTS" << std::endl;
  test<EPICK, Surface_mesh>();

  std::cout << "EPICK POLYHEDRON TESTS" << std::endl;
  test<EPICK, Polyhedron>();

  return EXIT_SUCCESS;
}
