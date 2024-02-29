#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <cassert>
#include <fstream>
#include <map>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef CGAL::Surface_mesh<K::Point_3>                            Surface_mesh;
typedef CGAL::Polyhedron_3<K>                                     Polyhedron;

typedef std::vector<std::vector<std::size_t> >                    Vertices_to_merge_container;

template <typename PolygonMesh>
void read_mesh(const std::string fname,
               PolygonMesh& mesh)
{
  std::ifstream input(fname);
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file." << std::endl;
    std::exit(1);
  }
}

// tests merge_and_duplication
template <typename PolygonMesh>
void merge_vertices(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_keep,
                    typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_rm,
                    PolygonMesh& mesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  assert(v_keep != v_rm);

  std::size_t ini_nv = static_cast<std::size_t>(vertices(mesh).size());

  halfedge_descriptor h = halfedge(v_rm, mesh);
  halfedge_descriptor start = h;
  do
  {
    set_target(h, v_keep, mesh);
    h = opposite(next(h, mesh), mesh);
  }
  while(h != start);

  remove_vertex(v_rm, mesh);

  assert(vertices(mesh).size() == ini_nv - 1);
}

template <typename PolygonMesh>
void merge_vertices(const Vertices_to_merge_container& all_vertices_to_merge,
                    std::map<typename boost::graph_traits<PolygonMesh>::vertex_descriptor, std::size_t>& merged_onto,
                    PolygonMesh& mesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;

  // int to vd
  std::vector<vertex_descriptor> vds(vertices(mesh).begin(), vertices(mesh).end());

  for(std::size_t i=0, avtmn=all_vertices_to_merge.size(); i<avtmn; ++i)
  {
    const std::vector<std::size_t>& vertices_to_merge = all_vertices_to_merge[i];
    if(vertices_to_merge.size() <= 1)
      continue;

    vertex_descriptor vd_to_merge_onto = vds[vertices_to_merge[0]];

    for(std::size_t j=1, vtmn=vertices_to_merge.size(); j<vtmn; ++j)
      merge_vertices(vd_to_merge_onto, vds[vertices_to_merge[j]], mesh);

    std::pair<typename std::map<vertex_descriptor, std::size_t>::iterator, bool > is_insert_successful =
      merged_onto.insert(std::make_pair(vd_to_merge_onto, vertices_to_merge.size() - 1));
    if(!is_insert_successful.second)
      is_insert_successful.first->second += vertices_to_merge.size() - 1;

    assert(CGAL::Polygon_mesh_processing::is_non_manifold_vertex(vd_to_merge_onto, mesh));
  }
}

template <typename PolygonMesh>
std::size_t test_nm_vertices_duplication(const Vertices_to_merge_container& all_merges,
                                         std::map<typename boost::graph_traits<PolygonMesh>::vertex_descriptor, std::size_t>& merged_onto,
                                         std::vector<std::vector<
                                           typename boost::graph_traits<PolygonMesh>::vertex_descriptor> >& duplicated_vertices,
                                         PolygonMesh& mesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor              halfedge_descriptor;

  merge_vertices(all_merges, merged_onto, mesh);

  std::size_t new_vertices_nb =
    CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(mesh,
      CGAL::parameters::output_iterator(std::back_inserter(duplicated_vertices)));

  std::vector<halfedge_descriptor> non_manifold_cones;
  CGAL::Polygon_mesh_processing::non_manifold_vertices(mesh, std::back_inserter(non_manifold_cones));
  assert(non_manifold_cones.empty());

  assert(CGAL::is_valid_polygon_mesh(mesh));

  return new_vertices_nb;
}

template <typename PolygonMesh>
std::size_t test_nm_vertices_duplication(const Vertices_to_merge_container& all_merges,
                                         std::vector<std::vector<
                                           typename boost::graph_traits<PolygonMesh>::vertex_descriptor> >& duplicated_vertices,
                                         PolygonMesh& mesh)
{
  std::map<typename boost::graph_traits<PolygonMesh>::vertex_descriptor, std::size_t> useless_map;

  return test_nm_vertices_duplication(all_merges, useless_map, duplicated_vertices, mesh);
}

template <typename PolygonMesh>
std::size_t test_nm_vertices_duplication(std::vector<std::vector<
                                           typename boost::graph_traits<PolygonMesh>::vertex_descriptor> >& duplicated_vertices,
                                         PolygonMesh& mesh)
{
  Vertices_to_merge_container all_merges;
  std::map<typename boost::graph_traits<PolygonMesh>::vertex_descriptor, std::size_t> useless_map;

  return test_nm_vertices_duplication(all_merges, useless_map, duplicated_vertices, mesh);
}

template <typename PolygonMesh>
void test_unpinched_mesh(const Vertices_to_merge_container& all_merges,
                         PolygonMesh& mesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;

  const std::size_t ini_vertices_size = num_vertices(mesh);

  // this is a nice, smooth, surface initially
  for(vertex_descriptor vd : vertices(mesh)) {
    assert(!CGAL::Polygon_mesh_processing::is_non_manifold_vertex(vd, mesh));
  }

  std::vector<std::vector<vertex_descriptor> > duplicated_vertices;
  std::map<vertex_descriptor, std::size_t> number_of_vertices_merged_onto;
  std::size_t nb = test_nm_vertices_duplication(all_merges,
                                                number_of_vertices_merged_onto,
                                                duplicated_vertices,
                                                mesh);

  const std::size_t final_vertices_size = vertices(mesh).size();
  std::cout << "    ini: " << ini_vertices_size << " final: " << final_vertices_size << std::endl;
  assert(final_vertices_size == ini_vertices_size);

  std::size_t expected_nb = 0;
  for(std::size_t i=0, n=all_merges.size(); i<n; ++i)
    expected_nb += all_merges[i].size() - 1;
  assert(nb == expected_nb);

  assert(duplicated_vertices.size() == all_merges.size());
  for(std::size_t i=0, n=duplicated_vertices.size(); i<n; ++i)
    assert(duplicated_vertices[i].size() - 1 == number_of_vertices_merged_onto[duplicated_vertices[i][0]]);
}

template <typename PolygonMesh>
void test_blobby()
{
  std::cout << "  test: data/blobby.off" << std::endl;

  PolygonMesh mesh;
  read_mesh(CGAL::data_file_path("meshes/blobby.off"), mesh);

  // non-manifold vertices
  Vertices_to_merge_container all_merges;
  std::vector<std::size_t> single_merge;
  single_merge.push_back(1); single_merge.push_back(7); single_merge.push_back(14); single_merge.push_back(21);
  all_merges.push_back(single_merge);

  single_merge.clear();
  single_merge.push_back(2); single_merge.push_back(8);
  all_merges.push_back(single_merge);

  test_unpinched_mesh(all_merges, mesh);
}

template <typename PolygonMesh>
void test_nm_cubes()
{
  std::cout << "  test: data_repair/nm_closed_cubes.off" << std::endl;

  PolygonMesh mesh;
  read_mesh("data_repair/nm_closed_cubes.off", mesh);

  // non-manifold vertices
  Vertices_to_merge_container all_merges;
  std::vector<std::size_t> single_merge;
  single_merge.push_back(5); single_merge.push_back(14);
  all_merges.push_back(single_merge);

  single_merge.clear();
  single_merge.push_back(6); single_merge.push_back(15);
  all_merges.push_back(single_merge);

  test_unpinched_mesh(all_merges, mesh);
}

template <typename PolygonMesh>
void test_pinched_triangles(const std::string filename,
                            const std::size_t expected_nb)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;

  std::cout << "  test: " << filename << " expected: " << expected_nb << std::endl;

  PolygonMesh mesh;
  read_mesh(filename, mesh);

  // in the triangles, the second (id==1) vertex is non-manifold because it is pinched
  int id = 0;
  for(vertex_descriptor vd : vertices(mesh)) {
    if(id++ == 1) {
      assert(CGAL::Polygon_mesh_processing::is_non_manifold_vertex(vd, mesh));
    } else {
      assert(!CGAL::Polygon_mesh_processing::is_non_manifold_vertex(vd, mesh));
    }
  }

  std::vector<std::vector<vertex_descriptor> > duplicated_vertices;
  std::size_t new_vertices_nb =
    CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(mesh,
      CGAL::parameters::output_iterator(std::back_inserter(duplicated_vertices)));
  std::cout << "    new_vertices_nb: " << new_vertices_nb << " vs expected: " << expected_nb << std::endl;
  assert(new_vertices_nb == expected_nb);
  assert(duplicated_vertices.size() == 1);
  assert(duplicated_vertices[0].size() == 1 + expected_nb);
}

template <typename PolygonMesh>
void test_many_umbrellas()
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;

  std::cout << "  test: data_repair/many_umbrellas.off" << std::endl;

  PolygonMesh mesh;
  read_mesh("data_repair/many_umbrellas.off", mesh);

  // non-manifold vertices
  Vertices_to_merge_container all_merges;
  std::vector<std::size_t> single_merge;
  single_merge.push_back(1); single_merge.push_back(9); single_merge.push_back(15);
  all_merges.push_back(single_merge);

  std::vector<std::vector<vertex_descriptor> > duplicated_vertices;
  std::size_t nb = test_nm_vertices_duplication(all_merges, duplicated_vertices, mesh);
  assert(nb == 5);

  const std::size_t final_vertices_size = vertices(mesh).size();

  assert(duplicated_vertices.size() == 1);
  assert(duplicated_vertices[0].size() == 6);
  assert(final_vertices_size == 19); // 5 new ones, but we merged 2 before, so +3 from '16' at the start
}

template <typename PolygonMesh>
void test_torso()
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;

  std::cout << "  test: data_repair/torso.off" << std::endl;

  PolygonMesh mesh;
  read_mesh("data_repair/torso.off", mesh);

  CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);

  std::vector<std::vector<vertex_descriptor> > duplicated_vertices;
  std::size_t nb = test_nm_vertices_duplication(duplicated_vertices, mesh);

  std::cout << "    new vertices: " << nb << std::endl;
}

int main(int /*argc*/, char** /*argv*/)
{
  std::cout << "Test Vertex Manifoldness Functions (SM)" << std::endl;

  test_blobby<Surface_mesh>(); // data/blobby.off
  test_nm_cubes<Surface_mesh>(); // data_repair/nm_closed_cubes.off
  test_pinched_triangles<Surface_mesh>("data_repair/two_triangles_sharing_a_vertex.off", 1);
  test_pinched_triangles<Surface_mesh>("data_repair/three_triangles_sharing_a_vertex.off", 2);
  test_many_umbrellas<Surface_mesh>(); // data_repair/many_umbrellas.off
  test_torso<Surface_mesh>(); // data_repair/torso.off (only for SM because Polyhedron cannot even read it)

  std::cout << "Test Vertex Manifoldness Functions (Polyhedron)" << std::endl;

  test_blobby<Polyhedron>(); // data/blobby.off
  test_nm_cubes<Polyhedron>(); // data_repair/nm_closed_cubes.off
  test_pinched_triangles<Polyhedron>("data_repair/two_triangles_sharing_a_vertex.off", 1);
  test_pinched_triangles<Polyhedron>("data_repair/three_triangles_sharing_a_vertex.off", 2);
  test_many_umbrellas<Polyhedron>(); // data_repair/many_umbrellas.off

  return EXIT_SUCCESS;
}
