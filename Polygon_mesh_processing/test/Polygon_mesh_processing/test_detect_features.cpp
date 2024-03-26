#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Real_timer.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;
typedef CGAL::Surface_mesh<K::Point_3>                               Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor                 vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor                   face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

void test_cube()
{
  std::cout << "Test cube..." << std::endl;

  std::ifstream input(CGAL::data_file_path("meshes/cube_quad.off"));
  Mesh mesh;
  if(!input || !(input >> mesh))
  {
    std::cerr << "Failed to read cube_quad.off input." << std::endl;
    assert(false);
    return;
  }

  typedef boost::property_map<Mesh,CGAL::face_patch_id_t<int> >::type PatchID;
  typedef boost::property_map<Mesh, CGAL::vertex_incident_patches_t<int> >::type VIP;
  typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIF_map;

  PatchID pid = get(CGAL::face_patch_id_t<int>(), mesh);
  VIP vip     = get(CGAL::vertex_incident_patches_t<int>(), mesh);
  EIF_map eif = get(CGAL::edge_is_feature, mesh);
  std::size_t number_of_patches = PMP::sharp_edges_segmentation(mesh, 90, eif, pid);

  std::size_t nb_sharp_edges = 0;
  for(boost::graph_traits<Mesh>::edge_descriptor e : edges(mesh))
  {
    if(get(eif, e))
      ++nb_sharp_edges;
  }

  assert(nb_sharp_edges == 12);
  assert(number_of_patches == 6);
  CGAL_USE(number_of_patches);

  number_of_patches = PMP::sharp_edges_segmentation(mesh, 90, eif, pid,
                                                    CGAL::parameters::first_index(1)
                                                                     .vertex_incident_patches_map(vip));

  assert(number_of_patches == 6);


  number_of_patches = PMP::internal::detect_surface_patches(mesh, pid, eif);
  PMP::detect_vertex_incident_patches(mesh, pid, vip, eif);

  nb_sharp_edges = 0;
  for(boost::graph_traits<Mesh>::edge_descriptor e : edges(mesh))
  {
    if(get(eif, e))
      ++nb_sharp_edges;
  }

  assert(nb_sharp_edges == 12);
  assert(number_of_patches == 6);

  Mesh::Property_map<face_descriptor,std::pair<int, int> > patch_id_map
    = mesh.add_property_map<face_descriptor,std::pair<int, int> >("f:pid",std::pair<int,int>()).first;
  Mesh::Property_map<vertex_descriptor,std::set<std::pair<int, int> > > vertex_incident_patch_map
    = mesh.add_property_map<vertex_descriptor,std::set<std::pair<int, int> > >("f:vip",std::set<std::pair<int, int> >()).first;
  PMP::detect_sharp_edges(mesh, 90, eif);
  number_of_patches = PMP::internal::detect_surface_patches(mesh, patch_id_map, eif,
                                                            CGAL::parameters::first_index(1));
  PMP::detect_vertex_incident_patches(mesh, patch_id_map, vertex_incident_patch_map, eif);

  nb_sharp_edges =0;
  for(boost::graph_traits<Mesh>::edge_descriptor e : edges(mesh))
  {
    if(get(eif, e))
      ++nb_sharp_edges;
  }

  assert(nb_sharp_edges == 12);
  assert(number_of_patches == 6);
}

void test_blobby()
{
  std::cout << "Test blobby..." << std::endl;

  std::ifstream input(CGAL::data_file_path("meshes/blobby.off"));
  Mesh mesh;
  if(!input || !(input >> mesh))
  {
    std::cerr << "Failed to read blobby.off input." << std::endl;
    assert(false);
    return;
  }

  typedef CGAL::dynamic_edge_property_t<bool> EIF_tag;
  typedef boost::property_map<Mesh, EIF_tag>::type EIF;

  EIF eif = get(EIF_tag(), mesh);

  CGAL::Real_timer timer;
  timer.start();

  PMP::detect_sharp_edges(mesh, 5, eif);

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << std::endl;

  unsigned int nb_sharp_edges = 0;
  for(auto e : edges(mesh))
    if(get(eif, e))
      ++nb_sharp_edges;

  std::cout << "Found " << nb_sharp_edges << " sharp edges" << std::endl;
  assert(nb_sharp_edges == 2565);
}

int main(int, char**)
{
  test_cube();
  test_blobby();

  return EXIT_SUCCESS;
}
