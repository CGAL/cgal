#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/P.off";

  Mesh mesh;
  if(!PMP::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIFMap;
  typedef boost::property_map<Mesh, CGAL::face_patch_id_t<int> >::type PIMap;
  typedef boost::property_map<Mesh, CGAL::vertex_incident_patches_t<int> >::type VIMap;
  typedef boost::property_map<Mesh, CGAL::vertex_is_feature_t>::type VIFMap;

  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  PIMap pid = get(CGAL::face_patch_id_t<int>(), mesh);
  VIMap vip = get(CGAL::vertex_incident_patches_t<int>(), mesh);
  VIFMap vif = get(CGAL::vertex_is_feature, mesh);

  std::size_t number_of_patches
    = PMP::sharp_edges_segmentation(mesh, 90, eif, pid,
                                    PMP::parameters::vertex_incident_patches_map(vip));
  
  PMP::detect_sharp_corners(60, vif, mesh);
  std::size_t sharp_corners_counter = 0;
  for(boost::graph_traits<Mesh>::vertex_descriptor v : vertices(mesh))
  {
    if(get(vif, v))
      ++sharp_corners_counter;
  }

  PMP::detect_sharp_edges(mesh, 60, eif);
  std::size_t nb_sharp_edges = 0;
  for(boost::graph_traits<Mesh>::edge_descriptor e : edges(mesh))
  {
    if(get(eif, e))
      ++nb_sharp_edges;
  }

  std::cout<<"This mesh contains "<<nb_sharp_edges<<" sharp edges"<<std::endl;
  std::cout<<sharp_corners_counter<<" sharp corners"<<std::endl;
  std::cout<<" and "<<number_of_patches<<" surface patches."<<std::endl;

  return 0;
}
