#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/Detect_features_in_polygon_mesh.h>
#include <CGAL/Mesh_3/properties_Surface_mesh.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/P.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }


  typedef typename boost::property_map<Mesh,CGAL::face_patch_id_t<int> >::type PatchID;
  PatchID pid = get(CGAL::face_patch_id_t<int>(), mesh);
  typedef typename boost::property_map<Mesh,CGAL::vertex_incident_patches_t<int> >::type VIP;
  VIP vip = get(CGAL::vertex_incident_patches_t<int>(), mesh);
  std::size_t number_of_patches = 1;
  PMP::detect_features(mesh, 90, pid, vip, PMP::parameters::maximum_number_of_patches(&number_of_patches));
  typedef boost::property_map<Mesh,CGAL::halfedge_is_feature_t>::type HIF_map;
  HIF_map hif = get(CGAL::halfedge_is_feature, mesh);
  int nb_sharp_edges =0;
  BOOST_FOREACH(boost::graph_traits<Mesh>::halfedge_descriptor h, halfedges(mesh))
  {
    if(get(hif, h))
      ++nb_sharp_edges;
  }


  std::cout<<"This mesh contains "<<nb_sharp_edges/2<<" sharp edges"<<std::endl;
  std::cout<<" and "<<number_of_patches - 1<<" surface patches."<<std::endl;

  return 0;
}
