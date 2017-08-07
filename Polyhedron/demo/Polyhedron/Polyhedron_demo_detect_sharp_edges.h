
#ifndef POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
#define POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H

#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Mesh_3/properties.h>

namespace CGAL
{
  template<typename Polyhedron>
  void reset_sharp_edges(Polyhedron* pMesh)
  {
    typename boost::property_map<Polyhedron,edge_is_feature_t>::type if_pm =
        get(CGAL::edge_is_feature, *pMesh);
    BOOST_FOREACH(typename boost::graph_traits<Polyhedron>::edge_descriptor ed, edges(*pMesh))
    {
      put(if_pm,ed,false);
    }
  }
  template<typename Polyhedron>
  void detect_sharp_edges(Polyhedron* pMesh, const double angle)
  {
    reset_sharp_edges(pMesh);

    // Detect edges in current polyhedron
    typedef typename boost::property_map<Polyhedron,CGAL::face_patch_id_t<int> >::type PIDMap;
    typedef typename boost::property_map<Polyhedron,CGAL::vertex_incident_patches_t<int> >::type VIPMap;
    PIDMap pid_map = get(face_patch_id_t<int>(), *pMesh);
    VIPMap vip_map = get(vertex_incident_patches_t<int>(), *pMesh);


    CGAL::Polygon_mesh_processing::detect_sharp_edges(*pMesh, angle);
  }

}//end namespace CGAL

#endif //POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
