
#ifndef POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
#define POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H

#include <CGAL/Polygon_mesh_processing/detect_features.h>

namespace CGAL
{
  template<typename Polyhedron>
  void reset_sharp_edges(Polyhedron* pMesh)
  {
    typename boost::property_map<Polyhedron,edge_is_feature_t>::type if_pm =
        get(CGAL::edge_is_feature, *pMesh);
    for(typename boost::graph_traits<Polyhedron>::edge_descriptor ed : edges(*pMesh))
    {
      put(if_pm,ed,false);
    }
  }
  template<typename Polyhedron>
  void detect_sharp_edges(Polyhedron* pMesh, const double angle)
  {
    reset_sharp_edges(pMesh);

    // Detect edges in current polyhedron
    typename boost::property_map<Polyhedron,edge_is_feature_t>::type eif =
        get(CGAL::edge_is_feature, *pMesh);
    CGAL::Polygon_mesh_processing::detect_sharp_edges(*pMesh, angle, eif);
  }

}//end namespace CGAL

#endif //POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
