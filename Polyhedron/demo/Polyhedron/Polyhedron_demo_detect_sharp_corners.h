#pragma once
#include "detect_sharp_corners.h"
#include <CGAL/Polygon_mesh_processing/detect_features.h>

namespace CGAL
{
  template<typename Polyhedron>
  void reset_sharp_corners(Polyhedron* pMesh)
  {
    typename boost::property_map<Polyhedron,vertex_is_feature_t>::type if_pm =
        get(CGAL::vertex_is_feature, *pMesh);
    for(typename boost::graph_traits<Polyhedron>::vertex_descriptor vd : vertices(*pMesh))
    {
      put(if_pm,vd,false);
    }
  }
  template<typename Polyhedron>
  void detect_sharp_corners(Polyhedron* pMesh, const double angle)
  {
    reset_sharp_edges(pMesh);

    // Detect corners in current polyhedron
    typename boost::property_map<Polyhedron,vertex_is_feature_t>::type vif =
        get(CGAL::vertex_is_feature, *pMesh);
    CGAL::Polygon_mesh_processing::detect_sharp_corners(*pMesh, angle, vif);
  }

}//end namespace CGAL

