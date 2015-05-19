
#ifndef POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
#define POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H

#include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>

namespace CGAL
{
  template<typename Polyhedron>
  void detect_sharp_edges(Polyhedron* pMesh, const double angle)
  {
    CGAL::Mesh_3::Detect_features_in_polyhedra<Polyhedron> detect_features;

    for (typename Polyhedron::Edge_iterator
          eit = pMesh->edges_begin(),
          end = pMesh->edges_end(); eit != end; ++eit)
    {
      eit->set_feature_edge(false);
    }

    // Detect edges in current polyhedron
    detect_features.detect_sharp_edges(*pMesh, angle);
  }
}//end namespace CGAL

#endif //POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
