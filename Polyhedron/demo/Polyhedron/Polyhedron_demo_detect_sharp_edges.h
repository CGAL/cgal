
#ifndef POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
#define POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H

#include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>

namespace CGAL
{
  template<typename Polyhedron>
  void reset_sharp_edges(Polyhedron* pMesh)
  {
    for (typename Polyhedron::Edge_iterator
      eit = pMesh->edges_begin(),
      end = pMesh->edges_end(); eit != end; ++eit)
    {
      eit->set_feature_edge(false);
    }
  }

  template<typename Polyhedron>
  void detect_sharp_edges(Polyhedron* pMesh, const double angle)
  {
    reset_sharp_edges(pMesh);

    // Detect edges in current polyhedron
    CGAL::Mesh_3::Detect_features_in_polyhedra<Polyhedron> detect_features;
    detect_features.detect_sharp_edges(*pMesh, angle);
  }

}//end namespace CGAL

#endif //POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
