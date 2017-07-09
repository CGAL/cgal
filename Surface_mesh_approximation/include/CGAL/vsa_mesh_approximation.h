#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/internal/Surface_mesh_approximation/VSA.h>
#include <CGAL/property_map.h>

namespace CGAL
{
/*!
\ingroup PkgTSMA
Main function
*/
template<typename TriangleMesh,
  typename SegmentPropertyMap,
  typename PointPropertyMap,
  typename AnchorIndexContainer,
  typename AnchorPositionContainer,
  typename AnchorVertexContainer,
  typename BoundaryContainer,
  typename GeomTraits>
  void vsa_mesh_approximation(const TriangleMesh &triangle_mesh,
    const std::size_t number_of_segments,
    const std::size_t number_of_iterations,
    SegmentPropertyMap segment_ids,
    PointPropertyMap ppmap,
    AnchorIndexContainer &tris,
    AnchorPositionContainer &pos,
    AnchorVertexContainer &vtx,
    BoundaryContainer &bdrs,
    GeomTraits traits) {
  typedef CGAL::internal::VSA<TriangleMesh, GeomTraits, PointPropertyMap> VSA;

  VSA algorithm(triangle_mesh, ppmap, traits);

  algorithm.partition(number_of_segments, number_of_iterations, segment_ids);

  algorithm.extract_mesh(segment_ids, tris);
  BOOST_FOREACH(const typename VSA::Anchor &a, algorithm.collect_anchors()) {
    vtx.push_back(a.vtx);
    pos.push_back(a.pos);
  }

  bdrs = algorithm.collect_borders(segment_ids);
}
}

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
