#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include "VSA.h"
#include <CGAL/property_map.h>

namespace CGAL
{
/*!
\ingroup PkgTSMA
Main function
*/
template <typename TriangleMesh,
  typename SegmentPropertyMap,
  typename FittingPropertyMap,
  typename PointPropertyMap
  = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type,
  typename GeomTraits
  = typename Kernel_traits<typename boost::property_traits<PointPropertyMap>::value>::Kernel>
  void vsa_mesh_approximation(const TriangleMesh &triangle_mesh,
    const std::size_t number_of_segments,
    SegmentPropertyMap segment_ids,
    FittingPropertyMap fit_error_map = FittingPropertyMap(),
    PointPropertyMap ppmap = PointPropertyMap(),
    GeomTraits traits = GeomTraits()) {
  internal::VSA<TriangleMesh, PointPropertyMap, GeomTraits>
    algorithm(triangle_mesh, ppmap, traits);
  algorithm.partition(number_of_segments, segment_ids);
}
}

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
