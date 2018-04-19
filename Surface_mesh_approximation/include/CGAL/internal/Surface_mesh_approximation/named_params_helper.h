#ifndef CGAL_SMA_NAMED_PARAMETERS_HELPERS_H
#define CGAL_SMA_NAMED_PARAMETERS_HELPERS_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>

#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/mpl/if.hpp>

namespace CGAL {
namespace Surface_mesh_approximation {
namespace internal_np {

// shortcut for accessing the value type of the property map
template <class Graph, class Property>
class property_map_value {
  typedef typename boost::property_map<Graph, Property>::const_type PMap;
public:
  typedef typename boost::property_traits<PMap>::value_type type;
};

template<typename PolygonMesh, typename NamedParameters>
class GetVertexPointMap
{
  typedef typename property_map_selector<PolygonMesh, boost::vertex_point_t>::const_type
    DefaultVPMap_const;
  typedef typename property_map_selector<PolygonMesh, boost::vertex_point_t>::type
    DefaultVPMap;
public:
  typedef typename boost::lookup_named_param_def<
    internal_np::vertex_point_t,
    NamedParameters,
    DefaultVPMap
  > ::type  type;
  typedef typename boost::lookup_named_param_def<
    internal_np::vertex_point_t,
    NamedParameters,
    DefaultVPMap_const
  > ::type  const_type;
};

template<typename PolygonMesh, typename NamedParameters>
class GetK
{
  typedef typename boost::property_traits<
    typename GetVertexPointMap<PolygonMesh, NamedParameters>::type
  >::value_type Point;
public:
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
};

template<typename PolygonMesh, typename NamedParameters>
class GetGeomTraits
{
  typedef typename boost::graph_has_property<PolygonMesh, boost::vertex_point_t>::type
    Has_internal_pmap;
  struct Fake_GT {};//to be used if there is no internal vertex_point_map in PolygonMesh

  typedef typename boost::mpl::if_c< Has_internal_pmap::value
                                   , typename GetK<PolygonMesh, NamedParameters>::Kernel
                                   , Fake_GT
  >::type DefaultKernel;

public:
  typedef typename boost::lookup_named_param_def <
    internal_np::geom_traits_t,
    NamedParameters,
    DefaultKernel
  > ::type  type;
};

// output helper functions
template <typename Approximation, typename FacetProxyMap>
void facet_proxy_map_helper(const Approximation &approx, FacetProxyMap fproxymap) {
  approx.proxy_map(fproxymap);
}

template <typename Approximation>
void facet_proxy_map_helper(const Approximation &, internal_np::dummy_output_t) {}

// proxies

template <typename Approximation, typename OutputIterator>
void proxies_helper(const Approximation &approx, OutputIterator out) 
{
  approx.proxies(out);
}

template <typename Approximation>
void proxies_helper(const Approximation &, internal_np::dummy_output_t) 
{}

// anchors 

template <typename Approximation, typename OutputIterator>
void anchors_helper(const Approximation &approx, OutputIterator out)
{
  approx.anchor_points(out);
}

template <typename Approximation>
void anchors_helper(const Approximation &, internal_np::dummy_output_t)
{}

// indexed triangles

template <typename Approximation, typename OutputIterator>
void triangles_helper(const Approximation &approx, OutputIterator out)
{
  approx.indexed_triangles(out);
}

template <typename Approximation>
void triangles_helper(const Approximation &, internal_np::dummy_output_t)
{}

} // namespace internal_np
} // namespace Surface_mesh_approximation
} // namespace CGAL

#endif //CGAL_SMA_NAMED_PARAMETERS_HELPERS_H
