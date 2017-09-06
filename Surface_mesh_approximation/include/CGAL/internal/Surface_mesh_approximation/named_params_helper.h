#ifndef CGAL_NAMED_PARAMETERS_HELPERS_H
#define CGAL_NAMED_PARAMETERS_HELPERS_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/internal/Surface_mesh_approximation/named_function_params.h>

#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/mpl/if.hpp>

namespace CGAL {

// shortcut for accessing the value type of the property map
template <class Graph, class Property>
class property_map_value {
  typedef typename boost::property_map<Graph, Property>::const_type PMap;
public:
  typedef typename boost::property_traits<PMap>::value_type type;
};

template<typename PolygonMesh, typename PropertyTag>
class property_map_selector
{
public:
  typedef typename boost::graph_has_property<PolygonMesh, PropertyTag>::type Has_internal_pmap;
  typedef typename boost::mpl::if_c< Has_internal_pmap::value
                          , typename boost::property_map<PolygonMesh, PropertyTag>::type
                          , typename boost::cgal_no_property::type
  >::type type;
  typedef typename boost::mpl::if_c< Has_internal_pmap::value
                          , typename boost::property_map<PolygonMesh, PropertyTag>::const_type
                          , typename boost::cgal_no_property::const_type
  >::type const_type;

  type get_pmap(const PropertyTag& p, PolygonMesh& pmesh)
  {
    return get_impl(p, pmesh, Has_internal_pmap());
  }

  const_type get_const_pmap(const PropertyTag& p, const PolygonMesh& pmesh)
  {
    return get_const_pmap_impl(p, pmesh, Has_internal_pmap());
  }

private:
  type get_impl(const PropertyTag&, PolygonMesh&, CGAL::Tag_false)
  {
    return type(); //boost::cgal_no_property::type
  }
  type get_impl(const PropertyTag& p, PolygonMesh& pmesh, CGAL::Tag_true)
  {
    return get(p, pmesh);
  }

  const_type get_const_pmap_impl(const PropertyTag&
                               , const PolygonMesh&, CGAL::Tag_false)
  {
    return const_type(); //boost::cgal_no_property::type
  }
  const_type get_const_pmap_impl(const PropertyTag& p
                               , const PolygonMesh& pmesh, CGAL::Tag_true)
  {
    return get(p, pmesh);
  }
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

template<typename PolygonMesh, typename PropertyTag>
typename property_map_selector<PolygonMesh, PropertyTag>::type
get_property_map(const PropertyTag& p, PolygonMesh& pmesh)
{
  property_map_selector<PolygonMesh, PropertyTag> pms;
  return pms.get_pmap(p, pmesh);
}

template<typename PolygonMesh, typename PropertyTag>
typename property_map_selector<PolygonMesh, PropertyTag>::const_type
get_const_property_map(const PropertyTag& p, const PolygonMesh& pmesh)
{
  property_map_selector<PolygonMesh, PropertyTag> pms;
  return pms.get_const_pmap(p, pmesh);
}

// output helper functions
template <typename Approximation, typename FacetProxyMap>
void get_proxy_map(const Approximation &approx, FacetProxyMap fproxymap) {
  approx.get_proxy_map(fproxymap);
}

template <typename Approximation>
void get_proxy_map(const Approximation &, internal_np::vsa_no_output_t) {}

template <typename Approximation, typename OutputIterator>
void get_anchor_vertices(const Approximation &approx, OutputIterator out_itr) {
  approx.get_anchor_vertices(out_itr);
}

template <typename Approximation>
void get_anchor_vertices(const Approximation &, internal_np::vsa_no_output_t) {}

template <typename Approximation, typename OutputIterator>
void get_anchor_points(const Approximation &approx, OutputIterator out_itr) {
  approx.get_anchor_points(out_itr);
}

template <typename Approximation>
void get_anchor_points(const Approximation &, internal_np::vsa_no_output_t) {}

template <typename Approximation, typename OutputIterator>
void get_indexed_triangles(const Approximation &approx, OutputIterator out_itr) {
  approx.get_indexed_triangles(out_itr);
}

template <typename Approximation>
void get_indexed_triangles(const Approximation &, internal_np::vsa_no_output_t) {}

template <typename Approximation, typename OutputIterator>
void get_proxies(const Approximation &approx, OutputIterator out_itr) {
  approx.get_proxies(out_itr);
}

template <typename Approximation>
void get_proxies(const Approximation &, internal_np::vsa_no_output_t) {}

} //end of namespace CGAL

#endif //CGAL_NAMED_PARAMETERS_HELPERS_H
