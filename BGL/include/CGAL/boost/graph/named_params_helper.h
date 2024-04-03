// Copyright (c) 2007-2015  GeometryFactory (France).  All rights reserved.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri, Fernando Cacciola, Jane Tournois

#ifndef CGAL_BOOST_GRAPH_NAMED_PARAMETERS_HELPERS_H
#define CGAL_BOOST_GRAPH_NAMED_PARAMETERS_HELPERS_H

#include <CGAL/boost/graph/internal/initialized_index_maps_helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/iterator.h>
#include <CGAL/Default.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>

#include <boost/mpl/has_xxx.hpp>

#include <fstream>
#include <iterator>
#include <type_traits>

namespace CGAL {

// forward declarations to avoid dependency to Solver_interface
template <typename FT, unsigned int dim>
class Default_diagonalize_traits;
class Eigen_svd;
class Lapack_svd;
struct Alpha_expansion_boost_adjacency_list_tag;
//

//helper classes
template<typename PolygonMesh, typename PropertyTag>
class property_map_selector
{
public:
  typedef typename graph_has_property<PolygonMesh, PropertyTag>::type Has_internal_pmap;
  typedef std::conditional_t<Has_internal_pmap::value,
                             typename boost::property_map<PolygonMesh, PropertyTag>::type,
                             typename boost::cgal_no_property::type> type;
  typedef std::conditional_t<Has_internal_pmap::value,
                             typename boost::property_map<PolygonMesh, PropertyTag>::const_type,
                             typename boost::cgal_no_property::const_type
                             > const_type;

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

// Shortcut for accessing the value type of the property map
template <class Graph, class Property>
class property_map_value
{
  typedef typename boost::property_map<Graph, Property>::const_type PMap;

public:
  typedef typename boost::property_traits<PMap>::value_type type;
};


template <typename PolygonMesh,
          typename VPM_from_NP>
struct GetVertexPointMap_impl
{
  typedef VPM_from_NP type;
  typedef VPM_from_NP const_type;

  template<class NamedParameters>
  static const_type
  get_const_map(const NamedParameters& np, const PolygonMesh&)
  {
    return parameters::get_parameter(np, internal_np::vertex_point);
  }

  template<class NamedParameters>
  static type
  get_map(const NamedParameters& np, PolygonMesh&)
  {
    return parameters::get_parameter(np, internal_np::vertex_point);
  }
};

template <typename PolygonMesh>
struct GetVertexPointMap_impl<PolygonMesh, internal_np::Param_not_found>
{
  typedef typename property_map_selector<PolygonMesh, boost::vertex_point_t>::const_type const_type;
  typedef typename property_map_selector<PolygonMesh, boost::vertex_point_t>::type type;

  template<class NamedParameters>
  static const_type
  get_const_map(const NamedParameters& /* np */, const PolygonMesh& pm)
  {
    return get_const_property_map(boost::vertex_point, pm);
  }

  template<class NamedParameters>
  static type
  get_map(const NamedParameters& /* np */, PolygonMesh& pm)
  {
    return get_property_map(boost::vertex_point, pm);
  }
};

template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
class GetVertexPointMap
{
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_point_t,
                                                       NamedParameters,
                                                       internal_np::Param_not_found>::type VPM_from_NP;

  typedef GetVertexPointMap_impl<PolygonMesh, VPM_from_NP> Impl;

public:
  typedef typename Impl::type type;
  typedef typename Impl::const_type const_type;

  static const_type
  get_const_map(const NamedParameters& np, const PolygonMesh& pm)
  {
    return Impl::get_const_map(np, pm);
  }

  static type
  get_map(const NamedParameters& np, PolygonMesh& pm)
  {
    return Impl::get_map(np, pm);
  }
};

template<typename PolygonMesh, typename NamedParameters>
class GetK
{
  typedef typename boost::property_traits<
    typename GetVertexPointMap<PolygonMesh, NamedParameters>::type>::value_type Point;

public:
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
};


template<typename PolygonMesh, class GT, class NamedParametersVPM>
struct GetGeomTraits_impl
{
  typedef GT type;
};

template<typename PolygonMesh, class NamedParametersVPM>
struct GetGeomTraits_impl<PolygonMesh, internal_np::Param_not_found, NamedParametersVPM>
{
  typedef typename CGAL::graph_has_property<PolygonMesh, boost::vertex_point_t>::type Has_internal_pmap;

  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_point_t,
                                                       NamedParametersVPM,
                                                       internal_np::Param_not_found>::type NP_vpm;

  struct Fake_GT {}; // to be used if there is no internal vertex_point_map in PolygonMesh

  typedef std::conditional_t<Has_internal_pmap::value ||
                             !std::is_same<internal_np::Param_not_found, NP_vpm>::value,
                             typename GetK<PolygonMesh, NamedParametersVPM>::Kernel,
                             Fake_GT> type;
};

template <typename PolygonMesh,
          typename NamedParametersGT = parameters::Default_named_parameters,
          typename NamedParametersVPM = NamedParametersGT>
struct GetGeomTraits
{
  typedef typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                       NamedParametersGT,
                                                       internal_np::Param_not_found>::type GT_from_NP;
  typedef typename GetGeomTraits_impl<PolygonMesh,
                                          GT_from_NP,
                                          NamedParametersVPM>::type type;
};

// Define the following structs:
//
// GetInitializedVertexIndexMap
// GetInitializedHalfedgeIndexMap
// GetInitializedEdgeIndexMap
// GetInitializedFaceIndexMap

#define CGAL_DEF_GET_INDEX_TYPE(CTYPE, DTYPE, STYPE)                                               \
template <typename Graph,                                                                          \
        typename NamedParameters = parameters::Default_named_parameters>                           \
struct GetInitialized##CTYPE##IndexMap                                                             \
  : public BGL::internal::GetInitializedIndexMap<internal_np::DTYPE##_index_t,                     \
                                                 boost::DTYPE##_index_t,                           \
                                                 CGAL::dynamic_##DTYPE##_property_t<STYPE>,        \
                                                 Graph, NamedParameters>                           \
{ };

CGAL_DEF_GET_INDEX_TYPE(Vertex, vertex, typename boost::graph_traits<Graph>::vertices_size_type)
CGAL_DEF_GET_INDEX_TYPE(Halfedge, halfedge, typename boost::graph_traits<Graph>::halfedges_size_type)
CGAL_DEF_GET_INDEX_TYPE(Edge, edge, typename boost::graph_traits<Graph>::edges_size_type)
CGAL_DEF_GET_INDEX_TYPE(Face, face, typename boost::graph_traits<Graph>::faces_size_type)

#undef CGAL_DEF_GET_INDEX_TYPE

// Define the following functions:
//
// get_initialized_vertex_index_map()
// get_initialized_halfedge_index_map()
// get_initialized_edge_index_map()
// get_initialized_face_index_map()
//
// The function returns:
// - the index property map passed in the NPs, if passed in the NPs; it must be initialized by the user;
// - the internal index property map if it is the graph has one. It is initialized if needed and possible;
// - an initialized dynamic pmap otherwise.

#define CGAL_DEF_GET_INITIALIZED_INDEX_MAP(DTYPE, STYPE)                                           \
template <typename Graph,                                                                          \
          typename NamedParameters = parameters::Default_named_parameters>                         \
typename BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                 \
                                               boost::DTYPE##_index_t,                             \
                                               CGAL::dynamic_##DTYPE##_property_t<STYPE>,          \
                                               Graph, NamedParameters>::const_type                 \
get_initialized_##DTYPE##_index_map(const Graph& g,                                                \
                                    const NamedParameters& np = parameters::default_values())      \
{                                                                                                  \
  typedef BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                \
                                                boost::DTYPE##_index_t,                            \
                                                CGAL::dynamic_##DTYPE##_property_t<STYPE>,         \
                                                Graph, NamedParameters>          Index_map_getter; \
  return Index_map_getter::get_const(CGAL::internal_np::DTYPE##_index_t{}, g, np);                 \
}                                                                                                  \
/* same as above, non-const version*/                                                              \
template <typename Graph,                                                                          \
          typename NamedParameters = parameters::Default_named_parameters,                         \
          /*otherwise compilers will try to use 'Graph := const PM' and things will go badly*/     \
          std::enable_if_t<                                                                        \
            !std::is_const<typename std::remove_reference<Graph>::type>::value, int> = 0>          \
typename BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                 \
                                               boost::DTYPE##_index_t,                             \
                                               CGAL::dynamic_##DTYPE##_property_t<STYPE>,          \
                                               Graph, NamedParameters>::type                       \
get_initialized_##DTYPE##_index_map(Graph& g,                                                      \
                                    const NamedParameters& np = parameters::default_values())      \
{                                                                                                  \
  typedef BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                \
                                                boost::DTYPE##_index_t,                            \
                                                CGAL::dynamic_##DTYPE##_property_t<STYPE>,         \
                                                Graph, NamedParameters>          Index_map_getter; \
  return Index_map_getter::get(CGAL::internal_np::DTYPE##_index_t{}, g, np);                       \
}                                                                                                  \

CGAL_DEF_GET_INITIALIZED_INDEX_MAP(vertex, typename boost::graph_traits<Graph>::vertices_size_type)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(halfedge, typename boost::graph_traits<Graph>::halfedges_size_type)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(edge, typename boost::graph_traits<Graph>::edges_size_type)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(face, typename boost::graph_traits<Graph>::faces_size_type)

#undef CGAL_DEF_GET_INITIALIZED_INDEX_MAP

namespace internal {

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_iterator, iterator, false)

} // namespace internal

template<typename PointRange,
         typename NamedParameters = parameters::Default_named_parameters,
         bool has_nested_iterator = internal::Has_nested_type_iterator<PointRange>::value,
         typename NP_TAG = internal_np::point_t>
class GetPointMap
{
  typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Point;
  typedef typename CGAL::Identity_property_map<Point> DefaultPMap;
  typedef typename CGAL::Identity_property_map<const Point> DefaultConstPMap;

public:
  typedef typename internal_np::Lookup_named_param_def<NP_TAG,
                                                       NamedParameters,
                                                       DefaultPMap>::type type;

  typedef typename internal_np::Lookup_named_param_def<NP_TAG,
                                                       NamedParameters,
                                                       DefaultConstPMap>::type const_type;
};

// to please compiler instantiating non valid overloads
template<typename PointRange, typename NamedParameters>
class GetPointMap<PointRange, NamedParameters, false>
{
  struct Dummy_point{};

public:
  typedef typename CGAL::Identity_property_map<Dummy_point> type;
  typedef typename CGAL::Identity_property_map<const Dummy_point> const_type;
};

template <typename PointRange, typename NamedParameters>
struct GetPolygonSoupGeomTraits
{
  typedef typename internal_np::Lookup_named_param_def <
                     internal_np::geom_traits_t,
                     NamedParameters,
                     typename CGAL::Kernel_traits<
                       typename boost::property_traits<
                         typename GetPointMap<PointRange, NamedParameters>::type
                       >::value_type
                     >::type
                   > ::type                                                         type;
};


template <class PointRange, class NamedParameters, class PointMap = Default, class NormalMap = Default>
struct Point_set_processing_3_np_helper
{
  typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Value_type;
  typedef typename Default::Get<PointMap, CGAL::Identity_property_map<Value_type>>::type DefaultPMap;
  typedef typename Default::Get<PointMap, CGAL::Identity_property_map<const Value_type>>::type DefaultConstPMap;

  typedef typename internal_np::Lookup_named_param_def<internal_np::point_t,
    NamedParameters,DefaultPMap> ::type  Point_map; // public
  typedef typename internal_np::Lookup_named_param_def<internal_np::point_t,
    NamedParameters,DefaultConstPMap> ::type  Const_point_map; // public

  typedef typename boost::property_traits<Point_map>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Default_geom_traits;

  typedef typename internal_np::Lookup_named_param_def <
      internal_np::geom_traits_t,
      NamedParameters,
      Default_geom_traits
    > ::type  Geom_traits; // public

  typedef typename Geom_traits::FT FT; // public

  typedef Constant_property_map<Value_type, typename Geom_traits::Vector_3> DummyNormalMap;
  typedef typename Default::Get<NormalMap, DummyNormalMap>::type DefaultNMap;

  typedef typename internal_np::Lookup_named_param_def<
    internal_np::normal_t,
    NamedParameters,
    DefaultNMap
    > ::type  Normal_map; // public

  static Point_map get_point_map(PointRange&, const NamedParameters& np)
  {
    return parameters::choose_parameter<Point_map>(parameters::get_parameter(np, internal_np::point_map));
  }

  static Point_map get_point_map(const NamedParameters& np)
  {
    return parameters::choose_parameter<Point_map>(parameters::get_parameter(np, internal_np::point_map));
  }

  static Const_point_map get_const_point_map(const PointRange&, const NamedParameters& np)
  {
    return parameters::choose_parameter<Const_point_map>(parameters::get_parameter(np, internal_np::point_map));
  }

  static Normal_map get_normal_map(const PointRange&, const NamedParameters& np)
  {
    return parameters::choose_parameter<Normal_map>(parameters::get_parameter(np, internal_np::normal_map));
  }

  static Normal_map get_normal_map(const NamedParameters& np)
  {
    return parameters::choose_parameter<Normal_map>(parameters::get_parameter(np, internal_np::normal_map));
  }

  static Geom_traits get_geom_traits(const PointRange&, const NamedParameters& np)
  {
    return parameters::choose_parameter<Geom_traits>(parameters::get_parameter(np, internal_np::geom_traits));
  }

  static constexpr bool has_normal_map(const PointRange&, const NamedParameters&)
  {
    using CGAL::parameters::is_default_parameter;
    return !(is_default_parameter<NamedParameters, internal_np::normal_t>::value);
  }
};

namespace Point_set_processing_3 {

template <typename ValueType>
struct Fake_point_range
{
  struct iterator
  {
    typedef ValueType value_type;
    typedef std::ptrdiff_t difference_type;
    typedef ValueType* pointer;
    typedef ValueType reference;
    typedef std::random_access_iterator_tag iterator_category;
  };
};

template<typename PlaneRange, typename NamedParameters>
class GetPlaneMap
{
  typedef typename PlaneRange::iterator::value_type Plane;
  typedef typename CGAL::Identity_property_map<Plane> DefaultPMap;
  typedef typename CGAL::Identity_property_map<const Plane> DefaultConstPMap;

public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::plane_t,
                                                       NamedParameters,
                                                       DefaultPMap>::type type;

  typedef typename internal_np::Lookup_named_param_def<internal_np::plane_t,
                                                       NamedParameters,
                                                       DefaultConstPMap>::type const_type;
};

template<typename NamedParameters>
class GetPlaneIndexMap
{
  typedef Constant_property_map<std::size_t, int> DummyPlaneIndexMap;

public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::plane_index_t,
                                                       NamedParameters,
                                                       DummyPlaneIndexMap>::type type;
};

template<typename PointRange, typename NamedParameters>
class GetIsConstrainedMap
{
  typedef Static_boolean_property_map<
    typename std::iterator_traits<typename PointRange::iterator>::value_type, false> Default_map;

public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::point_is_constrained_t,
                                                       NamedParameters,
                                                       Default_map>::type type;
};

template<typename PointRange, typename NamedParameters>
class GetAdjacencies
{
public:
  typedef Emptyset_iterator Empty;
  typedef typename internal_np::Lookup_named_param_def<internal_np::adjacencies_t,
                                                       NamedParameters,
                                                       Empty>::type type;
};

} // namespace Point_set_processing_3

template<typename NamedParameters, typename DefaultSolver>
class GetSolver
{
public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::sparse_linear_solver_t,
                                                       NamedParameters,
                                                       DefaultSolver>::type type;
};

template<typename NamedParameters, typename FT, unsigned int dim = 3>
class GetDiagonalizeTraits
{
public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::diagonalize_traits_t,
                                                       NamedParameters,
                                                       Default_diagonalize_traits<FT, dim> >::type type;
};

template<typename NamedParameters>
class GetSvdTraits
{
  struct DummySvdTraits
  {
    typedef double FT;
    typedef int Vector;
    typedef int Matrix;
    static FT solve (const Matrix&, Vector&) { return 0.; }
  };

public:
  typedef DummySvdTraits NoTraits;

  typedef typename internal_np::Lookup_named_param_def<internal_np::svd_traits_t,
                                                       NamedParameters,
#if defined(CGAL_EIGEN3_ENABLED)
                                                       Eigen_svd
#elif defined(CGAL_LAPACK_ENABLED)
                                                       Lapack_svd
#else
                                                       NoTraits
#endif
                                                       >::type type;
};

template<typename NamedParameters>
class GetImplementationTag
{
public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::implementation_tag_t,
                                                       NamedParameters,
                                                       Alpha_expansion_boost_adjacency_list_tag>::type type;
};

template<typename NP>
void set_stream_precision_from_NP(std::ostream& os, const NP& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;
  using parameters::is_default_parameter;

  if(!is_default_parameter<NP, internal_np::stream_precision_t>::value)
  {
    const int precision = choose_parameter<int>(get_parameter(np,
                            internal_np::stream_precision));
    os.precision(precision);
  }
}

} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_NAMED_PARAMETERS_HELPERS_H
