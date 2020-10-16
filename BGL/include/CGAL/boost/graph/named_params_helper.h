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
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/iterator.h>

#include <CGAL/property_map.h>

#include <boost/mpl/if.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/type_traits/is_same.hpp>

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
  class property_map_value {
    typedef typename boost::property_map<Graph, Property>::const_type PMap;
  public:
    typedef typename boost::property_traits<PMap>::value_type type;
  };

  template<typename PolygonMesh,
           typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
  class GetVertexPointMap
  {
    typedef typename property_map_selector<PolygonMesh, boost::vertex_point_t>::const_type
    DefaultVPMap_const;
    typedef typename property_map_selector<PolygonMesh, boost::vertex_point_t>::type
    DefaultVPMap;
  public:
    typedef typename internal_np::Lookup_named_param_def<
    internal_np::vertex_point_t,
    NamedParameters,
    DefaultVPMap
    > ::type  type;
    typedef typename internal_np::Lookup_named_param_def<
      internal_np::vertex_point_t,
      NamedParameters,
      DefaultVPMap_const
      > ::type  const_type;
  };

  namespace Polygon_mesh_processing {

  template<typename PolygonMesh, typename NamedParameters>
  class GetK
  {
    typedef typename boost::property_traits<
      typename GetVertexPointMap<PolygonMesh, NamedParameters>::type
      >::value_type Point;
  public:
    typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
  };

  } // namespace Polygon_mesh_processing

  template<typename PolygonMesh,
           typename NamedParametersGT = Named_function_parameters<bool, internal_np::all_default_t>,
           typename NamedParametersVPM = NamedParametersGT >
  class GetGeomTraits
  {
    typedef typename CGAL::graph_has_property<PolygonMesh, boost::vertex_point_t>::type
      Has_internal_pmap;

    typedef typename internal_np::Lookup_named_param_def <
      internal_np::vertex_point_t,
      NamedParametersVPM,
      internal_np::Param_not_found
    > ::type  NP_vpm;

    struct Fake_GT {};//to be used if there is no internal vertex_point_map in PolygonMesh

    typedef typename boost::mpl::if_c<Has_internal_pmap::value || !boost::is_same<internal_np::Param_not_found, NP_vpm>::value,
                                     typename Polygon_mesh_processing::GetK<PolygonMesh, NamedParametersVPM>::Kernel,
                                     Fake_GT>::type DefaultKernel;

  public:
    typedef typename internal_np::Lookup_named_param_def <
      internal_np::geom_traits_t,
      NamedParametersGT,
      DefaultKernel
    > ::type  type;
  };

// Define the following structs:
//
// GetInitializedVertexIndexMap
// GetInitializedHalfedgeIndexMap
// GetInitializedEdgeIndexMap
// GetInitializedFaceIndexMap

#define CGAL_DEF_GET_INDEX_TYPE(CTYPE, DTYPE, STYPE)                                               \
template <typename Graph,                                                                          \
          typename NamedParameters =                                                               \
            CGAL::Named_function_parameters<bool, CGAL::internal_np::all_default_t> >              \
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
          typename NamedParameters>                                                                \
typename BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                 \
                                               boost::DTYPE##_index_t,                             \
                                               CGAL::dynamic_##DTYPE##_property_t<STYPE>,          \
                                               Graph, NamedParameters>::const_type                 \
get_initialized_##DTYPE##_index_map(const Graph& g,                                                \
                                    const NamedParameters& np)                                     \
{                                                                                                  \
  typedef BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                \
                                                boost::DTYPE##_index_t,                            \
                                                CGAL::dynamic_##DTYPE##_property_t<STYPE>,         \
                                                Graph, NamedParameters>          Index_map_getter; \
  return Index_map_getter::get_const(CGAL::internal_np::DTYPE##_index_t{}, g, np);                 \
}                                                                                                  \
template <typename Graph>                                                                          \
typename BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                 \
                                               boost::DTYPE##_index_t,                             \
                                               CGAL::dynamic_##DTYPE##_property_t<STYPE>,          \
                                               Graph>::const_type                                  \
get_initialized_##DTYPE##_index_map(const Graph& g)                                                \
{                                                                                                  \
  return get_initialized_##DTYPE##_index_map(g, CGAL::parameters::all_default());                  \
}                                                                                                  \
/* same as above, non-const version*/                                                              \
template <typename Graph,                                                                          \
          typename NamedParameters,                                                                \
          /*otherwise compilers will try to use 'Graph := const PM' and things will go badly*/     \
          std::enable_if_t<                                                                        \
            !std::is_const<typename std::remove_reference<Graph>::type>::value, int> = 0>          \
typename BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                 \
                                               boost::DTYPE##_index_t,                             \
                                               CGAL::dynamic_##DTYPE##_property_t<STYPE>,          \
                                               Graph, NamedParameters>::type                       \
get_initialized_##DTYPE##_index_map(Graph& g,                                                      \
                                    const NamedParameters& np)                                     \
{                                                                                                  \
  typedef BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                \
                                                boost::DTYPE##_index_t,                            \
                                                CGAL::dynamic_##DTYPE##_property_t<STYPE>,         \
                                                Graph, NamedParameters>          Index_map_getter; \
  return Index_map_getter::get(CGAL::internal_np::DTYPE##_index_t{}, g, np);                       \
}                                                                                                  \
template <typename Graph,                                                                          \
          std::enable_if_t<                                                                        \
            !std::is_const<typename std::remove_reference<Graph>::type>::value, int> = 0>          \
typename BGL::internal::GetInitializedIndexMap<CGAL::internal_np::DTYPE##_index_t,                 \
                                               boost::DTYPE##_index_t,                             \
                                               CGAL::dynamic_##DTYPE##_property_t<STYPE>,          \
                                               Graph>::type                                        \
get_initialized_##DTYPE##_index_map(Graph& g)                                                      \
{                                                                                                  \
  return get_initialized_##DTYPE##_index_map(g, CGAL::parameters::all_default());                  \
}

CGAL_DEF_GET_INITIALIZED_INDEX_MAP(vertex, typename boost::graph_traits<Graph>::vertices_size_type)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(halfedge, typename boost::graph_traits<Graph>::halfedges_size_type)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(edge, typename boost::graph_traits<Graph>::edges_size_type)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(face, typename boost::graph_traits<Graph>::faces_size_type)

#undef CGAL_DEF_GET_INITIALIZED_INDEX_MAP

  template<typename PolygonMesh, typename NamedParameters>
  class GetFaceNormalMap
  {
    struct DummyNormalPmap
    {
      typedef typename boost::graph_traits<PolygonMesh>::face_descriptor key_type;
      typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Vector_3 value_type;
      typedef value_type reference;
      typedef boost::readable_property_map_tag category;

      typedef DummyNormalPmap Self;
      friend reference get(const Self&, const key_type&) { return CGAL::NULL_VECTOR; }
    };

  public:
    typedef DummyNormalPmap NoMap;
    typedef typename internal_np::Lookup_named_param_def <
      internal_np::face_normal_t,
      NamedParameters,
      DummyNormalPmap//default
      > ::type  type;
  };

  namespace internal {
    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_iterator, iterator, false)
  }

  template<typename PointRange,
           typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t>,
           bool has_nested_iterator = internal::Has_nested_type_iterator<PointRange>::value>
  class GetPointMap
  {
    typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Point;
    typedef typename CGAL::Identity_property_map<Point> DefaultPMap;

  public:
    typedef typename internal_np::Lookup_named_param_def<
    internal_np::point_t,
    NamedParameters,
    DefaultPMap
    > ::type  type;

    typedef typename internal_np::Lookup_named_param_def<
    internal_np::point_t,
    NamedParameters,
    DefaultPMap
    > ::type  const_type;
  };

  // to please compiler instantiating non valid overloads
  template<typename PointRange, typename NamedParameters>
  class GetPointMap<PointRange, NamedParameters, false>
  {
    struct Dummy_point{};
  public:
    typedef typename CGAL::Identity_property_map<Dummy_point> type;
    typedef typename CGAL::Identity_property_map<Dummy_point> const_type;
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

    namespace parameters
    {
      template <typename PointRange>
      Named_function_parameters<bool, internal_np::all_default_t>
      inline all_default(const PointRange&)
      {
        return CGAL::parameters::all_default();
      }
    }

    template<typename PointRange>
    class GetFT
    {
    public:
      typedef typename Kernel_traits<
        typename std::iterator_traits<
          typename PointRange::iterator
          >::value_type
        >::Kernel::FT type;
    };

    template<typename PointRange, typename NamedParameters>
    class GetQueryPointMap
    {
      typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Point;
      typedef typename CGAL::Identity_property_map<Point> DefaultPMap;

    public:
      typedef typename internal_np::Lookup_named_param_def<
      internal_np::query_point_t,
      NamedParameters,
      DefaultPMap
      > ::type  type;

      typedef typename internal_np::Lookup_named_param_def<
      internal_np::query_point_t,
      NamedParameters,
      DefaultPMap
      > ::type  const_type;
    };

    template<typename PointRange, typename NamedParameters>
    class GetK
    {
      typedef typename GetPointMap<PointRange, NamedParameters>::type Vpm;
      typedef typename Kernel_traits<
        typename boost::property_traits<Vpm>::value_type
      >::Kernel Default_kernel;

    public:
      typedef typename internal_np::Lookup_named_param_def <
        internal_np::geom_traits_t,
        NamedParameters,
        Default_kernel
      > ::type  Kernel;
    };

    template<typename PointRange, typename NamedParameters>
    class GetNormalMap
    {
      struct DummyNormalMap
      {
        typedef typename std::iterator_traits<typename PointRange::iterator>::value_type key_type;
        typedef typename GetK<PointRange, NamedParameters>::Kernel::Vector_3 value_type;
        typedef value_type reference;
        typedef boost::read_write_property_map_tag category;

        typedef DummyNormalMap Self;
        friend reference get(const Self&, const key_type&) { return CGAL::NULL_VECTOR; }
        friend void put(const Self&, const key_type&, const value_type&) { }
      };

    public:
      typedef DummyNormalMap NoMap;
      typedef typename internal_np::Lookup_named_param_def <
        internal_np::normal_t,
        NamedParameters,
        DummyNormalMap//default
        > ::type  type;
    };

    template<typename PlaneRange, typename NamedParameters>
    class GetPlaneMap
    {
      typedef typename PlaneRange::iterator::value_type Plane;
      typedef typename CGAL::Identity_property_map<Plane> DefaultPMap;

    public:
      typedef typename internal_np::Lookup_named_param_def<
      internal_np::plane_t,
      NamedParameters,
      DefaultPMap
      > ::type  type;

      typedef typename internal_np::Lookup_named_param_def<
      internal_np::plane_t,
      NamedParameters,
      DefaultPMap
      > ::type  const_type;
    };

    template<typename NamedParameters>
    class GetPlaneIndexMap
    {
      struct DummyPlaneIndexMap
      {
        typedef std::size_t key_type;
        typedef int value_type;
        typedef value_type reference;
        typedef boost::readable_property_map_tag category;

        typedef DummyPlaneIndexMap Self;
        friend reference get(const Self&, const key_type&) { return -1; }
      };

    public:
      typedef DummyPlaneIndexMap NoMap;
      typedef typename internal_np::Lookup_named_param_def <
        internal_np::plane_index_t,
        NamedParameters,
        DummyPlaneIndexMap//default
        > ::type  type;
    };

    template<typename PointRange, typename NamedParameters>
    class GetIsConstrainedMap
    {
      struct DummyConstrainedMap
      {
        typedef typename std::iterator_traits<typename PointRange::iterator>::value_type key_type;
        typedef bool value_type;
        typedef value_type reference;
        typedef boost::readable_property_map_tag category;

        typedef DummyConstrainedMap Self;
        friend reference get(const Self&, const key_type&) { return false; }
      };

    public:
      typedef DummyConstrainedMap NoMap;
      typedef typename internal_np::Lookup_named_param_def <
        internal_np::point_is_constrained_t,
        NamedParameters,
        DummyConstrainedMap //default
        > ::type  type;
    };

    template<typename PointRange, typename NamedParameters>
    class GetAdjacencies
    {
    public:
      typedef Emptyset_iterator Empty;
      typedef typename internal_np::Lookup_named_param_def <
        internal_np::adjacencies_t,
        NamedParameters,
        Empty//default
        > ::type  type;
    };

  } // namespace Point_set_processing_3

  template<typename NamedParameters, typename DefaultSolver>
  class GetSolver
  {
  public:
    typedef typename internal_np::Lookup_named_param_def <
    internal_np::sparse_linear_solver_t,
    NamedParameters,
    DefaultSolver
    > ::type type;
  };

  template<typename NamedParameters, typename FT, unsigned int dim = 3>
  class GetDiagonalizeTraits
  {
  public:
    typedef typename internal_np::Lookup_named_param_def <
    internal_np::diagonalize_traits_t,
    NamedParameters,
    Default_diagonalize_traits<FT, dim>
    > ::type type;
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

    typedef typename internal_np::Lookup_named_param_def <
    internal_np::svd_traits_t,
    NamedParameters,
#if defined(CGAL_EIGEN3_ENABLED)
    Eigen_svd
#elif defined(CGAL_LAPACK_ENABLED)
    Lapack_svd
#else
    NoTraits
#endif
    > ::type type;
  };

  template<typename NamedParameters>
  class GetImplementationTag
  {
  public:
    typedef typename internal_np::Lookup_named_param_def <
    internal_np::implementation_tag_t,
    NamedParameters,
    Alpha_expansion_boost_adjacency_list_tag
    >::type type;
  };
} //namespace CGAL


#endif // CGAL_BOOST_GRAPH_NAMED_PARAMETERS_HELPERS_H
