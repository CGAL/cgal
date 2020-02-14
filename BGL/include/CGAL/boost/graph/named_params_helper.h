//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// This file is part of the Boost Graph Library
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// https://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
// Copyright (c) 2007-2015  GeometryFactory (France).  All rights reserved.
//
// $URL$
// $Id$
// SPDX-License-Identifier: BSL-1.0
//
// Author(s)     : Andreas Fabri, Fernando Cacciola, Jane Tournois

#ifndef CGAL_BOOST_GRAPH_NAMED_PARAMETERS_HELPERS_H
#define CGAL_BOOST_GRAPH_NAMED_PARAMETERS_HELPERS_H

#include <CGAL/boost/graph/Named_function_parameters.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>

#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/mpl/if.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <CGAL/Dynamic_property_map.h>

#include <boost/type_traits/is_same.hpp>

#include <type_traits>

namespace CGAL {

  // forward declarations to avoid dependency to Solver_interface
  template <typename FT, unsigned int dim>
  class Default_diagonalize_traits;
  class Eigen_svd;
  class Lapack_svd;
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
// shortcut for accessing the value type of the property map
  template <class Graph, class Property>
  class property_map_value {
    typedef typename boost::property_map<Graph, Property>::const_type PMap;
  public:
    typedef typename boost::property_traits<PMap>::value_type type;
  };

  namespace internal_np
  {
  template<typename key_value, class PMap, class Graph, typename Tag>
  struct Index_map_initializer{
    void operator()(PMap, const Graph& )
    {}
  };

  template< class PMap, class Graph>
  struct Index_map_initializer<
      typename boost::graph_traits<Graph>::vertex_descriptor,
      PMap, Graph,
      CGAL::Tag_true>{
    void operator()(PMap map, const Graph& g)
    {
      CGAL::helpers::init_vertex_indices(g, map);
    }
  };

  template< class PMap, class Graph>
  struct Index_map_initializer<
      typename boost::graph_traits<Graph>::halfedge_descriptor,
      PMap, Graph,
      CGAL::Tag_true>{
    void operator()(PMap map, const Graph& g)
    {
      CGAL::helpers::init_halfedge_indices(g, map);
    }
  };

  template< class PMap, class Graph>
  struct Index_map_initializer<
      typename boost::graph_traits<Graph>::face_descriptor,
      PMap, Graph,
      CGAL::Tag_true>{
    void operator()(PMap map,const Graph& g)
    {
      CGAL::helpers::init_face_indices(g, map);
    }
  };

#define CGAL_IS_PMAP_WRITABLE(TAG) template<typename IsRefConst> \
  struct Is_pmap_writable<TAG, IsRefConst>{ typedef CGAL::Tag_true result; }; \
  template<> \
  struct Is_pmap_writable<TAG,std::integral_constant<bool, true> >{ typedef CGAL::Tag_false result; };

  template<typename PMapCategory, typename IsRefConst>
  struct Is_pmap_writable{
    typedef CGAL::Tag_false result;
  };

  CGAL_IS_PMAP_WRITABLE(boost::read_write_property_map_tag)
  CGAL_IS_PMAP_WRITABLE(boost::writable_property_map_tag)
  CGAL_IS_PMAP_WRITABLE(boost::lvalue_property_map_tag)
#undef CGAL_IS_PMAP_WRITABLE

  //overloads used to select a default map:
  // use the one passed in the named parameters (user must have initialized it)
  template <class MapFromNP, class Default_tag, class Dynamic_tag, class Mesh>
  MapFromNP
  get_ndi_map(MapFromNP m, Default_tag, Dynamic_tag, const Mesh&)
  {
    return m;
  }

  // use the one internal to the mesh (it will be init if writable)
  template <class Default_tag, class Dynamic_tag, class Mesh>
  typename boost::property_map<Mesh, Default_tag >::const_type
  get_ndi_map(CGAL::internal_np::Param_not_found, Default_tag t, Dynamic_tag , const Mesh& m)
  {
    typename boost::property_map<Mesh, Default_tag >::const_type map = get(t, m);
    Index_map_initializer<
        typename boost::property_traits<typename boost::property_map<Mesh, Default_tag >::const_type>::key_type,
        typename boost::property_map<Mesh, Default_tag >::const_type,
        Mesh,
        typename Is_pmap_writable<
        typename boost::property_traits
        <typename boost::property_map<Mesh, Default_tag >
        ::const_type>::category,
        typename boost::property_traits
        <typename boost::property_map<Mesh, Default_tag >
        ::const_type>::reference
        >::result>
        ()(map, m);
    return map;
  }

  // create a dynamic property and initialize it
  template <class Dynamic_tag, class Mesh>
  typename boost::property_map<Mesh, Dynamic_tag >::const_type
  get_ndi_map(CGAL::internal_np::Param_not_found, Dynamic_tag t, Dynamic_tag , const Mesh& m)
  {
    typename boost::property_map<Mesh, Dynamic_tag >::const_type map = get(t,m);
    Index_map_initializer<
        typename boost::property_traits<typename boost::property_map<Mesh, Dynamic_tag >::const_type>::key_type,
        typename boost::property_map<Mesh, Dynamic_tag >::const_type,
        Mesh,
        CGAL::Tag_true>()(map, m);
    return map;
  }

  }//end of internal_np

  namespace Polygon_mesh_processing
  {

  //define types for maps :
  //struct Default_face_index_map
  //struct Default_vertex_index_map
  //struct Default_halfedge_index_map
#define CGAL_DEF_MAP_TYPE(TYPE)                                    \
  template<typename NP, typename TM>                          \
  struct Default_##TYPE##_index_map{                          \
  typedef typename boost::mpl::if_c<                          \
  CGAL::graph_has_property<TM, boost::TYPE##_index_t>::value  \
  , boost::TYPE##_index_t                                     \
  , CGAL::dynamic_##TYPE##_property_t<int>                    \
  >::type Final_tag;                                          \
  typedef typename internal_np::Lookup_named_param_def<       \
  internal_np::TYPE##_index_t,                                  \
  NP,                                                         \
  typename boost::property_map<TM, Final_tag >::const_type    \
  > ::type  type;                                             \
  };

  CGAL_DEF_MAP_TYPE(face)
  CGAL_DEF_MAP_TYPE(vertex)
  CGAL_DEF_MAP_TYPE(halfedge)
#undef CGAL_DEF_MAP_TYPE


  template<typename Tag, typename Dynamic_tag, typename Mesh,
  typename NamedParameters, typename Parameter>
  class Get_index_map_from_NP {
  private :
    const Dynamic_tag dtag;
    const Mesh& m;
    const NamedParameters& np;
    const Parameter p;

  public:
    //get the Default tag :
    //if Mesh has an internal property map for Tag, use Tag, else use the Dynamic_tag.
    typedef typename boost::mpl::if_c<CGAL::graph_has_property<Mesh, Tag>::value
    , Tag
    , Dynamic_tag
    >::type Final_tag;

    //If Parameter is in NamedParameters, take the NP map.
    //Else, take the default map.
    typedef typename internal_np::Lookup_named_param_def<
    Parameter,
    NamedParameters,
    typename boost::property_map<Mesh, Final_tag >::const_type
    > ::type  PropertyMapType;


    Get_index_map_from_NP(const Tag,
                 const Dynamic_tag dtag,
                 const Mesh& m,
                 const NamedParameters& np,
                 const Parameter p)
      : dtag(dtag), m(m), np(np), p(p) {}


    PropertyMapType property_map()
    {
      return internal_np::get_ndi_map(
            parameters::get_parameter(np, p),
            Final_tag(),
            dtag,
            m);
    }
  };

  //define the
  // get_initialized_face_index_map(), get_initialized_vertex_index_map(), get_initialized_halfedge_index_map()
  // functions.
  //This comment is here to make it easier to find the definition of the functions with a grep.

#define CGAL_DEF_GET_INIT_ID_MAP(TYPE) template<class PolygonMesh, class NamedParameters> \
  typename Default_##TYPE##_index_map<NamedParameters, PolygonMesh>::type            \
  get_initialized_##TYPE##_index_map(const PolygonMesh& pmesh, const NamedParameters& np){ \
  typedef Get_index_map_from_NP<boost::TYPE##_index_t,                               \
  CGAL::dynamic_##TYPE##_property_t<int>,                                            \
  PolygonMesh, NamedParameters, internal_np::TYPE##_index_t> MapGetter;              \
  MapGetter get_map(boost::TYPE##_index_t(),                                         \
  CGAL::dynamic_##TYPE##_property_t<int>(),                                          \
  pmesh, np, internal_np::TYPE##_index);                                             \
  return get_map.property_map();                                                     \
  }
  CGAL_DEF_GET_INIT_ID_MAP(face)
  CGAL_DEF_GET_INIT_ID_MAP(vertex)
  CGAL_DEF_GET_INIT_ID_MAP(halfedge)

#undef CGAL_DEF_GET_INIT_ID_MAP
} //end Polygon_mesh_processing
  
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

  template<typename PolygonMesh, typename NamedParameters>
  class GetFaceIndexMap
  {
    typedef typename property_map_selector<PolygonMesh, boost::face_index_t>::type DefaultMap;
    typedef typename property_map_selector<PolygonMesh, boost::face_index_t>::const_type DefaultMap_const;
  public:
    typedef typename internal_np::Lookup_named_param_def <
    internal_np::face_index_t,
    NamedParameters,
    DefaultMap
    > ::type  type;
    typedef typename internal_np::Lookup_named_param_def <
      internal_np::face_index_t,
      NamedParameters,
      DefaultMap_const
      > ::type  const_type;
    typedef typename boost::is_same<type, DefaultMap>::type Is_internal_map;
    typedef typename boost::is_same<const_type, DefaultMap_const>::type Is_internal_map_const;
  };

  template<typename PolygonMesh, typename NamedParameters>
  class GetVertexIndexMap
  {
    typedef typename property_map_selector<PolygonMesh, boost::vertex_index_t>::type DefaultMap;
  public:
    typedef typename internal_np::Lookup_named_param_def <
    internal_np::vertex_index_t,
    NamedParameters,
    DefaultMap
    > ::type  type;
  };

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
  

  namespace Point_set_processing_3
  {
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

    namespace internal{
      BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_iterator, iterator, false)
    }
    
    template<typename PointRange, typename NamedParameters,
             bool has_nested_iterator=internal::Has_nested_type_iterator<PointRange>::value>
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
  
} //namespace CGAL


#endif // CGAL_BOOST_GRAPH_NAMED_PARAMETERS_HELPERS_H
