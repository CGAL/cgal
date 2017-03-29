// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Jane Tournois

#ifndef CGAL_NAMED_PARAMETERS_HELPERS_H
#define CGAL_NAMED_PARAMETERS_HELPERS_H

#include <CGAL/license/Polygon_mesh_processing/core.h>


#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>

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
    boost::vertex_point_t,
    NamedParameters,
    DefaultVPMap
  > ::type  type;
  typedef typename boost::lookup_named_param_def<
    boost::vertex_point_t,
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

template<typename PolygonMesh, typename NamedParameters = pmp_bgl_named_params<bool, all_default_t> >
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
    CGAL::geom_traits_t,
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

template<typename PolygonMesh, typename NamedParameters>
class GetFaceIndexMap
{
  typedef typename property_map_selector<PolygonMesh, boost::face_index_t>::type DefaultMap;
  typedef typename property_map_selector<PolygonMesh, boost::face_index_t>::const_type DefaultMap_const;
public:
  typedef typename boost::lookup_named_param_def <
    boost::face_index_t,
    NamedParameters,
    DefaultMap
  > ::type  type;
  typedef typename boost::lookup_named_param_def <
    boost::face_index_t,
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
  typedef typename boost::lookup_named_param_def <
    boost::vertex_index_t,
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
  typedef typename boost::lookup_named_param_def <
    CGAL::face_normal_t,
    NamedParameters,
    DummyNormalPmap//default
  > ::type  type;
};

template<typename NamedParameters, typename DefaultSolver>
class GetSolver
{
public:
  typedef typename boost::lookup_named_param_def <
    CGAL::sparse_linear_solver_t,
    NamedParameters,
    DefaultSolver
  > ::type type;
};

} //end of namespace CGAL

#endif //CGAL_NAMED_PARAMETERS_HELPERS_H
