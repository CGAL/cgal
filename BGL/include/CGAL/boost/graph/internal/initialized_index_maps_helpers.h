// Copyright (c) 2020 GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
// 
// Author(s) : Mael Rouxel-Labbé
//             Maxime Gimeno

#ifndef CGAL_BOOST_GRAPH_INITIALIZED_INTERNAL_INDEX_MAPS_HELPERS
#define CGAL_BOOST_GRAPH_INITIALIZED_INTERNAL_INDEX_MAPS_HELPERS

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/use.h>

#include <vector>

namespace CGAL {
namespace BGL {
namespace internal {

// Check that an index map has been correctly initialized
template <typename DescriptorRange, typename IndexMap>
bool is_index_map_valid(const IndexMap idmap,
                        const std::size_t num_simplices,
                        const DescriptorRange& range)
{
  typedef typename boost::property_traits<IndexMap>::value_type Id_type;

  Id_type max_id = static_cast<Id_type>(num_simplices);
  std::vector<bool> indices(max_id);
  for(const auto& d : range)
  {
    const Id_type id = get(idmap, d);
    if(id >= 0 && id < max_id && !indices[id])
    {
      indices[id] = true;
    }
    else
    {
      std::cerr << "Invalid ID: " << id << " num_simplices: " << num_simplices << std::endl;
      return false;
    }
  }

  return true;
}

template <typename VertexIndexPropertyMap, typename Graph>
bool is_index_map_valid(const CGAL::internal_np::vertex_index_t, VertexIndexPropertyMap vertex_index_map, const Graph& g)
{
  return is_index_map_valid(vertex_index_map, num_vertices(g), vertices(g));
}

template <typename HalfedgeIndexPropertyMap, typename Graph>
bool is_index_map_valid(const CGAL::internal_np::halfedge_index_t, HalfedgeIndexPropertyMap halfedge_index_map, const Graph& g)
{
  return is_index_map_valid(halfedge_index_map, num_halfedges(g), halfedges(g));
}

template <typename EdgeIndexPropertyMap, typename Graph>
bool is_index_map_valid(const CGAL::internal_np::edge_index_t, EdgeIndexPropertyMap edge_index_map, const Graph& g)
{
  return is_index_map_valid(edge_index_map, num_edges(g), edges(g));
}

template <typename FaceIndexPropertyMap, typename Graph>
bool is_index_map_valid(const CGAL::internal_np::face_index_t, FaceIndexPropertyMap face_index_map, const Graph& g)
{
  return is_index_map_valid(face_index_map, num_faces(g), faces(g));
}

template <typename Parameter, typename IndexPropertyMap, typename Graph>
void initialize_index_map(const Parameter, IndexPropertyMap, const Graph&)
{
  // Unknown parameter; should never be here.
  CGAL_assertion(false);
}

template <typename IndexPropertyMap,
          typename Graph,
          bool is_writable = CGAL::internal::Is_writable_property_map<IndexPropertyMap>::value>
struct Index_map_initializer
{
  void operator()(const CGAL::internal_np::vertex_index_t, IndexPropertyMap vertex_index_map, const Graph& g)
  {
    typename boost::property_traits<IndexPropertyMap>::value_type i = 0;
    for(typename boost::graph_traits<Graph>::vertex_descriptor vd : vertices(g))
      put(vertex_index_map, vd, i++);
  }

  void operator()(const CGAL::internal_np::halfedge_index_t, IndexPropertyMap halfedge_index_map, const Graph& g)
  {
    typename boost::property_traits<IndexPropertyMap>::value_type i = 0;
    for(typename boost::graph_traits<Graph>::halfedge_descriptor hd : halfedges(g))
      put(halfedge_index_map, hd, i++);
  }

  void operator()(const CGAL::internal_np::edge_index_t, IndexPropertyMap edge_index_map, const Graph& g)
  {
    typename boost::property_traits<IndexPropertyMap>::value_type i = 0;
    for(typename boost::graph_traits<Graph>::edge_descriptor ed : edges(g))
      put(edge_index_map, ed, i++);
  }

  void operator()(const CGAL::internal_np::face_index_t, IndexPropertyMap face_index_map, const Graph& g)
  {
    typename boost::property_traits<IndexPropertyMap>::value_type i = 0;
    for(typename boost::graph_traits<Graph>::face_descriptor fd : faces(g))
      put(face_index_map, fd, i++);
  }

  template <typename Parameter>
  void operator()(const Parameter, IndexPropertyMap, const Graph&)
  {
    // Unknown parameter; should never be here.
    CGAL_assertion(false);
  }
};

template <typename IndexPropertyMap, typename Graph>
struct Index_map_initializer<IndexPropertyMap, Graph, false>
{
  template <typename Parameter>
  void operator()(const Parameter, IndexPropertyMap, const Graph&)
  {
    // The property map is not writable; should never be here.
    CGAL_assertion(false);
  }
};

// Just for convenience, define the following functions:
//
// BGL::internal::initialize_vertex_index_map()
// BGL::internal::initialize_halfedge_index_map()
// BGL::internal::initialize_edge_index_map()
// BGL::internal::initialize_face_index_map()

#define CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(TYPE)                                                  \
template <typename WritableIndexPropertyMap, typename Graph>                                       \
void initialize_##TYPE##_index_map(WritableIndexPropertyMap index_map,                             \
                                   const Graph& g)                                                 \
{                                                                                                  \
  Index_map_initializer<WritableIndexPropertyMap, Graph> initializer;                              \
  initializer(CGAL::internal_np::TYPE##_index_t(), index_map, g);                                  \
}

CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(vertex)
CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(halfedge)
CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(edge)
CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(face)

#undef CGAL_DEF_INITIALIZE_ID_FUCNTION

// Using the pmap passed in named parameters
template <typename IndexMap, typename Parameter, typename Tag, typename DynamicTag, typename Graph>
IndexMap get_initialized_index_map(const IndexMap index_map,
                                   const Parameter p, Tag, DynamicTag,
                                   const Graph& g)
{
  CGAL_USE(g);
  CGAL_USE(p);
  CGAL_assertion(is_index_map_valid(p, index_map, g));

  return index_map;
}

// Using the internal to the mesh
template <typename Parameter, typename Tag, typename DynamicTag, typename Graph>
typename boost::property_map<Graph, Tag>::const_type
get_initialized_index_map(CGAL::internal_np::Param_not_found,
                          const Parameter p, const Tag tag, DynamicTag,
                          const Graph& g) // @todo non-const
{
  typedef typename boost::property_map<Graph, Tag>::const_type          Index_map;
  Index_map index_map = get(tag, g);

  if(CGAL::internal::Is_writable_property_map<Index_map>::value)
  {
    if(!is_index_map_valid(p, index_map, g))
      Index_map_initializer<Index_map, Graph>{}(p, index_map, g);
  }
  else // not writable
  {
    CGAL_assertion(is_index_map_valid(p, index_map, g));
  }

  return index_map;
}

// Create a dynamic property and initialize it
template <typename Parameter, typename DynamicTag, typename Graph>
typename boost::property_map<Graph, DynamicTag>::const_type
get_initialized_index_map(CGAL::internal_np::Param_not_found,
                          const Parameter p, const DynamicTag tag, DynamicTag,
                          const Graph& g)
{
  typedef typename boost::property_map<Graph, DynamicTag>::const_type   Index_map;

  Index_map index_map = get(tag, g);
  Index_map_initializer<Index_map, Graph>{}(p, index_map, g);

  return index_map;
}

template <typename Parameter, typename Tag, typename DynamicTag,
         typename Graph,
         typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
class GetInitializedIndexMap
{
public:
  // Definition of the Tag that will be used if there is no named parameter
  typedef typename boost::mpl::if_c<
      CGAL::graph_has_property<Graph, Tag>::value, Tag, DynamicTag>::type    Final_tag;

  typedef typename internal_np::Lookup_named_param_def<
      Parameter,
      NamedParameters,
      typename boost::property_map<Graph, Final_tag>::type>::type            type;

  typedef typename internal_np::Lookup_named_param_def<
      Parameter,
      NamedParameters,
      typename boost::property_map<Graph, Final_tag>::const_type>::type      const_type;

  static const_type get(const Parameter p, const Graph& g, const NamedParameters& np)
  {
    return BGL::internal::get_initialized_index_map(parameters::get_parameter(np, p),
                                                    p, Final_tag(), DynamicTag(), g);
  }

  static type get(const Parameter p, Graph& g, const NamedParameters& np)
  {
    return BGL::internal::get_initialized_index_map(parameters::get_parameter(np, p),
                                                    p, Final_tag(), DynamicTag(), g);
  }
};

} // namespace internal
} // namespace BGL

// @todo move below to named_params_...

#define CGAL_DEF_GET_INDEX_TYPE(CTYPE, TYPE)                                                       \
template <typename Graph,                                                                          \
          typename NamedParameters =                                                               \
            CGAL::Named_function_parameters<bool, CGAL::internal_np::all_default_t> >              \
struct GetInitialized##CTYPE##IndexMap                                                             \
  : public BGL::internal::GetInitializedIndexMap<internal_np::TYPE##_index_t,                      \
                                                 boost::TYPE##_index_t,                            \
                                                 CGAL::dynamic_##TYPE##_property_t<int>,           \
                                                 Graph, NamedParameters>                           \
{ };

CGAL_DEF_GET_INDEX_TYPE(Vertex, vertex)
CGAL_DEF_GET_INDEX_TYPE(Halfedge, halfedge)
CGAL_DEF_GET_INDEX_TYPE(Edge, edge)
CGAL_DEF_GET_INDEX_TYPE(Face, face)

#undef CGAL_DEF_GET_INDEX_TYPE

// @todo move below to properties.h

// Define the following functions:
//
// get_initialized_vertex_index_map();
// get_initialized_halfedge_index_map();
// get_initialized_edge_index_map();
// get_initialized_face_index_map()
//
// The function returns:
// - the index property map passed in the NPs, if passed in the NPs; it must be initialized by the user;
// - the internal index property map if it is the graph has one. It is initialized if needed and possible;
// - an initialized dynamic pmap otherwise.

#define CGAL_DEF_GET_INITIALIZED_INDEX_MAP(TYPE)                                                   \
template <typename Graph,                                                                          \
          typename NamedParameters>                                                                \
typename BGL::internal::GetInitializedIndexMap<CGAL::internal_np::TYPE##_index_t,                  \
                                               boost::TYPE##_index_t,                              \
                                               CGAL::dynamic_##TYPE##_property_t<int>,             \
                                               Graph, NamedParameters>::const_type                 \
get_initialized_##TYPE##_index_map(const Graph& g,                                                 \
                                   const NamedParameters& np)                                      \
{                                                                                                  \
  typedef BGL::internal::GetInitializedIndexMap<CGAL::internal_np::TYPE##_index_t,                 \
                                                boost::TYPE##_index_t,                             \
                                                CGAL::dynamic_##TYPE##_property_t<int>,            \
                                                Graph, NamedParameters>          Index_map_getter; \
  return Index_map_getter::get(CGAL::internal_np::TYPE##_index_t(), g, np);                        \
}                                                                                                  \
template <typename Graph>                                                                          \
typename BGL::internal::GetInitializedIndexMap<CGAL::internal_np::TYPE##_index_t,                  \
                                               boost::TYPE##_index_t,                              \
                                               CGAL::dynamic_##TYPE##_property_t<int>,             \
                                               Graph>::const_type                                  \
get_initialized_##TYPE##_index_map(const Graph& g)                                                 \
{                                                                                                  \
  return get_initialized_##TYPE##_index_map(g, CGAL::parameters::all_default());                   \
}

// @todo add the non-const Graph& version

CGAL_DEF_GET_INITIALIZED_INDEX_MAP(vertex)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(halfedge)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(edge)
CGAL_DEF_GET_INITIALIZED_INDEX_MAP(face)

#undef CGAL_DEF_GET_INITIALIZED_INDEX_MAP

} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_INITIALIZED_INTERNAL_INDEX_MAPS_HELPERS
