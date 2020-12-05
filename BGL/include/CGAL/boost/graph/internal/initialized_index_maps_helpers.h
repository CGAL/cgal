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
#include <type_traits>

namespace CGAL {
namespace BGL {
namespace internal {

// Check that an index map has been correctly initialized
template <typename DescriptorRange, typename IndexMap>
bool is_index_map_valid(IndexMap idmap,
                        const std::size_t num_simplices,
                        const DescriptorRange& range)
{
  typedef typename boost::property_traits<IndexMap>::value_type Id_type;

  Id_type max_id = static_cast<Id_type>(num_simplices);
  std::vector<bool> indices(max_id);

  // According to concepts, the descriptor ranges such as 'vertices(g)' return a 'std::pair<it, it>'
  for(auto it = range.first; it != range.second; ++it)
  {
    const Id_type id = get(idmap, *it);
    if(id >= 0 && id < max_id && !indices[id])
    {
      indices[id] = true;
    }
    else
    {
#ifdef CGAL_BGL_INDEX_MAP_DEBUG
      std::cerr << "Invalid ID: " << id << " num_simplices: " << num_simplices << std::endl;
#endif
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

template <typename PropertyTag, typename IndexPropertyMap, typename Graph>
void initialize_index_map(const PropertyTag, IndexPropertyMap, const Graph&)
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

  template <typename PropertyTag>
  void operator()(const PropertyTag, IndexPropertyMap, const Graph&)
  {
    // Unknown parameter; should never be here.
    CGAL_assertion(false);
  }
};

template <typename IndexPropertyMap, typename Graph>
struct Index_map_initializer<IndexPropertyMap, Graph, false>
{
  template <typename PropertyTag>
  void operator()(const PropertyTag, IndexPropertyMap, const Graph&)
  {
    // The property map is not writable; should never be here.
    CGAL_assertion_msg(false, "Initialization of a non-writable property map is impossible");
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
  initializer(CGAL::internal_np::TYPE##_index_t{}, index_map, g);                                  \
}

CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(vertex)
CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(halfedge)
CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(edge)
CGAL_DEF_INITIALIZE_ID_MAP_FUNCTION(face)

#undef CGAL_DEF_INITIALIZE_ID_FUCNTION

// Using the pmap passed in named parameters -------------------------------------------------------
template <typename IndexMap, typename PropertyTag, typename Tag, typename DynamicTag, typename Graph>
IndexMap get_initialized_index_map_const(const IndexMap index_map,
                                         const PropertyTag p, Tag, DynamicTag,
                                         const Graph& g)
{
  CGAL_USE(g);
  CGAL_USE(p);

  // If a pmap is passed via NPs, it must be initialized
  CGAL_assertion(is_index_map_valid(p, index_map, g));

  return index_map;
}

template <typename IndexMap, typename PropertyTag, typename Tag, typename DynamicTag, typename Graph>
IndexMap get_initialized_index_map(const IndexMap index_map,
                                   const PropertyTag p, Tag, DynamicTag,
                                   Graph& g)
{
  CGAL_USE(g);
  CGAL_USE(p);

  // If a pmap is passed via NPs, it must be initialized
  CGAL_assertion(is_index_map_valid(p, index_map, g));

  return index_map;
}

// Using the internal to the mesh ------------------------------------------------------------------
template <typename InternalIndexMap, typename PropertyTag, typename Graph>
InternalIndexMap
get_initialized_internal_index_map(InternalIndexMap index_map,
                                   const PropertyTag p,
                                   const Graph& g)
{
  if(CGAL::internal::Is_writable_property_map<InternalIndexMap>::value)
  {
    if(!is_index_map_valid(p, index_map, g))
      Index_map_initializer<InternalIndexMap, Graph>{}(p, index_map, g);
  }
  else // not writable
  {
    CGAL_assertion(is_index_map_valid(p, index_map, g));
  }

  return index_map;
}

template <typename PropertyTag, typename Tag, typename DynamicTag, typename Graph>
typename boost::property_map<Graph, Tag>::const_type
get_initialized_index_map_const(CGAL::internal_np::Param_not_found,
                                const PropertyTag p, const Tag tag, DynamicTag,
                                const Graph& g)
{
  return get_initialized_internal_index_map(get(tag, g), p, g);
}

// same as above, non-const graph overload
template <typename PropertyTag, typename Tag, typename DynamicTag, typename Graph>
typename boost::property_map<Graph, Tag>::type
get_initialized_index_map(CGAL::internal_np::Param_not_found,
                          const PropertyTag p, const Tag tag, DynamicTag,
                          Graph& g)
{
  // From now on the correct property map has been acquired
  // and there is no need to distinguish between const and non-const mesh
  return get_initialized_internal_index_map(get(tag, g), p, g);
}

// Create a dynamic property and initialize it -----------------------------------------------------
template <typename DynamicIndexMap, typename PropertyTag, typename Graph>
DynamicIndexMap
get_initialized_dynamic_index_map(DynamicIndexMap index_map,
                                  const PropertyTag p,
                                  const Graph& g)
{
#ifdef CGAL_PERFORMANCE_WARNINGS
  std::cerr << "Warning: the automatically selected index map is a dynamic property map,"
            << " which might not have constant-time access complexity." << std::endl;
#endif

  Index_map_initializer<DynamicIndexMap, Graph>{}(p, index_map, g);
  return index_map;
}

template <typename PropertyTag, typename DynamicTag, typename Graph>
typename boost::property_map<Graph, DynamicTag>::const_type
get_initialized_index_map_const(CGAL::internal_np::Param_not_found,
                                const PropertyTag p, const DynamicTag dtag, DynamicTag,
                                const Graph& g)
{
  return get_initialized_dynamic_index_map(get(dtag, g), p, g);
}

// same as above, non-const graph overload
template <typename PropertyTag, typename DynamicTag, typename Graph>
typename boost::property_map<Graph, DynamicTag>::type
get_initialized_index_map(CGAL::internal_np::Param_not_found,
                          const PropertyTag p, const DynamicTag dtag, DynamicTag,
                          Graph& g)
{
  // From now on the correct property map has been acquired
  // and there is no need to distinguish between const and non-const mesh
  return get_initialized_dynamic_index_map(get(dtag, g), p, g);
}

template <typename PropertyTag, typename Tag, typename DynamicTag,
          typename Graph,
          typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
class GetInitializedIndexMap
{
public:
  // Check if there is an internal property map; if not, we must a dynamic property map
  typedef typename boost::mpl::if_c<
      CGAL::graph_has_property<Graph, Tag>::value, Tag, DynamicTag>::type    Final_tag;

  typedef typename internal_np::Lookup_named_param_def<
      PropertyTag,
      NamedParameters,
      typename boost::property_map<Graph, Final_tag>::const_type>::type      const_type;

  typedef typename internal_np::Lookup_named_param_def<
      PropertyTag,
      NamedParameters,
      typename boost::property_map<Graph, Final_tag>::type>::type            type;

  static const_type get_const(const PropertyTag p, const Graph& g, const NamedParameters& np)
  {
    return BGL::internal::get_initialized_index_map_const(parameters::get_parameter(np, p),
                                                          p, Final_tag{}, DynamicTag{}, g);
  }

  static type get(const PropertyTag p, Graph& g, const NamedParameters& np)
  {
    return BGL::internal::get_initialized_index_map(parameters::get_parameter(np, p),
                                                    p, Final_tag{}, DynamicTag{}, g);
  }
};

} // namespace internal
} // namespace BGL
} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_INITIALIZED_INTERNAL_INDEX_MAPS_HELPERS
