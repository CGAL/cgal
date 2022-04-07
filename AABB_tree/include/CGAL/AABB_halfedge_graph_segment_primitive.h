// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_AABB_HALFEDGE_GRAPH_SEGMENT_PRIMITIVE_H
#define CGAL_AABB_HALFEDGE_GRAPH_SEGMENT_PRIMITIVE_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/AABB_primitive.h>
#include <CGAL/boost/graph/property_maps.h>

#include <iterator>
#include <boost/mpl/and.hpp>
#include <CGAL/is_iterator.h>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>

#include <CGAL/Default.h>

namespace CGAL {


/*!
 * \ingroup PkgAABBTreeRef
 * Primitive type for a edge of a polyhedral surface.
 * It wraps an `edge_descriptor` into a 3D segment.
 * The class model of `HalfedgeGraph` from which the primitive is built should not be deleted
 * while the AABB tree holding the primitive is in use.
 * The type of the 3D segment is taken from the kernel of the point type which is the value type
 * of `VertexPointPMap` (using the `Kernel_traits` mechanism).
 * The segment type of the primitive (`Datum`) is `CGAL::Kernel_traits< boost::property_traits< VertexPointPMap >::%value_type >::%Kernel::Segment_3`.
 *
 * \cgalModels `AABBPrimitive` if `OneHalfedgeGraphPerTree` is `CGAL::Tag_false`,
 *    and `AABBPrimitiveWithSharedData` if `OneHalfedgeGraphPerTree` is `CGAL::Tag_true`.
 *
 * \tparam HalfedgeGraph is a model of the halfedge graph concept.
 *   as key type and a \cgal Kernel `Point_3` as value type.
 * \tparam VertexPointPMap is a property map with `boost::graph_traits<HalfedgeGraph>::%vertex_descriptor`.
 *                         The default is `typename boost::property_map< HalfedgeGraph,vertex_point_t>::%const_type`.
 * \tparam OneHalfedgeGraphPerTree is either `CGAL::Tag_true` or `CGAL::Tag_false`.
 * In the former case, we guarantee that all the primitives will be from a
 * common `HalfedgeGraph` and some data will be factorized so that the size of
 * the primitive is reduced. In the latter case, the primitives can be from
 * different graphs and extra storage is required in the primitives. The default is `CGAL::Tag_true`.
 * \tparam CacheDatum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case, the datum is
 * stored in the primitive, while in the latter it is constructed on the fly to reduce
 * the memory footprint. The default is `CGAL::Tag_false` (datum is not stored).
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatum>`
 * \sa `AABB_face_graph_triangle_primitive<FaceGraph,OneFaceGraphPerTree,CacheDatum>`
 * \sa \link BGLPolyGT `boost::graph_traits<Polyhedron>` \endlink
 */
template < class HalfedgeGraph,
           class VertexPointPMap = Default,
           class OneHalfedgeGraphPerTree = Tag_true,
           class CacheDatum = Tag_false >
class AABB_halfedge_graph_segment_primitive
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  typename boost::mpl::if_<OneHalfedgeGraphPerTree,
                                                     typename boost::graph_traits<HalfedgeGraph>::edge_descriptor,
                                                     std::pair<typename boost::graph_traits<HalfedgeGraph>::edge_descriptor, const HalfedgeGraph*> >::type,
                            Segment_from_edge_descriptor_map<
                              HalfedgeGraph,
                              typename Default::Get<VertexPointPMap,
                                                    typename boost::property_map< HalfedgeGraph,
                                                                                  vertex_point_t>::const_type >::type >,
                            Source_point_from_edge_descriptor_map<
                              HalfedgeGraph,
                              typename Default::Get<VertexPointPMap,
                                                    typename boost::property_map< HalfedgeGraph,
                                                                                  vertex_point_t>::const_type >::type >,
                            OneHalfedgeGraphPerTree,
                            CacheDatum >
#endif
{
  typedef typename Default::Get<VertexPointPMap,typename boost::property_map< HalfedgeGraph,vertex_point_t>::const_type >::type  VertexPointPMap_;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor ED;
  typedef typename boost::mpl::if_<OneHalfedgeGraphPerTree, ED, std::pair<ED, const HalfedgeGraph*> >::type Id_;

  typedef Segment_from_edge_descriptor_map<HalfedgeGraph,VertexPointPMap_>  Segment_property_map;
  typedef Source_point_from_edge_descriptor_map<HalfedgeGraph,VertexPointPMap_> Point_property_map;

  typedef AABB_primitive< Id_,
                          Segment_property_map,
                          Point_property_map,
                          OneHalfedgeGraphPerTree,
                          CacheDatum > Base;

  ED make_id(ED ed, const HalfedgeGraph&, Tag_true)
  {
    return ed;
  }

  std::pair<ED, const HalfedgeGraph*> make_id(ED ed, const HalfedgeGraph& fg, Tag_false)
  {
    return std::make_pair(ed, &fg);
  }

public:

#ifdef DOXYGEN_RUNNING
 /// \name Types
  /// @{
  /*!
  The point type.
  */
  typedef boost::property_traits< boost::property_map< HalfedgeGraph, vertex_point_t>::const_type >::value_type Point;
  /*!
  Geometric data type.
  */
  typedef Kernel_traits<Point>::Kernel::Segment_3 Datum;
  /*!
  Id type:
  - `boost::graph_traits<HalfedgeGraph>::%edge_descriptor` if `OneHalfedgeGraphPerTree` is `Tag_true`
  - `std::pair<boost::graph_traits<HalfedgeGraph>::%edge_descriptor, const HalfedgeGraph*>` if `OneHalfedgeGraphPerTree` is `Tag_false`
  */
  unspecified_type Id;
  /// @}

  /*!
  If `OneHalfedgeGraphPerTree` is CGAL::Tag_true, constructs a `Shared_data` object from a reference to the halfedge graph.
  */
  static unspecified_type construct_shared_data( HalfedgeGraph& graph );
#else
  typedef typename Base::Id Id;
#endif
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor edge_descriptor;

#ifdef DOXYGEN_RUNNING
  /*!
  constructs a primitive.

  \tparam Iterator is an input iterator with `Id` as value type.

  This \ref AABB_tree/AABB_halfedge_graph_edge_example.cpp "example" gives a way to call this constructor
  using the insert-by-range method of the class `AABB_tree<Traits>`.
  If `VertexPointPMap` is the default of the class, an additional constructor
  is available with `vppm` set to `boost::get(vertex_point, graph)`.
  */
  template <class Iterator>
  AABB_halfedge_graph_segment_primitive(Iterator it, const HalfedgeGraph& graph, VertexPointPMap vppm);

  /*!
  constructs a primitive.
  If `VertexPointPMap` is the default of the class, an additional constructor
  is available with `vppm` set to `boost::get(vertex_point, graph)`.
  */
  AABB_halfedge_graph_segment_primitive(edge_descriptor ed, const HalfedgeGraph& graph, VertexPointPMap vppm);
#else
  template <class Iterator>
  AABB_halfedge_graph_segment_primitive(Iterator it, const HalfedgeGraph& graph, VertexPointPMap_ vppm)
    : Base( Id_(make_id(*it, graph, OneHalfedgeGraphPerTree())),
            Segment_property_map(const_cast<HalfedgeGraph*>(&graph), vppm),
            Point_property_map(const_cast<HalfedgeGraph*>(&graph), vppm) )
  {}

  AABB_halfedge_graph_segment_primitive(edge_descriptor ed, const HalfedgeGraph& graph, VertexPointPMap_ vppm)
    : Base( Id_(make_id(ed, graph, OneHalfedgeGraphPerTree())),
            Segment_property_map(const_cast<HalfedgeGraph*>(&graph), vppm),
            Point_property_map(const_cast<HalfedgeGraph*>(&graph), vppm) )
  {}

  template <class Iterator>
  AABB_halfedge_graph_segment_primitive(Iterator it, const HalfedgeGraph& graph)
    : Base( Id_(make_id(*it, graph, OneHalfedgeGraphPerTree())),
            Segment_property_map(const_cast<HalfedgeGraph*>(&graph)),
            Point_property_map(const_cast<HalfedgeGraph*>(&graph)) ){}

  AABB_halfedge_graph_segment_primitive(edge_descriptor ed, const HalfedgeGraph& graph)
    : Base( Id_(make_id(ed, graph, OneHalfedgeGraphPerTree())),
            Segment_property_map(const_cast<HalfedgeGraph*>(&graph)),
            Point_property_map(const_cast<HalfedgeGraph*>(&graph)) ){}
#endif

  /// \internal
  typedef internal::Cstr_shared_data<HalfedgeGraph, Base, Segment_property_map, Point_property_map, OneHalfedgeGraphPerTree> Cstr_shared_data;
  /// \internal
  static
  typename Cstr_shared_data::Shared_data
  construct_shared_data(const HalfedgeGraph& graph)
  {
    return Cstr_shared_data::construct_shared_data(const_cast<HalfedgeGraph&>(graph));
  }

  static
  typename Cstr_shared_data::Shared_data
  construct_shared_data(const HalfedgeGraph& graph, const VertexPointPMap_& vpm)
  {
    return Cstr_shared_data::construct_shared_data(const_cast<HalfedgeGraph&>(graph), vpm);
  }
};

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_HALFEDGE_GRAPH_SEGMENT_PRIMITIVE_H
