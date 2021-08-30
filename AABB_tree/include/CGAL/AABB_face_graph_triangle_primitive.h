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

#ifndef CGAL_AABB_FACE_GRAPH_TRIANGLE_PRIMITIVE_H
#define CGAL_AABB_FACE_GRAPH_TRIANGLE_PRIMITIVE_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/AABB_primitive.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/Default.h>
#include <boost/mpl/if.hpp>

namespace CGAL {

/*!
 * \ingroup PkgAABBTreeRef
 * Primitive type for a facet of a polyhedral surface.
 * It wraps a handle to a facet of a polyhedron to a 3D triangle.
 * The polyhedron from which the primitive is built should not be deleted
 * while the AABB tree holding the primitive is in use.
 * The triangle type of the primitive (`Datum`) is `CGAL::Kernel_traits< boost::property_traits< VertexPointPMap >::%value_type >::%Kernel::Triangle_3`.
 *
 * \cgalModels `AABBPrimitiveWithSharedData`
 *
 *\tparam FaceGraph is a model of the face graph concept.
 *\tparam VertexPointPMap  is a property map with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
 *   as key type and a \cgal Kernel `Point_3` as value type.
 *                         The default is `typename boost::property_map< FaceGraph,vertex_point_t>::%const_type`.
 *\tparam OneFaceGraphPerTree is either `CGAL::Tag_true` or `CGAL::Tag_false`.
 * In the former case, we guarantee that all the primitives will be from a
 * common `FaceGraph` and some data will be factorized so that the size of
 * the primitive is reduced. In the latter case, the primitives can be from
 * different graphs and extra storage is required in the primitives. The default is `CGAL::Tag_true`.
 *\tparam CacheDatum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case, the datum is stored
 *        in the primitive, while in the latter it is constructed on the fly to reduce the memory footprint.
 *        The default is `CGAL::Tag_false` (datum is not stored).
 *\sa `AABBPrimitive`
 *\sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatum>`
 *\sa `AABB_halfedge_graph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,CacheDatum>`
 */
template < class FaceGraph,
           class VertexPointPMap = Default,
           class OneFaceGraphPerTree = Tag_true,
           class CacheDatum=Tag_false >
class AABB_face_graph_triangle_primitive
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<typename boost::mpl::if_<OneFaceGraphPerTree,
                                                   typename boost::graph_traits<FaceGraph>::face_descriptor,
                                                   std::pair<typename boost::graph_traits<FaceGraph>::face_descriptor, const FaceGraph*> >::type,
                        Triangle_from_face_descriptor_map<
                          FaceGraph,
                          typename Default::Get<VertexPointPMap,
                                                typename boost::property_map< FaceGraph,
                                                                              vertex_point_t>::const_type >::type>,
                        One_point_from_face_descriptor_map<
                          FaceGraph,
                          typename Default::Get<VertexPointPMap,
                                                typename boost::property_map< FaceGraph,
                                                                              vertex_point_t>::const_type >::type>,
                        OneFaceGraphPerTree,
                        CacheDatum >
#endif
{
  typedef typename Default::Get<VertexPointPMap, typename boost::property_map< FaceGraph, vertex_point_t>::const_type >::type VertexPointPMap_;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor FD;
  typedef typename boost::mpl::if_<OneFaceGraphPerTree, FD, std::pair<FD, const FaceGraph*> >::type Id_;

  typedef Triangle_from_face_descriptor_map<FaceGraph,VertexPointPMap_>  Triangle_property_map;
  typedef One_point_from_face_descriptor_map<FaceGraph,VertexPointPMap_> Point_property_map;

  typedef AABB_primitive< Id_,
                          Triangle_property_map,
                          Point_property_map,
                          OneFaceGraphPerTree,
                          CacheDatum > Base;

  FD make_id(FD fd, const FaceGraph&, Tag_true)
  {
    return fd;
  }

  std::pair<FD, const FaceGraph*> make_id(FD fd, const FaceGraph& fg, Tag_false)
  {
    return std::make_pair(fd, &fg);
  }

public:
#ifdef DOXYGEN_RUNNING
  /// \name Types
  /// @{
  /*!
  The point type.
  */
  typedef boost::property_traits<VertexPointPMap>::value_type Point;
  /*!
  Geometric data type.
  */
  typedef Kernel_traits<Point>::Kernel::Triangle_3 Datum;
  /*!
  Id type:
  - `boost::graph_traits<FaceGraph>::%face_descriptor` if `OneFaceGraphPerTree` is `CGAL::Tag_true`
  - `std::pair<boost::graph_traits<FaceGraph>::%face_descriptor, const FaceGraph*>` if `OneFaceGraphPerTree` is `CGAL::Tag_false`
  */
  unspecified_type Id;

  /// @}

  /*!
  If `OneFaceGraphPerTree` is CGAL::Tag_true, constructs a `Shared_data` object from a reference to the polyhedon `graph`.
  */
  static unspecified_type construct_shared_data( FaceGraph& graph );
#else
  typedef typename Base::Id Id;
#endif
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

  // constructors

#ifdef DOXYGEN_RUNNING
  /*!
    constructs a primitive.

    \tparam Iterator an input iterator with `Id` as value type.

    If `VertexPointPMap` is the default of the class, an additional constructor
    is available with `vppm` set to `get(vertex_point, graph)`.
  */
  template <class Iterator>
  AABB_face_graph_triangle_primitive(Iterator it, const FaceGraph& graph, VertexPointPMap vppm);

  /*!
    constructs a primitive.
    If `VertexPointPMap` is the default of the class, an additional constructor
    is available with `vppm` set to `get(vertex_point, graph)`.
  */
  AABB_face_graph_triangle_primitive(face_descriptor fd, const FaceGraph& graph, VertexPointPMap vppm);
#else
  template <class Iterator>
  AABB_face_graph_triangle_primitive(Iterator it, const FaceGraph& graph, VertexPointPMap_ vppm)
    : Base( Id_(make_id(*it, graph, OneFaceGraphPerTree())),
            Triangle_property_map(const_cast<FaceGraph*>(&graph),vppm),
            Point_property_map(const_cast<FaceGraph*>(&graph),vppm) )
  {}

  AABB_face_graph_triangle_primitive(face_descriptor fd, const FaceGraph& graph, VertexPointPMap_ vppm)
    : Base( Id_(make_id(fd, graph, OneFaceGraphPerTree())),
            Triangle_property_map(const_cast<FaceGraph*>(&graph),vppm),
            Point_property_map(const_cast<FaceGraph*>(&graph),vppm) )
  {}

  template <class Iterator>
  AABB_face_graph_triangle_primitive(Iterator it, const FaceGraph& graph)
    : Base( Id_(make_id(*it, graph, OneFaceGraphPerTree())),
            Triangle_property_map(const_cast<FaceGraph*>(&graph)),
            Point_property_map(const_cast<FaceGraph*>(&graph)) )
  {}

  AABB_face_graph_triangle_primitive(face_descriptor fd, const FaceGraph& graph)
    : Base( Id_(make_id(fd, graph, OneFaceGraphPerTree())),
            Triangle_property_map(const_cast<FaceGraph*>(&graph)),
            Point_property_map(const_cast<FaceGraph*>(&graph)) )
  {}
#endif

  /// \internal
  typedef internal::Cstr_shared_data<FaceGraph, Base, Triangle_property_map, Point_property_map, OneFaceGraphPerTree> Cstr_shared_data;
  /// \internal
  static
  typename Cstr_shared_data::Shared_data
  construct_shared_data(const FaceGraph& graph)
  {
    return Cstr_shared_data::construct_shared_data(const_cast<FaceGraph&>(graph));
  }

  static
  typename Cstr_shared_data::Shared_data
  construct_shared_data(const FaceGraph& graph, const VertexPointPMap_& vpm)
  {
    return Cstr_shared_data::construct_shared_data(const_cast<FaceGraph&>(graph), vpm);
  }

};

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_FACE_GRAPH_TRIANGLE_PRIMITIVE_H

