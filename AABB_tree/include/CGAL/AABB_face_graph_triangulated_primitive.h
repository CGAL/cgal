// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_AABB_face_graph_triangulated_primitive_H
#define CGAL_AABB_face_graph_triangulated_primitive_H

#include <CGAL/AABB_primitive.h>
#include <CGAL/internal/AABB_tree/Halfedge_and_face_graph_property_maps.h>
#include <CGAL/Default.h>

namespace CGAL {

/*!
 * \ingroup PkgAABB_tree
 * Primitive type for a facet of a polyhedral surface.
 * It wraps a handle to a facet of a polyhedron to a 3D triangle.
 * The polyhedron from which the primitive is built should not be deleted
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels `AABBPrimitiveWithSharedData`
 *
 *\tparam FaceGraph is a model of the face graph concept.
 *\tparam VertexPointPMap  is a property map with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
 *   as key type and a \cgal Kernel `Point_3` as value type.
 *                         The default is `typename boost::property_map< FaceGraph,vertex_point_t>::%type`.
 *\tparam OneFaceGraphPerTree is either `CGAL::Tag_true` or `CGAL::Tag_false`.
 * In the former case, we guarantee that all the primitives will be from a
 * common polyhedron and some data will be factorized so that the size of
 * the primitive is reduced. In the latter case, the primitives can be from
 * different polyhedra and extra storage is required in the primitives. The default is `CGAL::Tag_true`.
 *        This parameter is useless for the moment and will be useful in an upcoming release of \cgal.
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
           class CacheDatum=Tag_true>
class AABB_face_graph_triangulated_primitive
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<typename boost::graph_traits<FaceGraph>::face_descriptor,
                        Triangle_from_face_descriptor_property_map<
                          FaceGraph,
                          typename Default::Get<VertexPointPMap,
                                                typename boost::property_map< FaceGraph,
                                                                              vertex_point_t>::type >::type>,
                        One_point_from_face_descriptor_property_map<
                          FaceGraph,
                          typename Default::Get<VertexPointPMap,
                                                typename boost::property_map< FaceGraph,
                                                                              vertex_point_t>::type >::type>,
                        OneFaceGraphPerTree,
                        CacheDatum >
#endif
{
  typedef typename Default::Get<VertexPointPMap, typename boost::property_map< FaceGraph, vertex_point_t>::type >::type VertexPointPMap_;

  typedef typename boost::graph_traits<FaceGraph>::face_descriptor Id_;
  typedef Triangle_from_face_descriptor_property_map<FaceGraph,VertexPointPMap_>  Triangle_property_map;
  typedef One_point_from_face_descriptor_property_map<FaceGraph,VertexPointPMap_> Point_property_map;

  typedef AABB_primitive< Id_,
                          Triangle_property_map,
                          Point_property_map,
                          OneFaceGraphPerTree,
                          CacheDatum > Base;

public:

 typedef typename Base::Id Id;
 typedef typename Base::Datum Datum;

#ifndef DOXYGEN_RUNNING
  AABB_face_graph_triangulated_primitive(Datum tr, Id id, const FaceGraph& graph)
    : Base( tr,
            Id_(id),
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

#endif // CGAL_AABB_face_graph_triangulated_primitive_H

 
