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

#ifndef CGAL_AABB_FACE_GRAPH_TRIANGLE_PRIMITIVE_H
#define CGAL_AABB_FACE_GRAPH_TRIANGLE_PRIMITIVE_H

#include <CGAL/AABB_primitive.h>
#include <CGAL/internal/AABB_tree/Halfedge_and_face_graph_property_maps.h>
#include <CGAL/Default.h>
#include <boost/mpl/has_xxx.hpp>

namespace CGAL {


#ifndef CGAL_NO_DEPRECATED_CODE
namespace internal_aabb_tree{
  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_facet_const_handle,Facet_const_handle,false)
  template <class FaceGraph,
            bool has_facet_const_handle=Has_facet_const_handle<FaceGraph>::value>
  struct Get_facet_const_handle{
    typedef typename FaceGraph::Facet_const_handle type;
  };
  
  template <class FaceGraph>
  struct Get_facet_const_handle<FaceGraph,false>{
    typedef void* type;
  };
}
#endif

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
 *\tparam VertexPointPMap  is a property map with `boost::graph_traits<HalfedgeGraph>::%vertex_descriptor`
 *   as key type and a \cgal Kernel `Point_3` as value type.
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
           class CacheDatum=Tag_false >
class AABB_face_graph_triangle_primitive
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
  Id type.
  */
  typedef boost::graph_traits<FaceGraph>::face_descriptor Id;
  /// @}

  /*!
  If `OneFaceGraphPerTree` is CGAL::Tag_true, constructs a `Shared_data` object from a reference to the polyhedon `graph`.
  */
  static unspecified_type construct_shared_data( FaceGraph& graph );
  #endif

  // constructors
  /*!
    \tparam Iterator an input iterator with `Id` as value type.
    Constructs a primitive.
  */
  template <class Iterator>
  AABB_face_graph_triangle_primitive(Iterator it, const FaceGraph& graph, VertexPointPMap vppm)
    : Base( Id_(*it),
            Triangle_property_map(const_cast<FaceGraph*>(&graph),vppm),
            Point_property_map(const_cast<FaceGraph*>(&graph),vppm) )
  {}

#ifndef DOXYGEN_RUNNING
  template <class Iterator>
  AABB_face_graph_triangle_primitive(Iterator it, const FaceGraph& graph)
    : Base( Id_(*it),
            Triangle_property_map(const_cast<FaceGraph*>(&graph)),
            Point_property_map(const_cast<FaceGraph*>(&graph)) )
  {}
#ifndef CGAL_NO_DEPRECATED_CODE
  // for backward compatibility with Polyhedron::facets_begin()
  AABB_face_graph_triangle_primitive(typename boost::graph_traits<FaceGraph>::face_descriptor fd, FaceGraph& graph)
    : Base( Id_(fd),
            Triangle_property_map(&graph),
            Point_property_map(&graph) )
  {}

  AABB_face_graph_triangle_primitive(
      typename internal_aabb_tree::Get_facet_const_handle<FaceGraph>::type fd,
      FaceGraph& graph
  ) : Base( Id_(fd.remove_const()),
            Triangle_property_map(&graph),
            Point_property_map(&graph) )
  {}
#endif
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
  construct_shared_data(const FaceGraph& graph, const VertexPointPMap& vpm)
  {
    return Cstr_shared_data::construct_shared_data(const_cast<FaceGraph&>(graph), vpm);
  }

};

}  // end namespace CGAL

#endif // CGAL_AABB_FACE_GRAPH_TRIANGLE_PRIMITIVE_H

