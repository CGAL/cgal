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

#ifndef CGAL_AABB_FACEGRAPH_TRIANGLE_PRIMITIVE_H
#define CGAL_AABB_FACEGRAPH_TRIANGLE_PRIMITIVE_H

#include <CGAL/AABB_primitive.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3_property_map.h>

namespace CGAL {

/*!
 * \ingroup PkgAABB_tree
 * Primitive type for a facet of a polyhedral surface.
 * It wraps a `face_descriptor` of a class model of `FaceGraph` to a 3D triangle.
 * The class model of `FaceGraph` from which the primitive is built should not be deleted
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels `AABBPrimitive` if `OneFaceGraphPerTree` is `CGAL::Tag_false`,
 *    and `AABBPrimitiveWithSharedData` if `OneFaceGraphPerTree` is `CGAL::Tag_true`.
 *
 *\tparam FaceGraph is a model of the face graph concept.
 *\tparam OneFaceGraphPerTree is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *        we guarantee that all the primitives will be from a common `FaceGraph` and some data
 *        will be factorized so that the size of the primitive is reduced. In the latter case,
 *        the primitives can be from different graphs and extra storage is required in the primitives.
 *        The default is `CGAL::Tag_true`.
 *\tparam cache_datum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case, the datum is stored
 *        in the primitive, while in the latter it is constructed on the fly to reduce the memory footprint.
 *        The default is `CGAL::Tag_false` (datum is not stored).
 *\sa `AABBPrimitive`
 *\sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,cache_datum>`
 *\sa `AABB_HalfedgeGraph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,cache_datum>`
 */
template < class FaceGraph,
           class OneFaceGraphPerTree=Tag_true,
           class cache_datum=Tag_false,
           class Id_=typename FaceGraph::Face_handle //this one should be autodetected using face_descriptor
            >
class AABB_FaceGraph_triangle_primitive
#ifndef DOXYGEN_RUNNING
: public AABB_primitive< Id_,
                         Triangle_from_facet_handle_property_map<FaceGraph>,
                         One_point_from_facet_handle_property_map<FaceGraph>,
                         OneFaceGraphPerTree,
                         cache_datum >
#endif
{
  typedef Triangle_from_facet_handle_property_map<FaceGraph>  Triangle_property_map;
  typedef One_point_from_facet_handle_property_map<FaceGraph> Point_property_map;

  typedef AABB_primitive< Id_,
                          Triangle_property_map,
                          Point_property_map,
                          Tag_true,
                          cache_datum > Base;

public:
  #ifdef DOXYGEN_RUNNING
  /// \name Types
  /// @{
  /*!
  The point type.
  */
  typedef boost::property_traits< boost::property_map< FaceGraph, vertex_point_t>::type >::value_type Point;
  /*!
  Geometric data type.
  */
  typedef Kernel_traits<Point>::Kernel::Triangle_3 Datum;
  /*!
  Id type.
  */
  typedef boost::graph_traits<FaceGraph>::face_descriptor Id;
  /// @}
  #endif

  // constructors
  /*!
    \tparam Iterator an input iterator with `Id` as value type.
    Constructs a primitive.
  */
  template <class Iterator>
  AABB_FaceGraph_triangle_primitive(Iterator it, const FaceGraph& graph)
    : Base( Id_(*it),
            Triangle_property_map(&graph),
            Point_property_map(&graph) ){}

  //for backward-compatibility with AABB_polyhedron_triangle_primitive
  AABB_FaceGraph_triangle_primitive(Id_ id)
    : Base( id,
            Triangle_property_map(NULL),
            Point_property_map(NULL) ){}

  static typename Base::Shared_data construct_shared_data( const FaceGraph& graph )
  {
    return Base::construct_shared_data(Triangle_property_map(&graph), Point_property_map(&graph));
  }

  //for backward-compatibility with AABB_polyhedron_triangle_primitive
  static typename Base::Shared_data construct_shared_data()
  {
    return Base::construct_shared_data(Triangle_property_map(NULL), Point_property_map(NULL));
  }

};

}  // end namespace CGAL

#endif // CGAL_AABB_FACEGRAPH_TRIANGLE_PRIMITIVE_H

