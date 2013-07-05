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

#ifndef CGAL_AABB_HALFEDGEGRAPH_SEGMENT_PRIMITIVE_H
#define CGAL_AABB_HALFEDGEGRAPH_SEGMENT_PRIMITIVE_H

#include <CGAL/AABB_primitive.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3_property_map.h>

#include <iterator>
#include <boost/mpl/and.hpp>
#include <CGAL/is_iterator.h>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

namespace CGAL {


/*!
 * \ingroup PkgAABB_tree
 * Primitive type for a edge of a polyhedral surface.
 * It wraps an `edge_descriptor` into a 3D segment.
 * The class model of `HalfedgeGraph` from which the primitive is built should not be deleted
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels `AABBPrimitive` if `OneHalfedgeGraphPerTree` is `CGAL::Tag_false`,
 *    and `AABBPrimitiveWithSharedData` if `OneHalfedgeGraphPerTree` is `CGAL::Tag_true`.
 *
 * \tparam HalfedgeGraph is a model of the halfedge graph concept.
 * \tparam OneHalfedgeGraphPerTree is either `CGAL::Tag_true or `CGAL::Tag_false`.
 * In the former case, we guarantee that all the primitives will be from a
 * common `HalfedgeGraph` and some data will be factorized so that the size of
 * the primitive is reduced. In the latter case, the primitives can be from
 * different graphs and extra storage is required in the primitives. The default is `CGAL::Tag_true`.
 * \tparam cache_datum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case, the datum is
 * stored in the primitive, while in the latter it is constructed on the fly to reduce
 * the memory footprint. The default is `CGAL::Tag_false` (datum is not stored).
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,cache_datum>`
 * \sa `AABB_FaceGraph_triangle_primitive<FaceGraph,OneFaceGraphPerTree,cache_datum>`
 */
template < class HalfedgeGraph,
           class OneHalfedgeGraphPerTree=Tag_true,
           class cache_datum=Tag_false >
class AABB_HalfedgeGraph_segment_primitive
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  typename boost::graph_traits<HalfedgeGraph>::edge_descriptor,
                            Segment_from_edge_descriptor_property_map<HalfedgeGraph>,
                            Source_point_from_edge_descriptor<HalfedgeGraph>,
                            OneHalfedgeGraphPerTree,
                            cache_datum >
#endif
{
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor Id_;
  typedef Segment_from_edge_descriptor_property_map<HalfedgeGraph>  Segment_property_map;
  typedef Source_point_from_edge_descriptor<HalfedgeGraph> Point_property_map;

  typedef AABB_primitive< Id_,
                          Segment_property_map,
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
  typedef boost::property_traits< boost::property_map< HalfedgeGraph, vertex_point_t>::type >::value_type Point;
  /*!
  Geometric data type.
  */
  typedef Kernel_traits<Point>::Kernel::Segment_3 Datum;
  /*!
  Id type.
  */
  typedef boost::graph_traits<HalfedgeGraph>::edge_descriptor Id;
  /// @}
#endif

  /*!
  Constructs a primitive.
  \tparam Iterator is an input iterator with `Id` as value type.
  This \ref AABB_tree/AABB_HalfedgeGraph_edge_example.cpp "example" gives a way to call this constructor
  using the the insert-by-range method of the class `AABB_tree<Traits>`.
  */
  template <class Iterator>
  AABB_HalfedgeGraph_segment_primitive(Iterator it, HalfedgeGraph& graph)
    : Base( Id_(*it),
            Segment_property_map(&graph),
            Point_property_map(&graph) ){}

  /// For backward-compatibility with AABB_polyhedron_segment_primitive only
  AABB_HalfedgeGraph_segment_primitive(Id_ id)
    : Base( id,
            Segment_property_map(NULL),
            Point_property_map(NULL) ){}

  static typename Base::Shared_data construct_shared_data( HalfedgeGraph& graph )
  {
    return Base::construct_shared_data(Segment_property_map(&graph), Point_property_map(&graph));
  }

  ///For backward-compatibility with AABB_polyhedron_segment_primitive only
  static typename Base::Shared_data construct_shared_data()
  {
    return Base::construct_shared_data(Segment_property_map(NULL), Point_property_map(NULL));
  }
};

}  // end namespace CGAL


#endif // CGAL_AABB_HALFEDGEGRAPH_SEGMENT_PRIMITIVE_H

