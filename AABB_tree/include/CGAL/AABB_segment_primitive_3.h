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


#ifndef CGAL_AABB_SEGMENT_PRIMITIVE_3_H_
#define CGAL_AABB_SEGMENT_PRIMITIVE_3_H_

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/AABB_primitive.h>
#include <iterator>

namespace CGAL {

namespace internal {
  template <class GeomTraits, class Iterator>
  struct Source_of_segment_3_iterator_property_map{
    //classical typedefs
    typedef Iterator key_type;
    typedef typename GeomTraits::Point_3 value_type;
    // typedef decltype(
    //   std::declval<typename GeomTraits::Construct_source_3>()(
    //     std::declval<typename GeomTraits::Segment_3>())) reference;
    typedef decltype(
      typename GeomTraits::Construct_source_3()(
        *std::declval<key_type&>())) reference;
    typedef boost::readable_property_map_tag category;
    typedef Source_of_segment_3_iterator_property_map<GeomTraits, Iterator> Self;

    inline friend reference
    get(Self, key_type it)
    {
      return typename GeomTraits::Construct_source_3()( *it );
    }
  };
}//namespace internal


/*!
 * \ingroup PkgAABBTreeRef
 * Primitive type that uses as identifier an iterator with a 3D segment as `value_type`.
 * The iterator from which the primitive is built should not be invalided
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels{AABBPrimitive}
 *
 * \tparam GeomTraits is a traits class providing the nested type `Point_3` and `Segment_3`.
 *         It also provides the functor `Construct_source_3` that has an operator taking a `Segment_3`
 *         and returning its source as a type convertible to `Point_3`.
 * \tparam Iterator is a model of `ForwardIterator` with its value type convertible to `GeomTraits::Segment_3`
 * \tparam CacheDatum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *           the datum is stored in the primitive, while in the latter it is
 *           constructed on the fly to reduce the memory footprint.
 *           The default is `CGAL::Tag_false` (datum is not stored).
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatum>`
 * \sa `AABB_segment_primitive_2<GeomTraits,Iterator,CacheDatum>`
 * \sa `AABB_triangle_primitive_3<GeomTraits,Iterator,CacheDatum>`
 * \sa `AABB_halfedge_graph_segment_primitive<HalfedgeGraph,VertexPointPMap,OneHalfedgeGraphPerTree,CacheDatum>`
 */
template < class GeomTraits,
           class Iterator,
           class CacheDatum=Tag_false>
class AABB_segment_primitive_3
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  Iterator,
                            Input_iterator_property_map<Iterator>,
                            internal::Source_of_segment_3_iterator_property_map<GeomTraits, Iterator>,
                            Tag_false,
                            CacheDatum >
#endif
{
  typedef AABB_primitive< Iterator,
                          Input_iterator_property_map<Iterator>,
                          internal::Source_of_segment_3_iterator_property_map<GeomTraits, Iterator>,
                          Tag_false,
                          CacheDatum > Base;
public:
  ///constructor from an iterator
  AABB_segment_primitive_3(Iterator it) : Base(it){}
};

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_SEGMENT_PRIMITIVE_3_H_
