// Copyright (c) 2024 GeometryFactory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau
//


#ifndef CGAL_AABB_POLYLINE_SEGMENT_PRIMITIVE_2_H_
#define CGAL_AABB_POLYLINE_SEGMENT_PRIMITIVE_2_H_

#include <CGAL/license/AABB_tree.h>


#include <CGAL/AABB_primitive.h>
#include <iterator>

namespace CGAL {

namespace internal {
  template <class GeomTraits, class Iterator>
  struct Segment_2_from_point_iterator_property_map {
    //classical typedefs
    typedef Iterator key_type;
    typedef typename GeomTraits::Segment_2 value_type;
    typedef typename GeomTraits::Segment_2 reference; // The segments are created on the fly, so working with references is not possible.
    typedef boost::readable_property_map_tag category;
    typedef Segment_2_from_point_iterator_property_map<GeomTraits, Iterator> Self;

    Segment_2_from_point_iterator_property_map(Iterator b, Iterator e) : begin(b), end(e) {}
    Segment_2_from_point_iterator_property_map() {}

    inline friend reference // Cannot return reference as the Segment does not exist, only the points exist.
    get(Self s, key_type it)
    {
      typename Iterator::value_type& p = *it;
      it++;
      if (it == s.end)
        return typename GeomTraits::Construct_segment_2()(p, *s.begin);
      else
        return typename GeomTraits::Construct_segment_2()( p, *it );
    }

    Iterator begin, end;
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
 * \tparam GeomTraits is a traits class providing the nested type `Point_2` and `Segment_2`.
 *         It also provides the functor `Construct_source_2` that has an operator taking a `Segment_2`
 *         and returning its source as a type convertible to `Point_2`.
 * \tparam Iterator is a model of `ForwardIterator` with its value type convertible to `GeomTraits::Segment_2`
 * \tparam CacheDatum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *           the datum is stored in the primitive, while in the latter it is
 *           constructed on the fly to reduce the memory footprint.
 *           The default is `CGAL::Tag_false` (datum is not stored).
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatum>`
 * \sa `AABB_triangle_primitive_3<Iterator,CacheDatum>`
 * \sa `AABB_halfedge_graph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,CacheDatum>`
 * \sa `AABB_face_graph_triangle_primitive<FaceGraph,OneFaceGraphPerTree,CacheDatum>`
 */
template < class PointRange,
           class GeomTraits,
           class Iterator,
           class CacheDatum=Tag_false>
class AABB_polyline_segment_primitive_2
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  Iterator,
                            internal::Segment_2_from_point_iterator_property_map<GeomTraits, Iterator>,
                            Input_iterator_property_map<Iterator>,
                            Tag_true,
                            CacheDatum >
#endif
{
  typedef AABB_primitive< Iterator,
                          internal::Segment_2_from_point_iterator_property_map<GeomTraits, Iterator>,
                          Input_iterator_property_map<Iterator>,
                          Tag_true,
                          CacheDatum > Base;
public:
  AABB_polyline_segment_primitive_2(Iterator it, PointRange& poly) : Base(it) {}

  /// \internal
  static typename Base::Shared_data construct_shared_data(PointRange& range) {
    return std::make_pair(internal::Segment_2_from_point_iterator_property_map<GeomTraits, Iterator>(range.begin(), range.end()), Input_iterator_property_map<Iterator>());
  }
};

}  // end namespace CGAL


#endif // CGAL_AABB_POLYLINE_SEGMENT_PRIMITIVE_2_H_
