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
  template <class GeomTraits, class Iterator, class PointMap>
  struct Segment_2_from_point_iterator_property_map {
    //classical typedefs
    typedef Iterator key_type;
    typedef typename GeomTraits::Segment_2 value_type;
    typedef typename GeomTraits::Segment_2 reference; // The segments are created on the fly, so working with references is not possible.
    typedef boost::readable_property_map_tag category;
    typedef Segment_2_from_point_iterator_property_map<GeomTraits, Iterator, PointMap> Self;

    Segment_2_from_point_iterator_property_map(Iterator b, Iterator e, PointMap& pmap) : begin(b), end(e), pmap(pmap) {}
    Segment_2_from_point_iterator_property_map() {}

    inline friend reference // Cannot return reference as the Segment does not exist, only the points exist.
    get(Self s, key_type it)
    {
      Iterator it2 = std::next(it);
      if (it2 == s.end)
        return typename GeomTraits::Construct_segment_2()(get(s.pmap, *it), get(s.pmap, *s.begin));
      else
        return typename GeomTraits::Construct_segment_2()(get(s.pmap, *it), get(s.pmap, *it2));
    }

    Iterator begin, end;
    PointMap pmap;
  };

  template <class GeomTraits, class Iterator, class PointMap>
  struct Point_from_iterator_property_map {
    //classical typedefs
    typedef Iterator key_type;
    typedef typename PointMap::value_type value_type;
    typedef const value_type reference;

    typedef boost::readable_property_map_tag category;
    typedef Point_from_iterator_property_map<GeomTraits, Iterator, PointMap> Self;

    Point_from_iterator_property_map() {}
    Point_from_iterator_property_map(PointMap& pmap) :  pmap(pmap) {}

    inline friend reference
      get(Self s, key_type it)
    {
      return get(s.pmap, *it);
    }

    PointMap pmap;
  };
}//namespace internal


/*!
 * \ingroup PkgAABBTreeRef
 * Primitive type that uses as identifier an iterator with a 2D point as `value_type`.
 * The iterator from which the primitive is built should not be invalided
 * while the AABB tree holding the primitive is in use. The `Segment_2` is constructed on the fly using the `Point_2`
 * the identifier is pointing to as source and the next `Point_2` as target.
 *
 * \cgalModels{AABBPrimitive}
 *
 * \tparam GeomTraits is a traits class providing the nested type `Point_2` and `Segment_2`.
 *         It also provides the functor `Construct_segment_2` that has an operator taking two `Point_2`
 *         and returning a `Segment_2`.
 * \tparam Iterator is a model of `ForwardIterator` with its value type convertible to `GeomTraits::Point_2`
 * \tparam PointRange is a model of `ConstRange`. Its value type needs to be compatible to PointMap or `Point_2` in the default case.
 * \tparam CacheDatum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *         the datum is stored in the primitive, while in the latter it is
 *         constructed on the fly to reduce the memory footprint.
 *         The default is `CGAL::Tag_false` (datum is not stored).
 * \tparam PointMap is a model of `ReadablePropertyMap` with its key type being the value type of `PointRange` and the value type being a `Point_2`.
 *         The default is \link Identity_property_map `CGAL::Identity_property_map`\endlink<PointRange::value_type>.
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatum>`
 * \sa `AABB_segment_primitive_2<Iterator,CacheDatum>`
 * \sa `AABB_segment_primitive_3<Iterator,CacheDatum>`
 * \sa `AABB_triangle_primitive_2<Iterator,CacheDatum>`
 * \sa `AABB_triangle_primitive_3<Iterator,CacheDatum>`
 * \sa `AABB_halfedge_graph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,CacheDatum>`
 * \sa `AABB_face_graph_triangle_primitive<FaceGraph,OneFaceGraphPerTree,CacheDatum>`
 */
template < class GeomTraits,
           class Iterator,
           class PointRange,
           class CacheDatum = Tag_false,
           class PointMap = Identity_property_map<typename PointRange::value_type>>
class AABB_polyline_segment_primitive_2
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  Iterator,
                            internal::Segment_2_from_point_iterator_property_map<GeomTraits, Iterator, PointMap>,
                            internal::Point_from_iterator_property_map<GeomTraits, Iterator, PointMap>,
                            Tag_true,
                            CacheDatum >
#endif
{
  typedef AABB_primitive< Iterator,
                          internal::Segment_2_from_point_iterator_property_map<GeomTraits, Iterator, PointMap>,
                          internal::Point_from_iterator_property_map<GeomTraits, Iterator, PointMap>,
                          Tag_true,
                          CacheDatum > Base;
public:
  AABB_polyline_segment_primitive_2(Iterator it, PointRange& poly) : Base(it) {}

  /// \internal
  static typename Base::Shared_data construct_shared_data(PointRange& range, PointMap pmap = PointMap()) {
    return std::make_pair(internal::Segment_2_from_point_iterator_property_map<GeomTraits, Iterator, PointMap>(range.begin(), range.end(), pmap), internal::Point_from_iterator_property_map<GeomTraits, Iterator, PointMap>(pmap));
  }
};

}  // end namespace CGAL


#endif // CGAL_AABB_POLYLINE_SEGMENT_PRIMITIVE_2_H_