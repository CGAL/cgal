// Copyright (c) 2024 GeometryFactory (France).
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


#ifndef CGAL_AABB_INDEXED_TRIANGLE_PRIMITIVE_2_H_
#define CGAL_AABB_INDEXED_TRIANGLE_PRIMITIVE_2_H_

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_primitive.h>
#include <iterator>

namespace CGAL {

namespace internal {

template <class GeomTraits, class Iterator, class PointIterator, class PointMap>
struct Triangle_2_from_index_range_iterator_property_map {
  //classical typedefs
  typedef Iterator key_type;
  typedef typename GeomTraits::Triangle_2 value_type;
  typedef typename GeomTraits::Triangle_2 reference;

  typedef boost::readable_property_map_tag category;
  typedef Triangle_2_from_index_range_iterator_property_map<GeomTraits, Iterator, PointIterator, PointMap> Self;

  Triangle_2_from_index_range_iterator_property_map() {}
  Triangle_2_from_index_range_iterator_property_map(PointIterator b, PointMap& pmap) : begin(b), pmap(pmap) {}

  inline friend value_type
  get(Self s, key_type it)
  {
    return typename GeomTraits::Construct_triangle_2()(get(s.pmap, s.begin[(*it)[0]]), get(s.pmap, s.begin[(*it)[1]]), get(s.pmap, s.begin[(*it)[2]]));
  }

  PointIterator begin;
  PointMap pmap;
};

template <class GeomTraits, class Iterator, class PointIterator, class PointMap>
struct Point_from_indexed_triangle_2_iterator_property_map {
  //classical typedefs
  typedef Iterator key_type;
  typedef typename PointMap::value_type value_type;
  typedef const value_type reference;

  typedef boost::readable_property_map_tag category;
  typedef Point_from_indexed_triangle_2_iterator_property_map<GeomTraits, Iterator, PointIterator, PointMap> Self;

  Point_from_indexed_triangle_2_iterator_property_map() {}
  Point_from_indexed_triangle_2_iterator_property_map(PointIterator b, PointMap &pmap) : begin(b), pmap(pmap) {}

  inline friend reference
  get(Self s, key_type it)
  {
    return get(s.pmap, s.begin[((*it)[0])]);
  }

  PointIterator begin;
  PointMap pmap;
};
}//namespace internal


/*!
 * \ingroup PkgAABBTreeRef
 * Primitive type that uses as identifier an iterator with a range of three indices as `value_type`.
 * The iterator from which the primitive is built should not be invalided
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels{AABBPrimitive}
 *
 * \tparam GeomTraits is a traits class providing the nested type `Point_2` and `Triangle_2`.
 *         It also provides the functor `Construct_triangle_2` that has an operator taking three `Point_2` as
 *         parameters and returns a `Triangle_2`
 * \tparam IndexIterator is a model of `ForwardIterator` with its value type being a `RandomAccessRange` of size 3 with an index type as `value_type`, e.g., `uint8_t`, `uint16_t` or int.
 * \tparam PointRange is a model of `RandomAccessRange`. Its value type needs to be compatible to `PointMap` or `Point_2` in the default case.
 * \tparam CacheDatum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *         the datum is stored in the primitive, while in the latter it is
 *         constructed on the fly to reduce the memory footprint.
 *         The default is `CGAL::Tag_false` (datum is not stored).
 * \tparam PointMap is a model of `ReadablePropertyMap` with its key type being the value type of `PointRange` and the value type being a `Point_2`.
 *         The default is \link Identity_property_map `CGAL::Identity_property_map`\endlink<PointRange::value_type>.
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatum>`
 * \sa `AABB_segment_primitive_2<GeomTraits,Iterator,CacheDatum>`
 * \sa `AABB_triangle_primitive_2<GeomTraits,Iterator,CacheDatum>`
 * \sa `AABB_triangle_primitive_3<GeomTraits,Iterator,CacheDatum>`
 */
template < class GeomTraits,
           class IndexIterator,
           class PointRange,
           class CacheDatum = Tag_false,
           class PointMap = Identity_property_map<typename PointRange::value_type>>
class AABB_indexed_triangle_primitive_2
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  IndexIterator,
                            internal::Triangle_2_from_index_range_iterator_property_map<GeomTraits, IndexIterator, typename PointRange::iterator, PointMap>,
                            internal::Point_from_indexed_triangle_2_iterator_property_map<GeomTraits, IndexIterator, typename PointRange::iterator, PointMap>,
                            Tag_true,
                            CacheDatum >
#endif
{
  typedef AABB_primitive< IndexIterator,
                          internal::Triangle_2_from_index_range_iterator_property_map<GeomTraits, IndexIterator, typename PointRange::iterator, PointMap>,
                          internal::Point_from_indexed_triangle_2_iterator_property_map<GeomTraits, IndexIterator, typename PointRange::iterator, PointMap>,
                          Tag_true,
                          CacheDatum > Base;
public:
  ///constructor from an iterator
  AABB_indexed_triangle_primitive_2(IndexIterator it, PointRange&) : Base(it) {}

  /// \internal
  static typename Base::Shared_data construct_shared_data(PointRange &range, PointMap pmap = PointMap()) {
    return std::make_pair(
      internal::Triangle_2_from_index_range_iterator_property_map<GeomTraits, IndexIterator, typename PointRange::iterator, PointMap>(range.begin(), pmap),
      internal::Point_from_indexed_triangle_2_iterator_property_map<GeomTraits, IndexIterator, typename PointRange::iterator, PointMap>(range.begin(), pmap));
  }
};

}  // end namespace CGAL

#endif // CGAL_AABB_INDEXED_TRIANGLE_PRIMITIVE_2_H_
