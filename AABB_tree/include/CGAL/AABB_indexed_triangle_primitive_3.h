// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Leo Valque
//


#ifndef CGAL_AABB_INDEXED_TRIANGLE_PRIMITIVE_3_H_
#define CGAL_AABB_INDEXED_TRIANGLE_PRIMITIVE_3_H_

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_primitive.h>
#include <iterator>

namespace CGAL {

namespace internal {

template <class GeomTraits, class PointRange, class FaceRange, class PointMap>
struct Triangle_3_from_triangle_soup_property_map
{
  using key_type = std::size_t;
  using value_type = typename GeomTraits::Triangle_3;
  using reference = value_type;

  using category = boost::readable_property_map_tag;
  using Self = Triangle_3_from_triangle_soup_property_map<GeomTraits, PointRange, FaceRange, PointMap>;

  Triangle_3_from_triangle_soup_property_map(){}
  template <typename std::enable_if<std::is_default_constructible<PointMap>::value, int>::type = 0>
  Triangle_3_from_triangle_soup_property_map(const PointRange &pts_, const FaceRange &triangles_) : pts(&pts_), triangles(&triangles_){}
  Triangle_3_from_triangle_soup_property_map(const PointRange &pts_, const FaceRange &triangles_, PointMap pmap) : pts(&pts_), triangles(&triangles_), pmap(pmap) {}

  inline friend value_type
  get(const Self &s, const key_type &i)
  {
    return GeomTraits().construct_triangle_3_object()(get(s.pmap, (*s.pts)[(*s.triangles)[i][0]]),
                                                      get(s.pmap, (*s.pts)[(*s.triangles)[i][1]]),
                                                      get(s.pmap, (*s.pts)[(*s.triangles)[i][2]]));
  }

  const PointRange *pts;
  const FaceRange *triangles;
  PointMap pmap;
};

template <class GeomTraits, class PointRange, class FaceRange, class PointMap>
struct Reference_point_from_triangle_soup_property_map
{
  using key_type = std::size_t;
  using value_type = typename GeomTraits::Point_3;
  using reference = value_type;

  using category = boost::readable_property_map_tag;
  using Self = Reference_point_from_triangle_soup_property_map<GeomTraits, PointRange, FaceRange, PointMap>;

  Reference_point_from_triangle_soup_property_map(){}
  template <typename std::enable_if<std::is_default_constructible<PointMap>::value, int>::type = 0>
  Reference_point_from_triangle_soup_property_map(const PointRange &pts_, const FaceRange &triangles_) : pts(&pts_), triangles(&triangles_){}
  Reference_point_from_triangle_soup_property_map(const PointRange &pts_, const FaceRange &triangles_, PointMap pmap) : pts(&pts_), triangles(&triangles_), pmap(pmap) {}

  inline friend value_type
  get(const Self &s, const key_type &i)
  {
    return get(s.pmap, (*s.pts)[(*s.triangles)[i][0]]);
  }

  const PointRange *pts;
  const FaceRange *triangles;
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
 * \tparam GeomTraits is a traits class providing the nested type `Point_3` and `Triangle_3`.
 *         It also provides the functor `Construct_triangle_3` that has an operator taking three `Point_3` as
 *         parameters and returns a `Triangle_3`
 * \tparam PointRange is a model of `RandomAccessRange`. Its value type needs to be compatible to `PointMap` or `Point_3` in the default case.
 * \tparam FaceRange  is a model of `RandomAccessRange`. Its value type needs to a `RandomAccessRange` of size 3 with value_type begin `std::size_t`.
 * \tparam CacheDatum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *         the datum is stored in the primitive, while in the latter it is
 *         constructed on the fly to reduce the memory footprint.
 *         The default is `CGAL::Tag_false` (datum is not stored).
 * \tparam PointMap is a model of `ReadablePropertyMap` with its key type being the value type of `PointRange` and the value type being a `Point_3`.
 *         The default is \link Identity_property_map `CGAL::Identity_property_map`\endlink<PointRange::value_type>.
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatum>`
 * \sa `AABB_segment_primitive_2<GeomTraits,Iterator,CacheDatum>`
 * \sa `AABB_triangle_primitive_2<GeomTraits,Iterator,CacheDatum>`
 * \sa `AABB_triangle_primitive_3<GeomTraits,Iterator,CacheDatum>`
 */
template < class GeomTraits,
           class PointRange,
           class FaceRange,
           class CacheDatum = Tag_false,
           class PointMap = Identity_property_map<typename PointRange::value_type> >
class AABB_indexed_triangle_primitive_3
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  std::size_t,
                            internal::Triangle_3_from_triangle_soup_property_map<GeomTraits, PointRange, FaceRange, PointMap>,
                            internal::Reference_point_from_triangle_soup_property_map<GeomTraits, PointRange, FaceRange, PointMap>,
                            Tag_true,
                            CacheDatum >
#endif
{
  using Base = AABB_primitive<  std::size_t,
                                internal::Triangle_3_from_triangle_soup_property_map<GeomTraits, PointRange, FaceRange, PointMap>,
                                internal::Reference_point_from_triangle_soup_property_map<GeomTraits, PointRange, FaceRange, PointMap>,
                                Tag_true,
                                CacheDatum >;
  using Triangle_property_map = internal::Triangle_3_from_triangle_soup_property_map<GeomTraits, PointRange, FaceRange, PointMap>;
  using Point_property_map = internal::Reference_point_from_triangle_soup_property_map<GeomTraits, PointRange, FaceRange, PointMap>;
  using Face_iterator = typename FaceRange::iterator;
  using Face_const_iterator = typename FaceRange::const_iterator;
public:
  ///constructor from an iterator
  template<typename std::enable_if<std::is_default_constructible<PointMap>::value, int>::type = 0>
  AABB_indexed_triangle_primitive_3(Face_const_iterator it, const PointRange&, const FaceRange& triangles) : Base(std::distance(triangles.begin(), it)){}
  template<typename std::enable_if<std::is_default_constructible<PointMap>::value, int>::type = 0>
  AABB_indexed_triangle_primitive_3(Face_iterator it, const PointRange&, const FaceRange& triangles) : Base(std::size_t(std::distance(triangles.begin(), Face_const_iterator(it)))){}
  template<class IndexIterator, typename std::enable_if<std::is_default_constructible<PointMap>::value, int>::type = 0>
  AABB_indexed_triangle_primitive_3(IndexIterator it, const PointRange&, const FaceRange&) : Base(it){}
  template <typename std::enable_if<std::is_default_constructible<PointMap>::value, int>::type = 0>
  AABB_indexed_triangle_primitive_3(std::size_t i, const PointRange&, const FaceRange&) : Base(i){}

  AABB_indexed_triangle_primitive_3(Face_const_iterator it, const PointRange&, const FaceRange& triangles, PointMap) : Base(std::distance(triangles.begin(), it)){}
  AABB_indexed_triangle_primitive_3(Face_iterator it, const PointRange&, const FaceRange& triangles, PointMap) : Base(std::size_t(std::distance(triangles.begin(), Face_const_iterator(it)))){}
  template<class IndexIterator>
  AABB_indexed_triangle_primitive_3(IndexIterator it, const PointRange&, const FaceRange&, PointMap) : Base(it){}
  AABB_indexed_triangle_primitive_3(std::size_t i, const PointRange&, const FaceRange&, PointMap) : Base(i){}

  /// \internal
  static typename Base::Shared_data construct_shared_data(const PointRange &pts, const FaceRange &triangles, PointMap pmap) {
    return std::make_pair(
      Triangle_property_map(pts, triangles, pmap),
      Point_property_map(pts, triangles, pmap));
  }
  template <typename std::enable_if<std::is_default_constructible<PointMap>::value, int>::type = 0>
  static typename Base::Shared_data construct_shared_data(const PointRange &pts, const FaceRange &triangles) {
    return std::make_pair(
      Triangle_property_map(pts, triangles),
      Point_property_map(pts, triangles));
  }
};

}  // end namespace CGAL

#endif // CGAL_AABB_INDEXED_TRIANGLE_PRIMITIVE_3_H_
