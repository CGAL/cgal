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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot
//


#ifndef CGAL_AABB_SEGMENT_PRIMITIVE_H_
#define CGAL_AABB_SEGMENT_PRIMITIVE_H_

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/AABB_primitive.h>
#include <CGAL/result_of.h>
#include <iterator>

namespace CGAL {

namespace internal {
  template <class GeomTraits, class Iterator>
  struct Source_of_segment_3_iterator_property_map{
    //classical typedefs
    typedef Iterator key_type;
    typedef typename GeomTraits::Point_3 value_type;
    typedef typename cpp11::result_of<
      typename GeomTraits::Construct_source_3(typename GeomTraits::Segment_3)
    >::type reference;
    typedef boost::readable_property_map_tag category;

    inline friend
    typename Source_of_segment_3_iterator_property_map<GeomTraits,Iterator>::reference
    get(Source_of_segment_3_iterator_property_map<GeomTraits,Iterator>, Iterator it)
    {
      return typename GeomTraits::Construct_source_3()( *it );
    }
  };
}//namespace internal


/*!
 * \ingroup PkgAABB_tree
 * Primitive type that uses as identifier an iterator with a 3D segment as `value_type`.
 * The iterator from which the primitive is built should not be invalided
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels `AABBPrimitive`
 *
 * \tparam GeomTraits is a traits class providing the nested type `Point_3` and `Segment_3`.
 *         It also provides the functor `Construct_source_3` that has an operator taking a `Segment_3`
 *         and returning its source as a type convertible to `Point_3`.
 *         In addition `Construct_source_3` must support the result_of protocol.
 * \tparam Iterator is a model of `ForwardIterator` with its value type convertible to `GeomTraits::Segment_3`
 * \tparam CacheDatum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *           the datum is stored in the primitive, while in the latter it is
 *           constructed on the fly to reduce the memory footprint.
 *           The default is `CGAL::Tag_false` (datum is not stored).
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatum>`
 * \sa `AABB_triangle_primitive<Iterator,CacheDatum>`
 * \sa `AABB_halfedge_graph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,CacheDatum>`
 * \sa `AABB_face_graph_triangle_primitive<FaceGraph,OneFaceGraphPerTree,CacheDatum>`
 */
template < class GeomTraits,
           class Iterator,
           class CacheDatum=Tag_false>
class AABB_segment_primitive
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
  ///Constructor from an iterator
  AABB_segment_primitive(Iterator it) : Base(it){}
};

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_SEGMENT_PRIMITIVE_H_

