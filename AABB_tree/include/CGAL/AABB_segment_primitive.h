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


#ifndef CGAL_AABB_SEGMENT_PRIMITIVE_H_
#define CGAL_AABB_SEGMENT_PRIMITIVE_H_

#include <CGAL/AABB_primitive.h>
#include <CGAL/Kernel_traits.h>
#include <iterator>

namespace CGAL {

namespace internal {
  template <class Iterator>
  struct Source_of_segment_3_iterator_property_map{
    //classical typedefs
    typedef typename CGAL::Kernel_traits< typename std::iterator_traits<Iterator>::value_type >::Kernel GeomTraits;
    typedef Iterator key_type;
    typedef typename GeomTraits::Point_3 value_type;
    typedef const typename GeomTraits::Point_3& reference;
    typedef boost::readable_property_map_tag category;
  };

  //get function for property map
  template <class Iterator>
  inline
  typename Source_of_segment_3_iterator_property_map<Iterator>::reference
  get(Source_of_segment_3_iterator_property_map<Iterator>, Iterator it)
  {
    return it->source();
  }
}//namespace internal


/*!
 * \ingroup PkgAABB_tree
 * Primitive type that uses as identifier an iterator with a 3D segment as `value_type`.
 * The iterator from which the primitive is built should not be invalided 
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels `AABBPrimitive`
 *
 * \tparam Iterator is a model of `ForwardIterator`, with `Segment_3<Kernel>`
 *    as value type
 * \tparam cache_datum is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *           the datum is stored in the primitive, while in the latter it is
 *           constructed on the fly to reduce the memory footprint.
 *           The default is `CGAL::Tag_false` (datum is not stored).
 *
 * \sa `AABBPrimitive`
 * \sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,cache_datum>`
 * \sa `AABB_triangle_primitive<Iterator,cache_datum>`
 * \sa `AABB_HalfedgeGraph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,cache_datum>`
 * \sa `AABB_FaceGraph_triangle_primitive<FaceGraph,OneFaceGraphPerTree,cache_datum>`
 */
template < class Iterator,
           class cache_datum=Tag_false>
class AABB_segment_primitive
#ifndef DOXYGEN_RUNNING
  : public AABB_primitive<  Iterator,
                            Input_iterator_property_map<Iterator>,
                            internal::Source_of_segment_3_iterator_property_map<Iterator>,
                            Tag_true,
                            cache_datum >
#endif
{
  typedef AABB_primitive< Iterator,
                          Input_iterator_property_map<Iterator>,
                          internal::Source_of_segment_3_iterator_property_map<Iterator>,
                          Tag_true,
                          cache_datum > Base;
public:
  ///Constructor from an iterator
  AABB_segment_primitive(Iterator it) : Base(it){}

  static typename Base::Shared_data construct_shared_data() {return typename Base::Shared_data();}
};

}  // end namespace CGAL


#endif // CGAL_AABB_SEGMENT_PRIMITIVE_H_

