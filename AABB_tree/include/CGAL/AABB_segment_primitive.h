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
//******************************************************************************
// File Description :
//
//******************************************************************************

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


template < class Iterator,
           class cache_datum=Tag_false>
class AABB_segment_primitive : public AABB_primitive< Iterator,
                                                      Input_iterator_property_map<Iterator>,
                                                      internal::Source_of_segment_3_iterator_property_map<Iterator>,
                                                      cache_datum >
{
  typedef AABB_primitive< Iterator,
                          Input_iterator_property_map<Iterator>,
                          internal::Source_of_segment_3_iterator_property_map<Iterator>,
                          cache_datum > Base;
public:
  // constructors
  AABB_segment_primitive(Iterator it) : Base(it){}
};

}  // end namespace CGAL


#endif // CGAL_AABB_SEGMENT_PRIMITIVE_H_

