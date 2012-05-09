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

namespace CGAL {

namespace internal {
  template <class GeomTraits>
  struct First_point_of_segment_3_property_map{
    //classical typedefs
    typedef const typename GeomTraits::Segment_3& key_type;
    typedef typename GeomTraits::Point_3 value_type;
    typedef const typename GeomTraits::Point_3& reference;
    typedef boost::readable_property_map_tag category;
  };
  
  //get function for property map
  template <class GeomTraits>
  inline
  const typename GeomTraits::Point_3&
  get(First_point_of_segment_3_property_map<GeomTraits>,
      const typename GeomTraits::Segment_3& s)
  {
    return s.source();
  }
}//namespace internal


template <class GeomTraits, 
          class Iterator,
          class cache_primitive=Tag_false>
class AABB_segment_primitive : public AABB_primitive< Iterator,
                                                       boost::typed_identity_property_map<typename GeomTraits::Segment_3>,
                                                       internal::First_point_of_segment_3_property_map<GeomTraits>,
                                                       cache_primitive >
{
  typedef AABB_primitive< Iterator,
                          boost::typed_identity_property_map<typename GeomTraits::Segment_3>,
                          internal::First_point_of_segment_3_property_map<GeomTraits>,
                          cache_primitive > Base;
public:
  // constructors
  AABB_segment_primitive(Iterator it) : Base(it){}
};

}  // end namespace CGAL


#endif // CGAL_AABB_SEGMENT_PRIMITIVE_H_

