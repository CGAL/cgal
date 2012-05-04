// Copyright (c) 2009, 2011 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Pierre Alliez, Stephane Tayeb, Sebastien Loriot
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_SEGMENT_PRIMITIVE_H_
#define CGAL_AABB_SEGMENT_PRIMITIVE_H_

#include <CGAL/internal/AABB_tree/Primitive_caching.h>
#include <CGAL/property_map.h>

namespace CGAL {

template <class GeomTraits, 
          class Iterator,
          class SegmentPropertyMap=boost::typed_identity_property_map<typename GeomTraits::Segment_3>,
          bool cache_primitive=false>
class AABB_segment_primitive :
  public internal::Primitive_caching<typename GeomTraits::Segment_3,Iterator,SegmentPropertyMap,cache_primitive>
{
        // types
        typedef internal::Primitive_caching<typename GeomTraits::Segment_3,Iterator,SegmentPropertyMap,cache_primitive> Base;
public:
        typedef typename GeomTraits::Point_3 Point; // point type
        typedef typename GeomTraits::Segment_3 Datum; // datum type
        typedef Iterator Id; // Id type

        // member data
private:
        Id m_it;
public:
        // constructors
        AABB_segment_primitive() {}
        AABB_segment_primitive(Id it,SegmentPropertyMap pmap=SegmentPropertyMap())
                : m_it(it)
        {
          this->set_primitive(it,pmap);
        }
public:
        Id& id() { return m_it; }
        const Id& id() const { return m_it; }
        typename Base::result_type datum() const {
          return this->get_primitive(m_it);
        }

        /// Returns a point on the primitive
        Point reference_point() const { return datum().source(); }
};

}  // end namespace CGAL

#endif // CGAL_AABB_SEGMENT_PRIMITIVE_H_

