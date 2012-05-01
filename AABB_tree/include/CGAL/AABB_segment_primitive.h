// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Pierre Alliez, Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_SEGMENT_PRIMITIVE_H_
#define CGAL_AABB_SEGMENT_PRIMITIVE_H_

namespace CGAL {

template <class GeomTraits, class Iterator>
class AABB_segment_primitive
{
        // types
public:
        typedef typename GeomTraits::Point_3 Point; // point type
        typedef typename GeomTraits::Segment_3 Datum; // datum type
        typedef Iterator Id; // Id type

        // member data
private:
        Id m_it;
        Datum m_datum;

public:
        // constructors
        AABB_segment_primitive() {}
        AABB_segment_primitive(Id it)
                : m_it(it)
        {
                m_datum = *it; // copy segment
        }
        AABB_segment_primitive(const AABB_segment_primitive& primitive)
        {
                m_it = primitive.id();
                m_datum = primitive.datum();
        }
public:
        Id& id() { return m_it; }
        const Id& id() const { return m_it; }
        Datum& datum() { return m_datum; }
        const Datum& datum() const { return m_datum; }

        /// Returns a point on the primitive
        Point reference_point() const { return m_datum.source(); }
};

}  // end namespace CGAL

#endif // CGAL_AABB_SEGMENT_PRIMITIVE_H_

