// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.7-branch/AABB_tree/include/CGAL/AABB_segment_primitive.h $
// $Id: AABB_segment_primitive.h 52790 2009-10-29 10:42:53Z palliez $
//
//
// Author(s) : Pierre Alliez, Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_SEGMENT_2_PRIMITIVE_H
#define CGAL_AABB_SEGMENT_2_PRIMITIVE_H

namespace CGAL {

template <class GeomTraits, class Iterator_, class ContainerType>
class AABB_segment_2_primitive {
    // types
public:
    typedef typename GeomTraits::Point_2 Point; // point type
    typedef typename GeomTraits::Segment_2 Datum; // datum type
    typedef ContainerType Container;
    typedef Iterator_ Id; // Id type

    // member data
private:
    Id m_it;
    Datum m_datum;

public:
    // constructors
    AABB_segment_2_primitive() {}
    AABB_segment_2_primitive(Id it)
        : m_it(it) {
        m_datum = *it; // copy segment
    }
    AABB_segment_2_primitive(const AABB_segment_2_primitive &primitive) {
        m_it = primitive.id();
        m_datum = primitive.datum();
    }
public:
    Id &id() {
        return m_it;
    }
    const Id &id() const {
        return m_it;
    }
    Datum &datum() {
        return m_datum;
    }
    const Datum &datum() const {
        return m_datum;
    }

    /// Returns a point on the primitive
    Point reference_point() const {
        return m_datum.source();
    }
};

} // namespace CGAL

#endif
