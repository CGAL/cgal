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
// $URL: $
// $Id: $
//
//
// Author(s)     : Pierre Alliez, Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef AABB_TRIANGLE_PRIMITIVE_H_
#define AABB_TRIANGLE_PRIMITIVE_H_

namespace CGAL {

template <class GeomTraits, class Iterator>
class AABB_triangle_primitive
{
public:
        // types
        typedef Iterator Id; // Id type
        typedef typename GeomTraits::FT FT; // field number type
        typedef typename GeomTraits::Point_3 Point; // point type
        typedef typename GeomTraits::Triangle_3 Datum; // datum type

private:
        // member data
        Id m_it; // iterator
        Datum m_datum; // 3D triangle

        // constructor
public:
        AABB_triangle_primitive(Id it)
                : m_it(it)
        {
                m_datum = *it; // copy triangle
        }
public:
        const Datum& datum() const { return m_datum; }
        Datum& datum() { return m_datum; }
        Id id() { return m_it; }

        /// Returns the x/y/z reference coordinate for sorting
        /// here simply the first vertex of the triangle
        const FT xref() const { return m_datum.vertex(0).x(); }
        const FT yref() const { return m_datum.vertex(0).y(); }
        const FT zref() const { return m_datum.vertex(0).z(); }

        /// Returns a point on the primitive
        Point point_on() const { m_datum.vertex(0); }
};

}  // end namespace CGAL


#endif // AABB_TRIANGLE_PRIMITIVE_H_

