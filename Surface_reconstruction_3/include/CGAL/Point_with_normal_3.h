// Copyright (c) 2007  INRIA (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez

#ifndef CGAL_POINT_WITH_NORMAL_3_H
#define CGAL_POINT_WITH_NORMAL_3_H

#include <CGAL/Point_3.h>
#include <CGAL/Origin.h>
#include <CGAL/Oriented_normal_3.h>

CGAL_BEGIN_NAMESPACE


/// The Point_with_normal_3 class represents a 3D point with:
/// - a position,
/// - a normal (oriented or not).
/// The normal vector is allocated only when needed.
///
/// @heading Is Model for the Concepts: Model of the PointWithNormal_3 concept.
///
/// @heading Parameters:
/// @param Gt   Kernel's geometric traits.

template<class Gt>
class Point_with_normal_3 : public Gt::Point_3
{
// Private types
private:

  typedef typename Gt::Point_3  Base;

// Public types
public:

    typedef Gt Geom_traits; ///< Kernel's geometric traits
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3  Point;  ///< Kernel's Point_3 class.
    typedef Oriented_normal_3<Geom_traits> Normal; ///< Model of OrientedNormal_3 concept.

// Public methods
public:

    /// Point is (0,0,0) by default.
    /// Normal is (0,0,0) by default.
    /// Normal is oriented by default.
    Point_with_normal_3(const Origin& o = ORIGIN)
    : Base(o)
    {
    }
    Point_with_normal_3(FT x, FT y, FT z)
    : Base(x,y,z)
    {
    }
    Point_with_normal_3(const Point& point,
                        const Normal& normal = NULL_VECTOR)
    : Base(point)
    {
      m_normal = normal;
    }

    // Default copy constructor and operator =() are fine.

    /// Compare positions
    bool operator==(const Point_with_normal_3& that)
    {
      return Base::operator==(that);
    }
    bool operator!=(const Point_with_normal_3& that)
    {
      return ! (*this == that);
    }

    /// Get/set normal (vector + orientation).
    const Normal& normal() const { return m_normal; }
    Normal&       normal()       { return m_normal; }

// Data
private:

    Normal  m_normal;
    bool    m_oriented_normal;
};


CGAL_END_NAMESPACE

#endif //CGAL_POINT_WITH_NORMAL_3_H

