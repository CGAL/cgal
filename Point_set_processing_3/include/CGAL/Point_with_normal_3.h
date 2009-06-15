// Copyright (c) 2007-2009  INRIA (France).
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
#include <CGAL/Vector_3.h>
#include <CGAL/Origin.h>

CGAL_BEGIN_NAMESPACE


/// The Point_with_normal_3 class represents a 3D point with:
/// - a position,
/// - a normal (oriented).
///
/// @heading Parameters:
/// @param Gt       Geometric traits class.

template<class Gt>
class Point_with_normal_3 : public Gt::Point_3
{
// Private types
private:

  typedef typename Gt::Point_3  Base;

// Public types
public:

    typedef Gt Geom_traits; ///< Geometric traits class
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::RT RT;
    typedef typename Geom_traits::Point_3  Point;  ///< typedef to Geom_traits::Point_3
    typedef typename Geom_traits::Vector_3 Vector; ///< typedef to Geom_traits::Vector_3

// Public methods
public:

    /// Point is (0,0,0) by default.
    /// Normal is (0,0,0) by default.
    Point_with_normal_3(const Origin& o = ORIGIN)
    : Base(o)
    {
    }
    Point_with_normal_3(FT x, FT y, FT z,
                        const Vector& normal = NULL_VECTOR)
    : Base(x,y,z),
      m_normal(normal)
    {
    }
    Point_with_normal_3(RT hx, RT hy, RT hz, RT hw,
                        const Vector& normal = NULL_VECTOR)
    : Base(hx,hy,hz,hw),
      m_normal(normal)
    {
    }
    Point_with_normal_3(const Point& point,
                        const Vector& normal = NULL_VECTOR)
    : Base(point),
      m_normal(normal)
    {
    }

    /// Copy constructor
    Point_with_normal_3(const Point_with_normal_3& pwn)
    : Base(pwn),
      m_normal(pwn.normal())
    {
    }
    template <class K>
    Point_with_normal_3(const Point_with_normal_3<K>& pwn)
    : Base(pwn),
      m_normal(pwn.normal())
    {
    }
    /// Operator =()
    Point_with_normal_3& operator=(const Point_with_normal_3& pwn)
    {
      Base::operator=(pwn);
      m_normal = pwn.normal();
      return *this;
    }

    /// Gets/sets position.
    const Point& position() const { return *this; }
    Point&       position()       { return *this; }

    /// Gets/sets normal.
    const Vector& normal() const { return m_normal; }
    Vector&       normal()       { return m_normal; }

// Data
private:

    Vector  m_normal;
};


CGAL_END_NAMESPACE


#endif //CGAL_POINT_WITH_NORMAL_3_H

