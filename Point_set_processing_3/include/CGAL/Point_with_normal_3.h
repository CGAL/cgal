// Copyright (c) 2007-2009  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez

#ifndef CGAL_POINT_WITH_NORMAL_3_H
#define CGAL_POINT_WITH_NORMAL_3_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Origin.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/property_map.h>


namespace CGAL {

/// \cond SKIP_IN_MANUAL

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


//=========================================================================


/// Property map that accesses the normal vector from a Point_with_normal_3 object
///
/// @heading Is Model for the Concepts:
/// \cgalModels `LvaluePropertyMap`
///
/// @heading Parameters:
/// @param Gt Geometric traits class.

template <class Gt>
struct Normal_of_point_with_normal_map
{
  typedef Normal_of_point_with_normal_map<Gt> Self;

  typedef Point_with_normal_3<Gt> Point_with_normal; ///< Position + normal
  typedef typename Gt::Vector_3 Vector; /// normal

  typedef Point_with_normal key_type;
  typedef Vector value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;

  /// Access a property map element
  value_type& operator[](key_type& k) const { return k.normal(); }

  /// \name Put/get free functions
  /// @{
  friend reference get(const Self&, const key_type& k) { return k.normal(); }
  friend void put(const Self&, key_type& k, const value_type& v) { k.normal() = v; }
  /// @};}
};

/// Free function to create a Normal_of_point_with_normal_map property map.
///
/// @relates Normal_of_point_with_normal_map

template <class Point_with_normal> // Point_with_normal type
Normal_of_point_with_normal_map<
  typename CGAL::Kernel_traits<Point_with_normal>::Kernel>
  make_normal_of_point_with_normal_map(Point_with_normal)
{
  return Normal_of_point_with_normal_map<typename CGAL::Kernel_traits<Point_with_normal>::Kernel>();
}

/// \endcond

} //namespace CGAL


#endif //CGAL_POINT_WITH_NORMAL_3_H
