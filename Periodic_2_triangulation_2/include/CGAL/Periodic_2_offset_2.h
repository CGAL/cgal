// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_OFFSET_H
#define CGAL_PERIODIC_2_OFFSET_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Cartesian.h>

namespace CGAL
{

/// The class Periodic_2_offset_2 is a model of the concept Periodic_2Offset_2.
class Periodic_2_offset_2
{
  //  template <class K2>
  //  friend std::ostream & operator<<(std::ostream &os,
  //				   const Periodic_2_offset_2 &off);

public:
  /// Default constructor.
  Periodic_2_offset_2() : _offx(0), _offy(0) {}
  /// Constructs the offset (x,y).
  Periodic_2_offset_2(int x, int y) : _offx(x), _offy(y) {}

  /// Returns true if o is equal to (0,0).
  inline bool is_null() const
  {
    return is_zero();
  }
  /// Returns true if o is equal to (0,0).
  inline bool is_zero() const
  {
    return ((_offx | _offy) == 0);
  }

  /// Return the x-entry of o.
  int& x()
  {
    return _offx;
  }
  /// Return the x-entry of o.
  int x() const
  {
    return _offx;
  }
  /// Return the y-entry of o.
  int& y()
  {
    return _offy;
  }
  /// Return the y-entry of o.
  int y() const
  {
    return _offy;
  }

  /// Return the i-th entry of o.
  int &operator[](int i)
  {
    if (i == 0) return _offx;
    CGAL_triangulation_assertion(i == 1);
    return _offy;
  }
  /// Return the i-th entry of o.
  int operator[](int i) const
  {
    if (i == 0) return _offx;
    CGAL_triangulation_assertion(i == 1);
    return _offy;
  }
  /// Add o' to o using vector addition.
  void operator+=(const Periodic_2_offset_2 &other)
  {
    _offx += other._offx;
    _offy += other._offy;
  }
  /// Subtract o' from o using vector subtraction.
  void operator-=(const Periodic_2_offset_2 &other)
  {
    _offx -= other._offx;
    _offy -= other._offy;
  }
  /// Return the negative vector of o.
  Periodic_2_offset_2 operator-() const
  {
    return Periodic_2_offset_2(-_offx, -_offy);
  }
  /// Return true if o' and o represent the same vector.
  bool operator==(const Periodic_2_offset_2 &other) const
  {
    return ((_offx == other._offx) &&
            (_offy == other._offy));
  }
  /// Return true if o' and o do not represent the same vector.
  bool operator!=(const Periodic_2_offset_2 &other) const
  {
    return  ((_offx != other._offx) ||
             (_offy != other._offy));
  }
  /// Compare o and o' lexicographically.
  bool operator<(const Periodic_2_offset_2 &other) const
  {
    if (_offx != other._offx)
      return (_offx < other._offx);
    else
      {
        return (_offy < other._offy);
      }
  }

  /// Return the vector sum of o and o'.
  Periodic_2_offset_2 operator+(const Periodic_2_offset_2 &off2) const
  {
    return Periodic_2_offset_2(_offx + off2.x(), _offy + off2.y());
  }
  /// Return the vector difference of o and o'.
  Periodic_2_offset_2 operator-(const Periodic_2_offset_2 &off2) const
  {
    return Periodic_2_offset_2(_offx - off2.x(), _offy - off2.y());
  }

private:
  int _offx, _offy;
};

template <class K>
inline typename K::Point_2 operator+(const typename K::Point_2 &p, const Periodic_2_offset_2 &off)
{
  return (off.is_null() ? p : typename K::Point_2(p.x() + off.x(), p.y() + off.y()));
}

/// Inputs an Periodic_2_offset_2 from is.
inline std::ostream
&operator<<(std::ostream &os, const Periodic_2_offset_2 &off)
{
  if (is_ascii(os))
    os << off.x() << " " << off.y();
  else
    {
      write(os, off.x());
      write(os, off.y());
    }
  return os;
}

/// Outputs an Periodic_2_offset_2 to os.
inline std::istream
&operator>>(std::istream &is, Periodic_2_offset_2 &off)
{
  int x = 0, y = 0;
  if (is_ascii(is))
    is >> x >> y;
  else
    {
      read(is, x);
      read(is, y);
    }
  off = Periodic_2_offset_2(x, y);
  return is;
}

} //namespace CGAL

#endif // CGAL_PERIODIC_2_OFFSET_H
