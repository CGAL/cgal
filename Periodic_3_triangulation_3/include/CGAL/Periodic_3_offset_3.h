// Copyright (c) 2006-2009   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_OFFSET_H
#define CGAL_PERIODIC_3_OFFSET_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Cartesian.h>

namespace CGAL { 

class Periodic_3_offset_3 {
  template <class K2>
  friend std::ostream & operator<<(std::ostream &os, 
				   const Periodic_3_offset_3 &off);
  
public:
  Periodic_3_offset_3() : _offx(0), _offy(0), _offz(0) {}
  Periodic_3_offset_3(int x, int y, int z) : _offx(x), _offy(y), _offz(z) {}

  bool is_null() const {
    return ((_offx | _offy | _offz) == 0);
  }

  int& x() { return _offx; }
  int x() const { return _offx; }
  int& y() { return _offy; }
  int y() const { return _offy; }
  int& z() { return _offz; }
  int z() const { return _offz; }

  int &operator[](int i) {
    if (i==0) return _offx;
    if (i==1) return _offy;
    CGAL_triangulation_assertion(i==2);
    return _offz;
  }
  int operator[](int i) const {
    if (i==0) return _offx;
    if (i==1) return _offy;
    CGAL_triangulation_assertion(i==2);
    return _offz;
  }
  void operator+=(const Periodic_3_offset_3 &other) {
    _offx += other._offx;
    _offy += other._offy;
    _offz += other._offz;
  }
  void operator-=(const Periodic_3_offset_3 &other) {
    _offx -= other._offx;
    _offy -= other._offy;
    _offz -= other._offz;
  }
  Periodic_3_offset_3 operator-() const {
    return Periodic_3_offset_3(-_offx,-_offy,-_offz);
  }
  bool operator==(const Periodic_3_offset_3 &other) const {
    return ((_offx == other._offx) &&
	    (_offy == other._offy) &&
	    (_offz == other._offz));
  }
  bool operator!=(const Periodic_3_offset_3 &other) const {
    return  ((_offx != other._offx) ||
	     (_offy != other._offy) ||
	     (_offz != other._offz));
  }
  bool operator<(const Periodic_3_offset_3 &other) const {
    if (_offx != other._offx)
      return (_offx < other._offx);
    else {
      if (_offy != other._offy)
        return (_offy < other._offy);
      else {
        return (_offz < other._offz);
      }
    }
  }

  Periodic_3_offset_3 operator+(const Periodic_3_offset_3 &off2) const {
    return Periodic_3_offset_3(_offx+off2.x(), _offy+off2.y(), _offz+off2.z());
  }
  Periodic_3_offset_3 operator-(const Periodic_3_offset_3 &off2) const {
    return Periodic_3_offset_3(_offx-off2.x(), _offy-off2.y(), _offz-off2.z());
  }
    
private:    
  int _offx, _offy, _offz;
};

template <class K>
inline Point_3<K> operator+(const Point_3<K> &p, const Periodic_3_offset_3 &off) {
  return (off.is_null() ? p : Point_3<K>(p.x()+off.x(), p.y()+off.y(), p.z()+off.z()));
}

inline std::ostream 
&operator<<(std::ostream &os, const Periodic_3_offset_3 &off) {
  if (is_ascii(os))
    os << off.x() << " " << off.y() << " " << off.z();
  else {
    write(os,off.x());
    write(os,off.y());
    write(os,off.z());
  }
  return os;
}

inline std::istream
&operator>>(std::istream &is, Periodic_3_offset_3 &off) {
  int x=0,y=0,z=0;
  if (is_ascii(is))
    is >> x >> y >> z;
  else {
    read(is,x);
    read(is,y);
    read(is,z);
  }
  off = Periodic_3_offset_3(x,y,z);
  return is;
}

} //namespace CGAL

#endif // CGAL_PERIODIC_3_OFFSET_H
