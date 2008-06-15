// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: 
// $Id: 
// 
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_MS3_POINTMARK_H
#define CGAL_MS3_POINTMARK_H

CGAL_BEGIN_NAMESPACE

template <class K>
class PointMark {

  typedef typename K::Point_3 Point_3;
  typedef PointMark<K> Self;

  Point_3 p;
  bool b;

 public:
  PointMark() : p(0,0,0), b(true) {}
  PointMark(const Self& pm) { p = pm.p; b = pm.b; }
  PointMark(const Point_3& p_, bool b_) : p(p_), b(b_) {}

  Self& operator=(const Self& pm) {
    p = pm.p;
    b = pm.b;
    return *this;
  }

  Self& operator+=(const Self& pm) {
    p = p + (pm.p - CGAL::ORIGIN);
    b = b && pm.b;
    return *this;
  }

  Point_3 point() const {
    return p;
  }

  bool boolean() const {
    return b;
  }

  void set_boolean(bool b_) {
    b = b_;
  }
};

template <typename Kernel>
std::ostream& operator<<(std::ostream& out, 
			 const PointMark<Kernel>& pm) {
  out << pm.point() << "/" << pm.boolean();
  return out;
}

template <typename Kernel>
bool operator==(const PointMark<Kernel>& pm1,
		const PointMark<Kernel>& pm2) {
  return 
    pm1.point() == pm2.point() &&
    pm1.boolean() == pm2.boolean();
}

template <typename Kernel>
const PointMark<Kernel> operator+(const PointMark<Kernel>& pm1,
				  const PointMark<Kernel>& pm2) {
  PointMark<Kernel> ret(pm1);
  ret += pm2;
  return ret;
}

CGAL_END_NAMESPACE
#endif
