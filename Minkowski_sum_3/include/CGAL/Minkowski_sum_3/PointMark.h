// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_MS3_POINTMARK_H
#define CGAL_MS3_POINTMARK_H

#include <CGAL/license/Minkowski_sum_3.h>


namespace CGAL {

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

} //namespace CGAL
#endif
