// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_SPHERE_3_H
#define CGAL_SPHERE_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Sphere_3 : public R_::Kernel_base::Sphere_3
{
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Kernel_base::Sphere_3  RSphere_3;
public:
  typedef          R_                       R;

      Sphere_3() {}

      Sphere_3(const RSphere_3& s)
      : RSphere_3(s) {}

      Sphere_3(const Point_3& p, const FT& sq_rad,
               const Orientation& o = COUNTERCLOCKWISE)
       : RSphere_3(p, sq_rad, o) {}

      Sphere_3(const Point_3& p, const Point_3& q,
               const Point_3& r, const Point_3& u)
       : RSphere_3(p, q, r, u) {}

      Sphere_3(const Point_3& p, const Point_3& q, const Point_3& r,
               const Orientation& o = COUNTERCLOCKWISE)
       : RSphere_3(p, q, r, o) {}

      Sphere_3(const Point_3& p, const Point_3&  q,
               const Orientation& o = COUNTERCLOCKWISE)
       : RSphere_3(p, q, o) {}

      Sphere_3(const Point_3& p, const Orientation& o = COUNTERCLOCKWISE)
       : RSphere_3(p, o) {}
};

CGAL_END_NAMESPACE

#endif // CGAL_SPHERE_3_H
