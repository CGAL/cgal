// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : Sphere_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

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

      Sphere_3()
      {}

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
