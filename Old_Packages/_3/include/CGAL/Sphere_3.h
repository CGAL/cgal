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
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_SPHERE_3_H
#define CGAL_SPHERE_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/SphereH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/Sphere_3.h>
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/SimpleCartesian/SphereS3.h>
#endif // CGAL_SIMPLE_CARTESIAN_H


namespace CGAL {


template <class R_>
class Sphere_3 : public R_::Sphere_3_base

{
  public:

  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Sphere_3_base  RSphere_3;

      Sphere_3()
      {}

      Sphere_3(const RSphere_3& s)
      : RSphere_3(s)
      {}

      Sphere_3(const CGAL::Point_3<R>& p, const FT& sq_rad,
               const Orientation& o = COUNTERCLOCKWISE)
       : RSphere_3(p, sq_rad, o)
      {}

      Sphere_3(const CGAL::Point_3<R>& p, const CGAL::Point_3<R>& q,
               const CGAL::Point_3<R>& r, const CGAL::Point_3<R>& u)
       : RSphere_3(p, q, r, u)
      {}

      Sphere_3(const CGAL::Point_3<R>& p, const CGAL::Point_3<R>& q,
               const CGAL::Point_3<R>& r,
               const Orientation& o = COUNTERCLOCKWISE)
       : RSphere_3(p, q, r, o)
      {}

      Sphere_3(const CGAL::Point_3<R>&  p, const CGAL::Point_3<R>&  q,
               const Orientation& o = COUNTERCLOCKWISE)
       : RSphere_3(p, q, o)
      {}

      Sphere_3(const CGAL::Point_3<R>&  p,
               const Orientation& o = COUNTERCLOCKWISE)
       : RSphere_3(p, o)
      {}

};

} // namespace CGAL

#endif // CGAL_SPHERE_3_H
