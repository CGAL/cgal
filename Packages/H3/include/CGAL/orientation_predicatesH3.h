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
// file          : orientation_predicatesH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_ORIENTATION_PREDICATESH3_H
#define CGAL_ORIENTATION_PREDICATESH3_H

#include <CGAL/PVDH3.h>
#include <CGAL/predicates_on_rtH2.h>
#include <CGAL/predicates/sign_of_determinant.h>

CGAL_BEGIN_NAMESPACE

template < class R>
CGAL_KERNEL_INLINE
Orientation
orientation( const PointH3<R>& p,
             const PointH3<R>& q,
             const PointH3<R>& r,
             const PointH3<R>& s)
{
  // Two rows are switched, because of the homogeneous column.
  return (Orientation) sign_of_determinant4x4( p.hx(), p.hy(), p.hz(), p.hw(),
                                               r.hx(), r.hy(), r.hz(), r.hw(),
                                               q.hx(), q.hy(), q.hz(), q.hw(),
                                               s.hx(), s.hy(), s.hz(), s.hw());
}

template < class R>
inline
bool
are_positive_oriented( const PointH3<R>& p,
                       const PointH3<R>& q,
                       const PointH3<R>& r,
                       const PointH3<R>& s)
{ return (orientation(p,q,r,s) == POSITIVE); }

template < class R>
inline
bool
are_negative_oriented( const PointH3<R>& p,
                       const PointH3<R>& q,
                       const PointH3<R>& r,
                       const PointH3<R>& s)
{ return (orientation(p,q,r,s) == NEGATIVE); }

template < class R>
inline
bool
coplanar( const PointH3<R>& p,
          const PointH3<R>& q,
          const PointH3<R>& r,
          const PointH3<R>& s)
{ return (orientation(p,q,r,s) == COPLANAR); }


template <class R>
Orientation
coplanar_orientation(const PointH3<R>& p,
                     const PointH3<R>& q,
                     const PointH3<R>& r,
                     const PointH3<R>& s)
// p,q,r,s supposed to be coplanar
// p, q, r supposed to be non collinear
// tests whether s is on the same side of p,q as r
// returns :
// COLLINEAR if pqs collinear
// POSITIVE if pqr and pqs have the same orientation
// NEGATIVE if pqr and pqs have opposite orientations
{
  CGAL_kernel_precondition( coplanar( p, q, r, s));
  // p, q, r define a plane P:
  CGAL_kernel_precondition( !collinear( p, q, r));
  // compute orientation of p,q,s in this plane P:
  Orientation save;
  if ( (save = orientationH2( p.hy(), p.hz(), p.hw(),
                              q.hy(), q.hz(), q.hw(),
                              r.hy(), r.hz(), r.hw())) != COLLINEAR)
  { return
      static_cast<Orientation>(
        static_cast<int>( save)
      * static_cast<int>( orientationH2( p.hy(), p.hz(), p.hw(),
                                         q.hy(), q.hz(), q.hw(),
                                         s.hy(), s.hz(), s.hw())) );
  }
  if ( (save = orientationH2( p.hx(), p.hz(), p.hw(),
                              q.hx(), q.hz(), q.hw(),
                              r.hx(), r.hz(), r.hw())) != COLLINEAR)
  { return
      static_cast<Orientation>(
        static_cast<int>( save)
      * static_cast<int>( orientationH2( p.hx(), p.hz(), p.hw(),
                                         q.hx(), q.hz(), q.hw(),
                                         s.hx(), s.hz(), s.hw())) );
  }
  if ( (save = orientationH2( p.hx(), p.hy(), p.hw(),
                              q.hx(), q.hy(), q.hw(),
                              r.hx(), r.hy(), r.hw())) != COLLINEAR)
  { return
      static_cast<Orientation>(
        static_cast<int>( save)
      * static_cast<int>( orientationH2( p.hx(), p.hy(), p.hw(),
                                         q.hx(), q.hy(), q.hw(),
                                         s.hx(), s.hy(), s.hw())) );
  }
  CGAL_kernel_assertion( false);
  return COLLINEAR;
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
Orientation
coplanar_orientation(const PointH3<R>& p,
                     const PointH3<R>& q,
                     const PointH3<R>& r)
// Returns an Orientation which is coherent for all (p,q,r) chosen in a same
// plane.
{
  Orientation oxy_pqr = orientationH2(p.hx(), p.hy(), p.hw(),
	                              q.hx(), q.hy(), q.hw(),
				      r.hx(), r.hy(), r.hw());
  if (oxy_pqr != COLLINEAR)
      return oxy_pqr;

  Orientation oyz_pqr = orientationH2(p.hy(), p.hz(), p.hw(),
	                              q.hy(), q.hz(), q.hw(),
				      r.hy(), r.hz(), r.hw());
  if (oyz_pqr != COLLINEAR)
      return oyz_pqr;

  return orientationH2(p.hx(), p.hz(), p.hw(),
	               q.hx(), q.hz(), q.hw(),
		       r.hx(), r.hz(), r.hw());
}

CGAL_END_NAMESPACE

#endif // CGAL_ORIENTATION_PREDICATESH3_H
