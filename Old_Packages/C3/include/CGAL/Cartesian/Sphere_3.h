// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 4.3
// release_date  :  6 Apr 2000
//
// file          : include/CGAL/Cartesian/Sphere_3.h
// package       : C3 (4.3)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_SPHERE_3_H
#define CGAL_CARTESIAN_SPHERE_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Sphere_rep_3.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class SphereC3
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
{
public:
  typedef _R                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef SphereC3<R CGAL_CTAG>                 Self;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Aff_transformation_3      Aff_transformation_3;
#else
  typedef SphereC3<R>                           Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  Sphere_3(const Point_3 &p, const FT &s,
           const Orientation &o = COUNTERCLOCKWISE);
  // Sphere with center p, squared radius s, orientation o
  Sphere_3(const R::Point_3 &p, const R::Point_3 &q,
           const R::Point_3 &r, const R::Point_3 &u);
  // Sphere passing through p,q,r,u, oriented by p, q, r, u
  Sphere_3(const R::Point_3 &p, const R::Point_3 &q, const R::Point_3 &r,
	   const Orientation &o = COUNTERCLOCKWISE);
  // Sphere with great circle passing through p,q,r, oriented by o
  Sphere_3(const Point_3 & p, const Point_3 & q,
           const Orientation &o = COUNTERCLOCKWISE);
  // Sphere with diameter pq and orientation o
  Sphere_3(const Point_3 & p,
           const Orientation& o = COUNTERCLOCKWISE);
  // Sphere centered at p, radius 0, orientation o

  Point  center() const;
  // Returns the center of c
  FT     squared_radius() const;
  // Returns the square of the radius (instead of the radius itself,
  // which would require square roots)
  Orientation orientation() const;
  // Returns the orientation of c

  Self   orthogonal_transform(const Aff_transformation_3 &t) const;
  //! precond: t.is_orthogonal() (*UNDEFINED*)
  // Returns the image of c by t. Since t is orthogonal, the image is
  // always a circle

  bool   is_degenerate() const;
  // A circle is degenerate if its (squared) radius is null or negative
  Self   opposite() const;
  // Returns a circle with opposite orientation

  Oriented_side  oriented_side(const Point_3 &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_POSITIVE_SIDE, R::ON_ORIENTED_BOUNDARY or
  // R::ON_NEGATIVE_SIDE
  bool   has_on_boundary(const Point_3 &p) const
  { return oriented_side(p)==ON_ORIENTED_BOUNDARY; }
  bool   has_on_positive_side(const Point_3 &p) const;
  { return oriented_side(p)==ON_POSITIVE_SIDE; }
  bool   has_on_negative_side(const Point_3 &p) const;
  { return oriented_side(p)==ON_NEGATIVE_SIDE; }

  Bounded_side   bounded_side(const Point_3 &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_BOUNDED_SIDE, R::ON_BOUNDARY or R::ON_UNBOUNDED_SIDE
  bool   has_on_bounded_side(const Point_3 &p) const;
  { return bounded_side(p)==ON_BOUNDED_SIDE; }
  bool   has_on_unbounded_side(const Point_3 &p) const;
  { return bounded_side(p)==ON_UNBOUNDED_SIDE; }

protected:
  _Sphere_repC3<R> *ptr();
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SPHERE_3_H
