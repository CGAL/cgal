// ==========================================================================
//
// Copyright (c) 1998, 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Tetrahedron_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================

#ifndef CGAL_CARTESIAN_TETRAHEDRON_3_H
#define CGAL_CARTESIAN_TETRAHEDRON_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Fourtuple.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class TetrahedronC3
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef _R                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef TetrahedronC3<R CGAL_CTAG>       Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Plane_3              Plane_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef TetrahedronC3<R>                      Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Plane_3_base              Plane_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  TetrahedronC3();
  TetrahedronC3(const Self &t);
  TetrahedronC3(const Point_3 &p,
                const Point_3 &q,
                const Point_3 &r,
                const Point_3 &s);
  ~TetrahedronC3();

  Self &operator=(const Self &t);

  Point_3    vertex(int i) const;
  Point_3    operator[](int i) const;

  bool       operator==(const Self &t) const;
  bool       operator!=(const Self &t) const;
  long       id() const;

  Bbox_3     bbox() const;

  Self       transform(const Aff_transformation_3 &t) const;

  Orientation    orientation() const;
  Oriented_side  oriented_side(const Point_3 &p) const;
  Bounded_side   bounded_side(const Point_3 &p) const;

  bool       has_on_boundary(const Point_3 &p) const;
  bool       has_on_positive_side(const Point_3 &p) const;
  bool       has_on_negative_side(const Point_3 &p) const;
  bool       has_on_bounded_side(const Point_3 &p) const;
  bool       has_on_unbounded_side(const Point_3 &p) const;

  bool       is_degenerate() const;

private:
  _Fourtuple< Point_3 >*   ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Tetrahedron_3.C>
#endif 

#endif // CGAL_CARTESIAN_TETRAHEDRON_3_H
