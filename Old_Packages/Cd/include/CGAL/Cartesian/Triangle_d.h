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
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Triangle_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_TRIANGLE_D_H
#define CGAL_CARTESIAN_TRIANGLE_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class TriangleCd CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef TriangleCd<R CGAL_CTAG>          Self;
  typedef typename R::Point_d              Point_d;
  typedef typename R::Vector_d             Vector_d;
  typedef typename R::Plane_d              Plane_d;
  typedef typename R::Aff_transformation_d Aff_transformation_d;
#else
  typedef TriangleCd<R>                         Self;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Vector_d_base             Vector_d;
  typedef typename R::Plane_d_base              Plane_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  TriangleCd();
  TriangleCd(const Self &t);
  TriangleCd(const Point_d &p, const Point_d &q, const Point_d &r);
  ~TriangleCd();

  Self       &operator=(const Self &t);

  Point_d    vertex(int i) const;
  Point_d    operator[](int i) const;

  bool       operator==(const Self &t) const;
  bool       operator!=(const Self &t) const;
  long       id() const;
  int        dimension() const;

  // Bbox_d     bbox() const;

  Self       transform(const Aff_transformation_d &t) const;

  // Only makes sense for 3D
  Plane_d    supporting_plane() const;
  // Makes sense for any dimension, but only implemented in 3D
  bool       has_on(const Point_d &p) const;
  // End of 3D section

  bool       is_degenerate() const;

private:
  _Threetuple< Point_d >*   ptr() const;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRIANGLE_D_H
