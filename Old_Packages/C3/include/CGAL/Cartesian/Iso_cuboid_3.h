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
// file          : include/CGAL/Cartesian/Iso_cuboid_3.h
// package       : C3 (4.3)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_ISO_CUBOID_3_H
#define CGAL_CARTESIAN_ISO_CUBOID_3_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class _R >
class Iso_cuboidC3
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
  typedef Iso_cuboidC3<R CGAL_CTAG>        Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef Iso_cuboidC3<R>                  Self;
  typedef typename R::Point_3_base         Point_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

public:
  Iso_cuboidC3();
  Iso_cuboidC3(const Self& );
  Iso_cuboidC3(const Point_3& p, const Point_3& q);
  ~Iso_cuboidC3();

  Self&        operator=(const Self& r);

  bool         operator==(const Self& s) const;
  bool         operator!=(const Self& s) const;
  bool         identical(const Self& s) const;
  int          id() const;

  Point_3         min() const;
  Point_3         max() const;
  Point_3         vertex(int i) const;
  Point_3         operator[](int i) const;

  Self         transform(const Aff_transformation_3& t) const;
  Bounded_side bounded_side(const Point_3& p) const;
  bool         has_on(const Point_3& p) const;
  bool         has_on_boundary(const Point_3& p) const;
  bool         has_on_bounded_side(const Point_3& p) const;
  bool         has_on_unbounded_side(const Point_3& p) const;
  bool         is_degenerate() const;
  Bbox_3       bbox() const;
  FT           xmin() const;
  FT           ymin() const;
  FT           zmin() const;
  FT           xmax() const;
  FT           ymax() const;
  FT           zmax() const;

protected:
  _Twotuple<Point_3>*  ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Iso_cuboid_3.C>
#endif

#endif // CGAL_CARTESIAN_ISO_CUBOID_3_C
