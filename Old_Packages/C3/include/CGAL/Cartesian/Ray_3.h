// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Ray_3.h
// source        : include/CGAL/Cartesian/Ray_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================



#ifndef CGAL_CARTESIAN_RAY_3_H
#define CGAL_CARTESIAN_RAY_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#ifndef CGAL_CARTESIAN_TWOTUPLE_H
#include <CGAL/Twotuple.h>
#endif // CGAL_CARTESIAN_TWOTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class _R >
class RayC3
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
  typedef RayC3<R CGAL_CTAG>               Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Direction_3          Direction_3;
  typedef typename R::Line_3               Line_3;
  typedef typename R::Plane_3              Plane_3;
  typedef typename R::Triangle_3           Triangle_3;
  typedef typename R::Segment_3            Segment_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef RayC3<R>                              Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Vector_3_base             Vector_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Line_3_base               Line_3;
  typedef typename R::Plane_3_base              Plane_3;
  typedef typename R::Triangle_3_base           Triangle_3;
  typedef typename R::Segment_3_base            Segment_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  RayC3();
  RayC3(const Self &r);
  RayC3(const Point_3 &sp, const Point_3 &secondp);
  RayC3(const Point_3 &sp, const Direction_3 &d);
  ~RayC3();

  Self        &operator=(const Self &r);

  bool        operator==(const Self &r) const;
  bool        operator!=(const Self &r) const;
  long        id() const;

  Point_3     start() const;
  Point_3     source() const;
  Point_3     second_point() const;
  Point_3     point(int i) const;

  Direction_3 direction() const;
  Line_3      supporting_line() const;
  RayC3       opposite() const;

  RayC3       transform(const Aff_transformation_3 &t) const;

  bool        is_degenerate() const;
  bool        has_on(const Point_3 &p) const;
  bool        collinear_has_on(const Point_3 &p) const;

private:
  _Twotuple< Point_3 > *ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_RAY_3_C
#include <CGAL/Cartesian/Ray_3.C>
#endif // CGAL_CARTESIAN_RAY_3_C
#endif 

#endif
