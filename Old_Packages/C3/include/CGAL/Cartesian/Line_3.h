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
// file          : include/CGAL/Cartesian/Line_3.h
// source        : include/CGAL/Cartesian/Line_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_LINE_3_H
#define CGAL_CARTESIAN_LINE_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#ifndef CGAL_CARTESIAN_TWOTUPLE_H
#include <CGAL/Twotuple.h>
#endif // CGAL_CARTESIAN_TWOTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class _R >
class LineC3
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
  typedef LineC3<R CGAL_CTAG>              Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Direction_3          Direction_3;
  typedef typename R::Plane_3              Plane_3;
  typedef typename R::Ray_3                Ray_3;
  typedef typename R::Segment_3            Segment_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef LineC3<R>                             Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Vector_3_base             Vector_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Plane_3_base              Plane_3;
  typedef typename R::Ray_3_base                Ray_3;
  typedef typename R::Segment_3_base            Segment_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  LineC3();
  LineC3(const Self  &l);
  LineC3(const Point_3 &p, const Point_3 &q);
  LineC3(const Segment_3 &s);
  LineC3(const Ray_3 &r);
  LineC3(const Point_3 &p, const Direction_3 &d);
  ~LineC3();

  Self        &operator=(const Self &l);

  bool        operator==(const Self &l) const;
  bool        operator!=(const Self &l) const;
  long        id() const;

  Plane_3     perpendicular_plane(const Point_3 &p) const;
  Self        opposite() const;

  Point_3     point() const;
  Point_3     second_point() const;
  Point_3     point(int i) const;

  Point_3     projection(const Point_3 &p) const;

  Direction_3 direction() const;

  bool        has_on(const Point_3 &p) const;
  bool        is_degenerate() const;

  Self        transform(const Aff_transformation_3 &t) const;


private:
  _Twotuple< Point_3 >* ptr() const;
  void            new_rep(const Point_3 &p,
                          const Vector_3 &v);
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_LINE_3_C
#include <CGAL/Cartesian/Line_3.C>
#endif // CGAL_CARTESIAN_LINE_3_C
#endif 

#endif // CGAL_CARTESIAN_LINE_3_H
