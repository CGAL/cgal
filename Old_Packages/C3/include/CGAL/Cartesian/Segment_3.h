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
// file          : include/CGAL/Cartesian/Segment_3.h
// source        : include/CGAL/Cartesian/Segment_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_SEGMENT_3_H
#define CGAL_CARTESIAN_SEGMENT_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#ifndef CGAL_CARTESIAN_TWOTUPLE_H
#include <CGAL/Twotuple.h>
#endif // CGAL_CARTESIAN_TWOTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class _R >
class SegmentC3
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
  typedef SegmentC3<R CGAL_CTAG>           Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Direction_3          Direction_3;
  typedef typename R::Line_3               Line_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef SegmentC3<R>                          Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Line_3_base               Line_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  SegmentC3();
  SegmentC3(const Self  &s);
  SegmentC3(const Point_3 &sp, const Point_3 &ep);
  ~SegmentC3();

  Self        &operator=(const Self &s);

  bool        has_on(const Point_3 &p) const;
  bool        collinear_has_on(const Point_3 &p) const;

  bool        operator==(const Self &s) const;
  bool        operator!=(const Self &s) const;
  long        id() const;

  Point_3     start() const;
  Point_3     end() const;

  Point_3     source() const;
  Point_3     target() const;

  Point_3     min() const;
  Point_3     max() const;
  Point_3     vertex(int i) const;
  Point_3     point(int i) const;
  Point_3     operator[](int i) const;

  FT          squared_length() const;

  Direction_3 direction() const;
  Line_3      supporting_line() const;
  Self        opposite() const;
  Self        transform(const Aff_transformation_3 &t) const;

  bool        is_degenerate() const;
  Bbox_3      bbox() const;

private:
  _Twotuple< Point_3 >*   ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_SEGMENT_3_C
#include <CGAL/Cartesian/Segment_3.C>
#endif // CGAL_CARTESIAN_SEGMENT_3_C
#endif 

#endif // CGAL_CARTESIAN_SEGMENT_3_C
