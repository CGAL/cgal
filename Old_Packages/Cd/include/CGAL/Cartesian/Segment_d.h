
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Segment_d.h
// source        : include/CGAL/Cartesian/Segment_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_SEGMENT_D_H
#define CGAL_CARTESIAN_SEGMENT_D_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#include <CGAL/Cartesian/redefine_names_d.h>
#endif

#ifndef CGAL_CARTESIAN_TWOTUPLE_H
#include <CGAL/Twotuple.h>
#endif // CGAL_CARTESIAN_TWOTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class _R >
class SegmentCd
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
  typedef SegmentCd<R CGAL_CTAG>           Self;
  typedef typename R::Point_d              Point_d;
  typedef typename R::Direction_d          Direction_d;
  typedef typename R::Line_d               Line_d;
  typedef typename R::Aff_transformation_d Aff_transformation_d;
#else
  typedef SegmentCd<R>                          Self;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Direction_d_base          Direction_d;
  typedef typename R::Line_d_base               Line_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  SegmentCd(int d = 0);
  SegmentCd(const Self  &s);
  SegmentCd(const Point_d &sp, const Point_d &ep);
  ~SegmentCd();

  Self        &operator=(const Self &s);

  bool        has_on(const Point_d &p) const;
  bool        collinear_has_on(const Point_d &p) const;

  bool        operator==(const Self &s) const;
  bool        operator!=(const Self &s) const;
  long        id() const;

  Point_d     start() const;
  Point_d     end() const;

  Point_d     source() const;
  Point_d     target() const;

  Point_d     min() const;
  Point_d     max() const;
  Point_d     vertex(int i) const;
  Point_d     point(int i) const;
  Point_d     operator[](int i) const;

  FT          squared_length() const;

  Direction_d direction() const;
  Line_d      supporting_line() const;
  Self        opposite() const;
  Self        transform(const Aff_transformation_d &t) const;

  bool        is_degenerate() const;
  Bbox_d      bbox() const;

private:
  _Twotuple< Point_d >*   ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_SEGMENT_D_C
#include <CGAL/Cartesian/Segment_d.C>
#endif // CGAL_CARTESIAN_SEGMENT_D_C
#endif 

#endif // CGAL_CARTESIAN_SEGMENT_D_C
