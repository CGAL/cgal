// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Vector_2.h
// source        : include/CGAL/Cartesian/Vector_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr, Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_CARTESIAN_VECTOR_2_H
#define CGAL_CARTESIAN_VECTOR_2_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#include <CGAL/Cartesian/redefine_names_2.h>
#endif

#ifndef CGAL_TWOTUPLE_H
#include <CGAL/Twotuple.h>
#endif // CGAL_TWOTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class _R >
class VectorC2
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef _R                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef VectorC2<R,Cartesian_tag>             Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef VectorC2<R>                           Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  VectorC2();
  VectorC2(const Self &v);
  VectorC2(const Null_vector &);
  VectorC2(const Point_2 &p);
  VectorC2(const Direction_2 &d);
  VectorC2(const FT &hx, const FT &hy, const FT &hw);
  VectorC2(const FT &x, const FT &y);
  ~VectorC2();

  Self            &operator=(const Self &v);

  bool            operator==(const Self &v) const;
  bool            operator!=(const Self &v) const;
  bool            operator==(const Null_vector &) const;
  bool            operator!=(const Null_vector &p) const;

  FT              x() const;
  FT              y() const;
  FT              cartesian(int i) const;
  FT              operator[](int i) const;

  FT              hx() const;
  FT              hy() const;
  FT              hw() const;
  FT              homogeneous(int i) const;

  int             id() const;
  int             dimension() const;

  Self            operator+(const Self &w) const;
  Self            operator-(const Self &w) const;
  Self            operator-() const;
  FT              operator*(const Self &w) const;
  Self            operator/(const FT &c) const;
  Direction_2     direction() const;

  Self            perpendicular(const Orientation &o) const;
  Self            transform(const Aff_transformation_2 &) const;

private:
  _Twotuple<FT>*  ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_VECTOR_2_C
#include <CGAL/Cartesian/Vector_2.C>
#endif // CGAL_CARTESIAN_VECTOR_2_C
#endif 

#endif // CGAL_CARTESIAN_VECTOR_2_H
