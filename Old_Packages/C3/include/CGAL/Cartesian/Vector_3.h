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
// file          : include/CGAL/Cartesian/Vector_3.h
// source        : include/CGAL/Cartesian/Vector_3.h
// revision      : $Revision$
// revision_date : $Date$
// author        : Andreas Fabri
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_VECTOR_3_H
#define CGAL_CARTESIAN_VECTOR_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#ifndef CGAL_CARTESIAN_THREETUPLE_H
#include <CGAL/Threetuple.h>
#endif // CGAL_CARTESIAN_THREETUPLE_H

CGAL_BEGIN_NAMESPACE

template < class _R >
class VectorC3
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
  typedef VectorC3<R CGAL_CTAG>            Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Direction_3          Direction_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef VectorC3<R>                           Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  VectorC3();
  VectorC3(const Self &v);
  VectorC3(const Null_vector &);
  VectorC3(const Point_3 &p);
  VectorC3(const Direction_3 &p);
  VectorC3(const FT &x, const FT &y, const FT &z);
  VectorC3(const FT &x, const FT &y, const FT &z, const FT &w);
  ~VectorC3();

  Self            &operator=(const Self &v);

  bool            operator==(const Self &p) const;
  bool            operator!=(const Self &p) const;

  bool            operator==(const Null_vector &) const;
  bool            operator!=(const Null_vector &) const;

  long            id() const;

  FT              x() const;
  FT              y() const;
  FT              z() const;
  FT              cartesian(int i) const;
  FT              operator[](int i) const;

  FT              hx() const;
  FT              hy() const;
  FT              hz() const;
  FT              hw() const;
  FT              homogeneous(int i) const;

  int             dimension() const;

  Self            operator+(const Self &w) const;
  Self            operator-(const Self &w) const;
  Self            operator-() const;
  FT              operator*(const Self &w) const;
  Self            operator/(const FT &c) const;
  Direction_3 direction() const;
  Self            transform(const Aff_transformation_3 &) const;

private:
  _Threetuple<FT>*   ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_VECTOR_3_C
#include <CGAL/Cartesian/Vector_3.C>
#endif // CGAL_CARTESIAN_VECTOR_3_C
#endif 

#endif // CGAL_CARTESIAN_VECTOR_3_H
