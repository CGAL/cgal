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
// file          : include/CGAL/Cartesian/Point_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri and Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_POINT_3_H
#define CGAL_CARTESIAN_POINT_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle_for< Threetuple<typename R_::FT> >
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef PointC3<R CGAL_CTAG>             Self;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef PointC3<R>                       Self;
  typedef typename R::Vector_3_base        Vector_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  PointC3();
  PointC3(const Origin &o);
  PointC3(const Self &p);
  PointC3(const Vector_3 &v);
  PointC3(const FT &x, const FT &y, const FT &z);
  PointC3(const FT &x, const FT &y, const FT &z, const FT &hw);
  ~PointC3();

  bool        operator==(const Self &p) const;
  bool        operator!=(const Self &p) const;

  FT          x() const;
  FT          y() const;
  FT          z() const;

  FT          hx() const;
  FT          hy() const;
  FT          hz() const;
  FT          hw() const;

  FT          cartesian(int i) const;
  FT          operator[](int i) const;

  FT          homogeneous(int i) const;

  int         dimension() const;
  Bbox_3      bbox() const;

  Self        transform( const Aff_transformation_3 &) const;

};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Point_3.C>
#endif 

#endif // CGAL_CARTESIAN_POINT_3_H
