// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// source        : ddim_points.fw
// file          : PointCd.h
// revision      : 2.2.3
// revision_date : 14 Sep 1999 
// author(s)     : Sven Schoenherr
//                 Bernd Gaertner
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
// ============================================================================
 
#ifndef CGAL_POINTCD_H
#define CGAL_POINTCD_H

#ifndef D_TUPLE_H
#include <CGAL/d_tuple.h>
#endif // D_TUPLE_H

CGAL_BEGIN_NAMESPACE

template < class R >
class DACd;

template < class R >
class PointCd
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<R,Cartesian_tag>
#endif
  : public Handle
{
  friend class DACd<FT>;

public:
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef PointCd<R,Cartesian_tag>              Self;
#else
  typedef PointCd<R>                            Self;
#endif

  PointCd ();
  PointCd (int dim, const Origin&);
  PointCd (const Self &p);

  template <class InputIterator>
  PointCd (int dim, InputIterator first, InputIterator last);

  ~PointCd();

  Self &operator=(const Self &p);

  bool operator==(const Self &p) const;
  bool operator!=(const Self &p) const;

  int id() const;

  FT homogeneous (int i) const;
  FT cartesian (int i) const;
  FT operator[] (int i) const;
  const FT* begin() const;
  const FT* end() const;

  int dimension () const;

  private:
   const _d_tuple<FT>* ptr() const;
};

CGAL_END_NAMESPACE

#endif // CGAL_POINTCD_H
