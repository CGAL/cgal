// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Kernel_d/Aff_transformation_d.h
// package       : Kernel_d
// chapter       : Basic
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// implementation: affine transformation interface type
// ============================================================================

#ifndef CGAL_AFF_TRANSFORMATION_D_H
#define CGAL_AFF_TRANSFORMATION_D_H

CGAL_BEGIN_NAMESPACE

template <class pR>
class Aff_transformation_d : public pR::Aff_transformation_d_base
{ public:
  typedef typename pR::Aff_transformation_d_base Base;
  typedef Aff_transformation_d<pR>               Self;
  typedef pR R;
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;

  Aff_transformation_d(int d = 0) : Base(d) {}
  Aff_transformation_d(int d, Identity_transformation tag) : Base(d,tag) {}
  Aff_transformation_d(const typename LA::Matrix& M) : Base(M) {}
  Aff_transformation_d(Translation tag, const Vector_d<R>& v) : Base(tag,v) {}
  Aff_transformation_d(int d, Scaling tag, const RT& num, const RT& den) 
    : Base(d,tag,num,den) {}
  Aff_transformation_d(int d, Rotation tag, 
		       const RT& sin_num, const RT& cos_num, 
                       const RT& den, int e1 = 0, int e2 = 1)
    : Base(d,tag,sin_num,cos_num,den,e1,e2) {}
  Aff_transformation_d(int d, Rotation tag, const Direction_d<R>& dir, 
                       const RT& num, const RT& den, 
                       int e1 = 0, int e2 = 1)
    : Base(d,tag,dir,num,den,e1,e2) {}
  Aff_transformation_d(const Base& a) : Base(a) {}
  Aff_transformation_d(const Self& a) : Base(a) {}

#ifndef CGAL_SIMPLE_INTERFACE

  template <typename Forward_iterator>
  Aff_transformation_d(Scaling tag,
    Forward_iterator start, Forward_iterator end) : Base(tag,start,end) {}

#else 
#define FIXATD(I) \
Aff_transformation_d(Scaling tag, I start, I end) : Base(tag,start,end) {}
FIXATD(int*)
FIXATD(const int*)
FIXATD(RT*)
FIXATD(const RT*)
#undef FIXATD
#endif // CGAL_SIMPLE_INTERFACE

  Self operator*(const Self& a)
  { return Base::operator*(a); }
  Self inverse() const { return Base::inverse(); }

  bool operator==(const Self& a) const
  { return Base::operator==(a); }
  bool operator!=(const Self& a) const
  { return Base::operator!=(a); } 

};

CGAL_END_NAMESPACE
#endif //CGAL_AFF_TRANSFORMATION_D_H
