
// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/point_vector_declarationsH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINT_VECTOR_DECLARATIONSH2_H
#define CGAL_POINT_VECTOR_DECLARATIONSH2_H

#include <CGAL/homogeneous_classes.h>

CGAL_BEGIN_NAMESPACE

template <class R>
inline
PointH2<R>
origin_plus_vector(const VectorH2<R>& v);

template <class R>
inline
PointH2<R>
origin_minus_vector(const VectorH2<R>& v);

template <class R>
inline
VectorH2<R>
point_minus_origin(const PointH2<R>& p);

template <class R>
inline
VectorH2<R>
origin_minus_point(const PointH2<R>& p);

template <class R>
CGAL_KERNEL_INLINE
PointH2<R>
operator-(const PointH2<R>& p, const VectorH2<R>& v);

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator-(const PointH2<R>& p, const PointH2<R>& q);


template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator+(const VectorH2<R>& u, const VectorH2<R>& v);

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator-(const VectorH2<R>& u, const VectorH2<R>& v);

template <class R>
CGAL_KERNEL_INLINE
typename VectorH2<R>::FT
operator*(const VectorH2<R>& u, const VectorH2<R>& v);

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator/(const VectorH2<R>& v, const typename R::RT& f);

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator*(const VectorH2<R>& v, const typename R::RT& f);

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator*(const typename R::RT& f, const VectorH2<R>& v);

CGAL_END_NAMESPACE

#endif // CGAL_POINT_VECTOR_DECLARATIONS_2_H
