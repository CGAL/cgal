
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
// release_date  : 2000, October 11
// 
// source        : PV_decl_2.fw
// file          : include/CGAL/point_vector_declarationsH2.h
// package       : H2 (2.13)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.13
// revision_date : 11 Oct 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINT_VECTOR_DECLARATIONSH2_H
#define CGAL_POINT_VECTOR_DECLARATIONSH2_H

#ifndef CGAL_HOMOGENEOUS_CLASSES_H
#include <CGAL/homogeneous_classes.h>
#endif // #ifndef CGAL_HOMOGENEOUS_CLASSES_H


CGAL_BEGIN_NAMESPACE

template <class FT, class RT>
inline
PointH2<FT,RT>
origin_plus_vector(const VectorH2<FT,RT>& v);

template <class FT, class RT>
inline
PointH2<FT,RT>
origin_minus_vector(const VectorH2<FT,RT>& v);

template <class FT, class RT>
inline
VectorH2<FT,RT>
point_minus_origin(const PointH2<FT,RT>& p);

template <class FT, class RT>
inline
VectorH2<FT,RT>
origin_minus_point(const PointH2<FT,RT>& p);

template <class FT, class RT>
CGAL_KERNEL_INLINE
PointH2<FT,RT>
operator-(const PointH2<FT,RT>& p, const VectorH2<FT,RT>& v);

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator-(const PointH2<FT,RT>& p, const PointH2<FT,RT>& q);


template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator+(const VectorH2<FT,RT>& u, const VectorH2<FT,RT>& v);

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator-(const VectorH2<FT,RT>& u, const VectorH2<FT,RT>& v);

template <class FT, class RT>
CGAL_KERNEL_INLINE
FT
operator*(const VectorH2<FT,RT>& u, const VectorH2<FT,RT>& v);

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator/(const VectorH2<FT,RT>& v, const RT& f);

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator*(const VectorH2<FT,RT>& v, const RT& f);

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator*(const RT& f, const VectorH2<FT,RT>& v);



CGAL_END_NAMESPACE


#endif // CGAL_POINT_VECTOR_DECLARATIONS_2_H
