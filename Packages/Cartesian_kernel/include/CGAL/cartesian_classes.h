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
// file          : include/CGAL/cartesian_classes.h
// revision      : $Revision$
// revision_date : $Date$
// authors       : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_CLASSES_H
#define CGAL_CARTESIAN_CLASSES_H

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_CFG_NO_ADVANCED_KERNEL

// The following scheme is proposed by Stefan Schirra
// It does not use partial specialization
// It is a partially extendible kernel, but you have to redefine and
// reinstantiate all the member classes in a new kernel
// even if only one differs

template < class R > class PointC2;
template < class R > class VectorC2;
template < class R > class DirectionC2;
template < class R > class LineC2;
template < class R > class RayC2;
template < class R > class SegmentC2;
template < class R > class TriangleC2;
template < class R > class CircleC2;
template < class R > class Data_accessorC2;
template < class PT, class DA > class ConicCPA2;
template < class R > class Iso_rectangleC2;
template < class R > class Aff_transformationC2;

template < class R > class PlaneC3;
template < class R > class PointC3;
template < class R > class VectorC3;
template < class R > class DirectionC3;
template < class R > class LineC3;
template < class R > class RayC3;
template < class R > class SegmentC3;
template < class R > class TriangleC3;
template < class R > class TetrahedronC3;
template < class R > class Iso_cuboidC3;
template < class R > class SphereC3;
template < class R > class Aff_transformationC3;

template < class R > class PointCd;

#endif // CGAL_CFG_NO_ADVANCED_KERNEL

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CLASSES_H
