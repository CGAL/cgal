// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : Simple_cartesian.lw
// file          : include/CGAL/SimpleCartesian/simple_cartesian_classes.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.7
// revision_date : 11 Aug 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_SIMPLE_CARTESIAN_CLASSES_H
#define CGAL_SIMPLE_CARTESIAN_CLASSES_H
#include <CGAL/basic_classes.h>

CGAL_BEGIN_NAMESPACE

template < class FT > class PointS2;
template < class FT > class VectorS2;
template < class FT > class DirectionS2;
template < class FT > class LineS2;
template < class FT > class RayS2;
template < class FT > class SegmentS2;
template < class FT > class TriangleS2;
template < class FT > class CircleS2;
template < class FT > class ParabolaS2;
template < class FT > class Parabola_arcS2;
template < class PT, class DA > class ConicCPA2;
template < class FT > class Iso_rectangleS2;
template < class FT > class Iso_cuboidS3;
template < class FT > class Aff_transformation_baseS2;
template < class R >  class Aff_transformation_base_2;
template < class FT > class Aff_transformationS2;

template < class FT > class PlaneS3;
template < class FT > class PointS3;
template < class FT > class VectorS3;
template < class FT > class DirectionS3;
template < class FT > class LineS3;
template < class FT > class RayS3;
template < class FT > class SegmentS3;
template < class FT > class TriangleS3;
template < class FT > class TetrahedronS3;
template < class FT > class SphereS3;
template < class FT > class Aff_transformationS3;

template < class FT > class PointCd;

CGAL_END_NAMESPACE

#endif // CGAL_SIMPLE_CARTESIAN_CLASSES_H
