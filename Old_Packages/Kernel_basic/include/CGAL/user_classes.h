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
// file          : user_classes.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_USER_CLASSES_H
#define CGAL_USER_CLASSES_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Point_2;

template < class R >
class Vector_2;

template < class R >
class Direction_2;

template < class R >
class Line_2;

template < class R >
class Ray_2;

template < class R >
class Segment_2;

template < class R >
class Triangle_2;

template < class R >
class Iso_rectangle_2;

template < class R >
class Circle_2;

template < class R >
class Aff_transformation_base_2;

template < class R >
class Aff_transformation_2;

template < class R >
class Aff_transformation_3;

template < class R >
class Plane_3;

template < class R >
class Point_3;

template < class R >
class Vector_3;

template < class R >
class Direction_3;

template < class R >
class Line_3;

template < class R >
class Ray_3;

template < class R >
class Segment_3;

template < class R >
class Triangle_3;

template < class R >
class Tetrahedron_3;

template < class R >
class Iso_cuboid_3;

template < class R >
class Sphere_3;

CGAL_END_NAMESPACE

#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Aff_transformation_2.h>

#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Aff_transformation_3.h>

CGAL_BEGIN_NAMESPACE

template < class R>
class Point_d;

CGAL_END_NAMESPACE

#endif  // CGAL_USER_CLASSES_H
