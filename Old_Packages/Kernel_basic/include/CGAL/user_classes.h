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
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_USER_CLASSES_H
#define CGAL_USER_CLASSES_H

#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Conic_2.h>
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

#ifdef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#  define CGAL_ITERATOR_TRAITS_POINTER_SPEC_2(K) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Point_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Vector_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Direction_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Line_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Segment_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Ray_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Iso_rectangle_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Triangle_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Circle_2< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Aff_transformation_2< K >)
 
#  define CGAL_ITERATOR_TRAITS_POINTER_SPEC_3(K) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Point_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Vector_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Direction_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Plane_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Line_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Segment_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Ray_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Triangle_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Tetrahedron_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Sphere_3< K >) \
     CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Aff_transformation_3< K >)

#  define CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_2(K) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_3(K)

#if defined(LEDA_NAMESPACE)
namespace leda {
class real;
class integer;
class rational;
class bigfloat;
}
#else
class leda_real;
class leda_integer;
class leda_rational;
class leda_bigfloat;
#endif

namespace CGAL {
class Gmpz;
template <class NumberType> class Quotient;
}

// NB : we could add MP_Float, and maybe other useful NTs.

#  define CGAL_ITERATOR_TRAITS_POINTER_SPEC_TEMPLATE(K) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <int>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <long>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <float>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <double>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <leda_real>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <leda_integer>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <leda_rational>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <leda_bigfloat>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <CGAL::Gmpz>) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <CGAL::Quotient<int> >) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <CGAL::Quotient<long> >) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <CGAL::Quotient<float> >) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <CGAL::Quotient<double> >) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <CGAL::Quotient<leda_real> >) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <CGAL::Quotient<leda_integer> >) \
     CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K <CGAL::Quotient<CGAL::Gmpz> >)

#else
#  define CGAL_ITERATOR_TRAITS_POINTER_SPEC_(K)
#  define CGAL_ITERATOR_TRAITS_POINTER_SPEC_TEMPLATE(K)
#endif

#endif  // CGAL_USER_CLASSES_H
