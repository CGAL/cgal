// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri
//                 Stefan Schirra

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
