// Copyright (c) 2007-2010 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>



#ifndef CGAL_MPFI_COERCION_TRAITS_H
#define CGAL_MPFI_COERCION_TRAITS_H

#include <CGAL/config.h>

#ifdef CGAL_USE_MPFI

#include <CGAL/number_type_basic.h>
#include <CGAL/GMP/Gmpfr_type.h>
#include <CGAL/GMP/Gmpfi_type.h>
#include <CGAL/Coercion_traits.h>

namespace CGAL{

CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Gmpfi)

// coercion with native types
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long        ,Gmpfi)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(unsigned long       ,Gmpfi)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int         ,Gmpfi)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short       ,Gmpfi)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double      ,Gmpfi)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float       ,Gmpfi)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double ,Gmpfi)

// coercion with gmp and mpfr types
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpz        ,Gmpfi)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpq        ,Gmpfi)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpfr       ,Gmpfi)

}

#endif        // CGAL_USE_MPFI
#endif  // CGAL_MPFI_COERCION_TRAITS_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
