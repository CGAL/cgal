// Copyright (c) 2007-2010 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_MPFR_COERCION_TRAITS_H
#define CGAL_MPFR_COERCION_TRAITS_H

#include <CGAL/config.h>

#ifdef CGAL_USE_MPFR

#include <CGAL/number_type_basic.h>
#include <CGAL/GMP/Gmpfr_type.h>
#include <CGAL/Coercion_traits.h>

namespace CGAL{

CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Gmpfr)

// coercion with native types
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long        ,Gmpfr)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(unsigned long       ,Gmpfr)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int         ,Gmpfr)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short       ,Gmpfr)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double      ,Gmpfr)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float       ,Gmpfr)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double ,Gmpfr)

// coercion with gmp types
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpz        ,Gmpfr)

}

#endif  // CGAL_USE_MPFR
#endif  // CGAL_MPFR_COERCION_TRAITS_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
