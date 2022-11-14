// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


/*! \file NiX/Gmp/Coercion_traits.h
 *  \brief provides specializations of Coercion_traits for the Gmp types.
 */

#ifndef CGAL_GMP_COERCION_TRAITS_H
#define CGAL_GMP_COERCION_TRAITS_H 1

#include <CGAL/number_type_basic.h>
#include <CGAL/GMP/Gmpz_type.h>
#include <CGAL/GMP/Gmpzf_type.h>
#include <CGAL/GMP/Gmpfr_type.h>
#include <CGAL/GMP/Gmpq_type.h>
#include <CGAL/Coercion_traits.h>

#ifdef CGAL_USE_GMP

namespace CGAL {

//Gmp internal coercions:
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Gmpz)
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Gmpzf)
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Gmpq)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpz,Gmpzf)
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpzf,Gmpzq); // todo
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpz,Gmpq)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpfr,Gmpq)

// The following definitions reflect the interaction of the Gmp number types
// with the built in types,
// Gmpz:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,Gmpz)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,Gmpz)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,Gmpz)

// Gmpzf: not tested
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,Gmpzf)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,Gmpzf)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,Gmpzf)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float    ,Gmpzf)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double   ,Gmpzf)

// Gmpq:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,Gmpq)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,Gmpq)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,Gmpq)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,Gmpq)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,Gmpq)

} //namespace CGAL


#endif // CGAL_USE_GMP
#endif //CGAL_GMP_COERCION_TRAITS_H 1
//EOF
