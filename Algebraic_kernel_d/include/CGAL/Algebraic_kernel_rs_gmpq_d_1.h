// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_ALGEBRAIC_KERNEL_RS_GMPQ_D_1
#define CGAL_ALGEBRAIC_KERNEL_RS_GMPQ_D_1

#include <CGAL/disable_warnings.h>

#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Polynomial.h>
#include "RS/rs2_isolator_1.h"
#ifdef CGAL_USE_RS3
#include "RS/rs23_k_isolator_1.h"
#include "RS/rs3_refiner_1.h"
#include "RS/rs3_k_refiner_1.h"
#else
#include "RS/bisection_refiner_1.h"
#endif
#ifdef CGAL_RS_USE_NATIVE_GMPQ_KERNEL
#include "RS/ak_1.h"
#else
#include "RS/ak_z_1.h"
#endif

// The RS algebraic kernel for non-Gmpz types comes in two flavors. The
// native kernel converts, before each operation, the input polynomial to a
// Polynomial<Gmpz>. The z-kernel only converts to a Polynomial<Gmpz>
// before isolation and stores both polynomials in the algebraic number.
// The latter is the default behavior, but the former can be selected by
// setting the compilation flag CGAL_RS_USE_NATIVE_GMPQ_KERNEL.

namespace CGAL{

#ifdef CGAL_RS_USE_NATIVE_GMPQ_KERNEL // Use the native kernel.
// Choice of the isolator: RS default or RS-k.
#ifdef CGAL_RS_USE_K
typedef CGAL::RS23_k_isolator_1<CGAL::Polynomial<CGAL::Gmpq>,CGAL::Gmpfr>
                                        QIsolator;
#else
typedef CGAL::RS2::RS2_isolator_1<CGAL::Polynomial<CGAL::Gmpq>,CGAL::Gmpfr>
                                        QIsolator;
#endif

// Choice of the refiner: bisection, RS3 or RS3-k.
#ifdef CGAL_USE_RS3
#ifdef CGAL_RS_USE_K
typedef CGAL::RS3::RS3_k_refiner_1<CGAL::Polynomial<CGAL::Gmpq>,CGAL::Gmpfr>
                                        QRefiner;
#else
typedef CGAL::RS3::RS3_refiner_1<CGAL::Polynomial<CGAL::Gmpq>,CGAL::Gmpfr>
                                        QRefiner;
#endif
#else
typedef CGAL::Bisection_refiner_1<CGAL::Polynomial<CGAL::Gmpq>,CGAL::Gmpfr>
                                        QRefiner;
#endif

typedef CGAL::RS_AK1::Algebraic_kernel_1<
        CGAL::Polynomial<CGAL::Gmpq>,
        CGAL::Gmpfr,
        QIsolator,
        QRefiner>
                                        Algebraic_kernel_rs_gmpq_d_1;

#else // Use the z-kernel.

// Choice of the z-isolator: RS default or RS-k.
#ifdef CGAL_RS_USE_K
typedef CGAL::RS23_k_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZIsolator;
#else
typedef CGAL::RS2::RS2_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZIsolator;
#endif

// Choice of the z-refiner: bisection, RS3 or RS3-k.
#ifdef CGAL_USE_RS3
#ifdef CGAL_RS_USE_K
typedef CGAL::RS3::RS3_k_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZRefiner;
#else
typedef CGAL::RS3::RS3_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZRefiner;
#endif
#else
typedef CGAL::Bisection_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZRefiner;
#endif

typedef CGAL::RS_AK1::Algebraic_kernel_z_1<
        CGAL::Polynomial<CGAL::Gmpq>,
        CGAL::Polynomial<CGAL::Gmpz>,
        CGAL::RS_AK1::Polynomial_converter_1<CGAL::Polynomial<CGAL::Gmpq>,
                                             CGAL::Polynomial<CGAL::Gmpz> >,
        CGAL::Gmpfr,
        ZIsolator,
        ZRefiner>
                                                Algebraic_kernel_rs_gmpq_d_1;

#endif

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif  // CGAL_ALGEBRAIC_KERNEL_RS_GMPQ_D_1
