// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.

// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_ALGEBRAIC_KERNEL_RS_GMPZ_D_1
#define CGAL_ALGEBRAIC_KERNEL_RS_GMPZ_D_1

#include <CGAL/disable_warnings.h>

#include <CGAL/Gmpz.h>
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
#include "RS/ak_1.h"

namespace CGAL{

// Choice of the isolator: RS default or RS-k.
#ifdef CGAL_RS_USE_K
typedef CGAL::RS23_k_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                Isolator;
#else
typedef CGAL::RS2::RS2_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                Isolator;
#endif

// Choice of the refiner: bisection, RS3 or RS3-k.
#ifdef CGAL_USE_RS3
#ifdef CGAL_RS_USE_K
typedef CGAL::RS3::RS3_k_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                Refiner;
#else
typedef CGAL::RS3::RS3_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                Refiner;
#endif
#else
typedef CGAL::Bisection_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                Refiner;
#endif

typedef CGAL::RS_AK1::Algebraic_kernel_1<
        CGAL::Polynomial<CGAL::Gmpz>,
        CGAL::Gmpfr,
        Isolator,
        Refiner>
                                                Algebraic_kernel_rs_gmpz_d_1;

}

#include <CGAL/enable_warnings.h>

#endif // CGAL_ALGEBRAIC_KERNEL_RS_GMPZ_D_1
