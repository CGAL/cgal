// Copyright (c) 2001
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_PRECISE_NUMBERS_H
#define CGAL_PRECISE_NUMBERS_H

#include <CGAL/config.h>
#if defined CGAL_USE_GMPXX
#  include <CGAL/gmpxx.h>
typedef mpz_class                       Precise_integer;
typedef mpq_class                       Precise_rational;
#elif defined CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
#  include <CGAL/leda_rational.h>
typedef leda_integer                    Precise_integer;
typedef leda_rational                   Precise_rational;
#elif defined CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#  include <CGAL/Gmpq.h>
typedef CGAL::Gmpz                      Precise_integer;
typedef CGAL::Gmpq                      Precise_rational;
#else
#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
typedef CGAL::MP_Float                  Precise_integer;
typedef CGAL::Quotient<CGAL::MP_Float>  Precise_rational;
#endif

#endif // CGAL_PRECISE_NUMBERS_H
