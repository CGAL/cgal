// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_PRECISE_NUMBERS_H
#define CGAL_PRECISE_NUMBERS_H

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
