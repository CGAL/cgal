// Copyright (c) 2004  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

#ifdef CGAL_USE_GMP

#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpz.h>

CGAL_BEGIN_NAMESPACE

double to_double(const Quotient<Gmpz>& quot)
{
  mpq_t  mpQ;
  mpq_init(mpQ);
  const Gmpz& n = quot.numerator();
  const Gmpz& d = quot.denominator();
  mpz_set(mpq_numref(mpQ), n.mpz());
  mpz_set(mpq_denref(mpQ), d.mpz());
    
  mpq_canonicalize(mpQ);
  
  double ret = mpq_get_d(mpQ);
  mpq_clear(mpQ);
  return ret;
}

CGAL_END_NAMESPACE

#endif // CGAL_USE_GMP
