// Copyright (c) 2011 Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

#ifndef CGAL_BASE_RATIONAL_ARC_DS_D_1_H
#define CGAL_BASE_RATIONAL_ARC_DS_D_1_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <vector>
#include <ostream>
#include <CGAL/Arr_enums.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>

#include <CGAL/Fraction_traits.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>

#include <boost/type_traits/is_same.hpp>

namespace CGAL {
namespace Arr_rational_arc {

template <typename Algebraic_kernel_ >
class Base_rational_arc_ds_1
{
public:
  typedef Algebraic_kernel_                             Algebraic_kernel;
  typedef Base_rational_arc_ds_1<Algebraic_kernel>      Self;

  //typedef typename Algebraic_kernel::Multiplicity_type  Multiplicity;
  typedef unsigned int                                  Multiplicity;
  typedef typename Algebraic_kernel::Coefficient        Coefficient;

  typedef typename Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel
                                                        Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Rational          Rational;
  typedef typename Arithmetic_kernel::Integer           Integer;
  typedef typename Algebraic_kernel::Algebraic_real_1   Algebraic_real_1;

  typedef typename Algebraic_kernel::Polynomial_1       Polynomial_1;
  typedef Polynomial_traits_d<Polynomial_1>             Polynomial_traits_1;
  typedef Fraction_traits<Rational>                     FT_rat_1;
  typedef typename Algebraic_kernel::Solve_1            Solve_1;
  typedef typename Algebraic_kernel::Bound              Bound;
  typedef Algebraic_structure_traits<Polynomial_1>      AT_poly;

  typedef Polynomial<Rational>                          Poly_rat_1;
  typedef Polynomial_traits_d<Poly_rat_1>               PT_rat_1;
  typedef Fraction_traits <Poly_rat_1>                  FT_poly_rat_1;
  typedef std::vector<Algebraic_real_1>                 Algebraic_vector;
  typedef std::vector<Multiplicity>                     Multiplicity_vector;
  typedef std::vector<std::pair<Algebraic_real_1, Multiplicity> >
                                                        Root_multiplicity_vector;

  CGAL_static_assertion((boost::is_same<Integer,Coefficient>::value));
  CGAL_static_assertion((boost::is_same<Polynomial_1,
                       typename FT_poly_rat_1::Numerator_type>::value));

public:

  //---------------------------------------------------------------------
  // Print a polynomial nicely.

  static std::ostream& print_polynomial(std::ostream& os,
                                        const Polynomial_1& poly,
                                        char var)
  {
    // Get the degree.
    const int    deg = CGAL::degree(poly);

    Integer     coeff;
    CGAL::Sign  sgn;
    int         k;

    if (deg < 0)
    {
      os << '0';
      return (os);
    }

    for (k = deg; k >= 0; k--)
    {
      //coeff = pt::Get_coefficient()(poly, k);
      coeff = CGAL::get_coefficient(poly, k);

      if (k == deg)
        os << coeff;
      else if ((sgn = CGAL::sign (coeff)) == POSITIVE)
        os << " + " << coeff;
      else if (sgn == NEGATIVE)
        os << " - " << -coeff;
      else
        continue;

      if (k > 1)
        os << '*' << var << '^' << k;
      else if (k == 1)
        os << '*' << var;
    }

    return (os);
  }

}; //Base_rational_arc_ds_1

}   // namespace Arr_rational_arc
}   //namespace CGAL {

#endif //CGAL_BASE_RATIONAL_ARC_DS_D_1_H
