// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_H

#include <CGAL/basic.h>

namespace CGAL {
  namespace AlgebraicFunctors {

  template < class AK, class OutputIterator >
  inline
  OutputIterator
  solve( const typename AK::Polynomial_for_circles_2_2 & e1,
	 const typename AK::Polynomial_for_circles_2_2 & e2,
	 OutputIterator res )
  {
    CGAL_precondition( ! (e1 == e2) ); // polynomials of this type cannot be multiple
    // of one another if they are not equal

    typedef typename AK::FT FT;
    typedef typename AK::Root_of_2 Root_of_2;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
    const FT dx = e2.a() - e1.a();
    const FT dy = e2.b() - e1.b();

    const FT dx2 = CGAL::square(dx);
    const FT dy2 = CGAL::square(dy);
    const FT dist2 = dx2 + dy2; // squared distance between centers
    const FT diff_sqr_rad = e1.r_sq() - e2.r_sq();
    const FT disc = 2*dist2*(e1.r_sq() + e2.r_sq()) -
                      (CGAL::square(diff_sqr_rad) + CGAL::square(dist2));
    CGAL::Sign sign_disc = CGAL::sign(disc);

    if (sign_disc == NEGATIVE) return res;

    const FT x_base = ((e1.a() + e2.a()) + dx*diff_sqr_rad / dist2) / 2;
    const FT y_base = ((e1.b() + e2.b()) + dy*diff_sqr_rad / dist2) / 2;

    if (sign_disc == ZERO) {
      // one double root,
      // no need to care about the boolean of the Root_of
      *res++ = std::make_pair
	( Root_for_circles_2_2
	  (Root_of_2(x_base), Root_of_2(y_base)),
	  static_cast<unsigned>(2) ); // multiplicity = 2
      return res;
    }

    CGAL::Sign sign_dy = CGAL::sign (dy);
    CGAL::Sign sign_dx = CGAL::sign (dx);

    // else, 2 distinct roots
    if (sign_dy == ZERO) {
      const FT y_root_coeff = dx / (2 * dist2);
      if(sign_dx == NEGATIVE) {
        * res++ = std::make_pair
	( Root_for_circles_2_2(Root_of_2(x_base),
                               make_root_of_2(y_base, y_root_coeff, disc)),
	  static_cast<unsigned>(1) );

        * res++ = std::make_pair
	( Root_for_circles_2_2(Root_of_2(x_base),
                               make_root_of_2(y_base, -y_root_coeff, disc)),
	  static_cast<unsigned>(1) );
      } else {
        * res++ = std::make_pair
	( Root_for_circles_2_2(Root_of_2(x_base),
                               make_root_of_2(y_base, -y_root_coeff, disc)),
	  static_cast<unsigned>(1) );

        * res++ = std::make_pair
	( Root_for_circles_2_2(Root_of_2(x_base),
                               make_root_of_2(y_base, y_root_coeff, disc)),
	  static_cast<unsigned>(1) );
      }
      return res;
    }

    if (sign_dx == ZERO) {
      const FT x_root_coeff = dy / (2 * dist2);
      if(sign_dy == POSITIVE) {
        * res++ = std::make_pair
	( Root_for_circles_2_2(make_root_of_2(x_base, -x_root_coeff, disc),
                               Root_of_2(y_base)),
	  static_cast<unsigned>(1) );

        * res++ = std::make_pair
	( Root_for_circles_2_2(make_root_of_2(x_base, x_root_coeff, disc),
                               Root_of_2(y_base)),
	  static_cast<unsigned>(1) );
      } else {
        * res++ = std::make_pair
	( Root_for_circles_2_2(make_root_of_2(x_base, x_root_coeff, disc),
                               Root_of_2(y_base)),
	  static_cast<unsigned>(1) );

        * res++ = std::make_pair
	( Root_for_circles_2_2(make_root_of_2(x_base, -x_root_coeff, disc),
                               Root_of_2(y_base)),
	  static_cast<unsigned>(1) );
      }
      return res;
    }

    const FT   x_root_coeff = dy / (2 * dist2);
    const FT   y_root_coeff = dx / (2 * dist2);

    if (sign_dy == POSITIVE) {
      * res++ = std::make_pair
      ( Root_for_circles_2_2(make_root_of_2(x_base, -x_root_coeff, disc),
                             make_root_of_2(y_base, y_root_coeff, disc)),
	static_cast<unsigned>(1) );

      * res++ = std::make_pair
      ( Root_for_circles_2_2(make_root_of_2(x_base, x_root_coeff, disc),
                             make_root_of_2(y_base, -y_root_coeff, disc)),
	static_cast<unsigned>(1) );
    } else {
      * res++ = std::make_pair
      ( Root_for_circles_2_2(make_root_of_2(x_base, x_root_coeff, disc),
                             make_root_of_2(y_base, -y_root_coeff, disc)),
	static_cast<unsigned>(1) );
      * res++ = std::make_pair
      ( Root_for_circles_2_2(make_root_of_2(x_base, -x_root_coeff, disc),
                             make_root_of_2(y_base, y_root_coeff, disc)),
	static_cast<unsigned>(1) );
    }

    return res;
  }

  template < class AK >
  inline
  Sign sign_at( const typename AK::Polynomial_for_circles_2_2 & equation,
	 const typename AK::Root_for_circles_2_2 & r)
  {
    Comparison_result c = compare(square(r.x() - equation.a()),
                                  equation.r_sq() -
                                    square(r.y() - equation.b()));
    if(c == EQUAL) return ZERO;
    if(c == LARGER) return POSITIVE;
    return NEGATIVE;
  }


  template <class AK>
  typename AK::Root_for_circles_2_2
  x_critical_point(const typename AK::Polynomial_for_circles_2_2 & c,
		   bool i)
  {
    typedef typename AK::Root_of_2            Root_of_2;
    typedef typename AK::FT                   FT;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;

    const Root_of_2 a1 = make_root_of_2(c.a(),FT(i?-1:1),c.r_sq());
    return Root_for_circles_2_2(a1, c.b());
  }

  template <class AK, class OutputIterator>
  OutputIterator
  x_critical_points(const typename AK::Polynomial_for_circles_2_2 & c,
		    OutputIterator res)
  {
    typedef typename AK::FT                   FT;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;

    *res++ =  Root_for_circles_2_2(
            make_root_of_2(c.a(),FT(-1),c.r_sq()), c.b());
    *res++ =  Root_for_circles_2_2(
            make_root_of_2(c.a(),FT(1),c.r_sq()), c.b());

    return res;
  }

  template <class AK>
  typename AK::Root_for_circles_2_2
  y_critical_point(const typename AK::Polynomial_for_circles_2_2 &c,
		   bool i)
  {
    typedef typename AK::Root_of_2            Root_of_2;
    typedef typename AK::FT                   FT;
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;

    const Root_of_2 b1 = make_root_of_2(c.b(),FT(i?-1:1),c.r_sq());
    return Root_for_circles_2_2(c.a(),b1);
  }

  template <class AK, class OutputIterator>
  OutputIterator
  y_critical_points(const typename AK::Polynomial_for_circles_2_2 & c,
		    OutputIterator res)
  {
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;

    *res++ = Root_for_circles_2_2(c.a(),
      make_root_of_2(c.b(),-1,c.r_sq()));
    *res++ = Root_for_circles_2_2(c.a(),
      make_root_of_2(c.b(),1,c.r_sq()));

    return res;
  }



  } // namespace AlgebraicFunctors
} // namespace CGAL

#endif //  CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_H
