// Copyright (c) 2006-2009, 2011 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     :  Eric Berberich <eric@mpi-inf.mpg.de>
//                  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//                  Luis Peñaranda <luis.penaranda@gmx.com>
//
// ============================================================================

/*! \file RS/isolator.h
  \brief Defines class CGAL::RS_real_root_isolator

  Isolate real roots of polynomials with Fabrice Roullier's Rs.

  The polynomial has to be a univariate polynomial over any number type
  which is contained in the real numbers. The polynomial does not need to
  be square-free.

*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_RS_ISOLATOR_H
#define CGAL_ALGEBRAIC_KERNEL_D_RS_ISOLATOR_H

#include <CGAL/config.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>
#include <CGAL/RS/isole_1.h>

namespace CGAL {

namespace internal {

/*! \brief A model of concept RealRootIsolator.

   Polynomial_ must be Polynomial<Gmpz>, and Bound_ must be Gmpfr.

 */
template <class Polynomial_, class Bound_>
class RS_real_root_isolator {

public:
    //! First template parameter
    typedef Polynomial_ Polynomial;

    //! Second template parameter
    typedef Bound_ Bound;

private:

    //! Coefficient type of polynomial
    typedef typename Polynomial::NT Coefficient;

    typedef Gmpfi Interval;

public:
    /*! \brief Constructor from univariate square free polynomial.

    The RealRootIsolator provides isolating intervals for the real
    roots of the polynomial
    */
   RS_real_root_isolator(const Polynomial& p = Polynomial(Coefficient(0))) :
      _m_polynomial(p)
      //, _m_interval_given(false)
    {
      _m_real_roots=RS::isolator<Polynomial>()(p);
    }

    // @LUIS: add constructor from interval (maybe some time in the future)

public: // functions

    /*! \brief returns the defining polynomial*/
    Polynomial polynomial() const {
      return _m_polynomial;
    }

    //! returns the number of real roots
    int number_of_real_roots() const {
      return _m_real_roots.size();
    }

    /*! \brief returns true if the isolating interval is degenerated to a
      single point.

      If is_exact_root(i) is true,
      then left_bound(int i) equals  \f$root_i\f$. \n
      If is_exact_root(i) is true,
      then right_bound(int i) equals  \f$root_i\f$. \n
    */
    bool is_exact_root(int i) const {
      return(_m_real_roots[i].inf()==_m_real_roots[i].sup());
    }

public:

    /*! \brief returns  \f${l_i}\f$ the left bound of the isolating interval
      for root  \f$root_{i}\f$.

      In case is_exact_root(i) is true,  \f$l_i = root_{i}\f$,\n
      otherwise:  \f$l_i < root_{i}\f$.

      If  \f$i-1>=0\f$, then  \f$l_i > root_{i-1}\f$. \n
      If  \f$i-1>=0\f$, then  \f$l_i >= r_{i-1}\f$,
      the right bound of  \f$root_{i-1}\f$\n

      \pre 0 <= i < number_of_real_roots()
    */
    Bound left_bound(int i) const {
      CGAL_assertion(i >= 0);
      CGAL_assertion(i < this->number_of_real_roots());
      return _m_real_roots[i].inf();
    }

    /*! \brief returns  \f${r_i}\f$ the right bound of the isolating interval
      for root  \f$root_{i}\f$.

      In case is_exact_root(i) is true,  \f$r_i = root_{i}\f$,\n
      otherwise:  \f$r_i > root_{i}\f$.

      If  \f$i+1< n \f$, then  \f$r_i < root_{i+1}\f$,
      where \f$n\f$ is number of real roots.\n
      If  \f$i+1< n \f$, then  \f$r_i <= l_{i+1}\f$,
      the left bound of  \f$root_{i+1}\f$\n

      \pre 0 <= i < number_of_real_roots()
    */
    Bound right_bound(int i) const {
      CGAL_assertion(i >= 0);
      CGAL_assertion(i < this->number_of_real_roots());
      return _m_real_roots[i].sup();
    }

private:

    //! the input polynomial
    Polynomial _m_polynomial;

    //! the solutions
    std::vector<Interval> _m_real_roots;

    //! restricted interval?
    // TODO bool _m_interval_given;

};

} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_RS_ISOLATOR_H
