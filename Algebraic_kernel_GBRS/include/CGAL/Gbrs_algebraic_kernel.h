// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
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
// $URL$
// $Id$
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#ifndef CGAL_GBRS_ALGEBRAIC_KERNEL
#define CGAL_GBRS_ALGEBRAIC_KERNEL

#include <CGAL/Algebraic_1.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <CGAL/Gbrs_functors.h>

CGAL_BEGIN_NAMESPACE

template <class IntegralDomain_>
class GBRS_algebraic_kernel {
	typedef GBRS_algebraic_kernel<IntegralDomain_>	Self;

	public:

	// constructor: we must initialize RS just a time, so this is a good
	// time to do it
	GBRS_algebraic_kernel () {
		CGAL_assertion_msg (!(init_rs ()), "error initializing RS");
	};

	typedef IntegralDomain_				Coefficient;
	typedef Rational_polynomial_1			Polynomial_1;
	typedef Algebraic_1				Algebraic_real_1;
	typedef AlgebraicFunctors::Construct_polynomial_1<Self>
							Construct_polynomial_1;
	typedef AlgebraicFunctors::Solve_1<Self>	Solve_1;
	typedef AlgebraicFunctors::SignAt_1<Self>	SignAt_1;
	typedef AlgebraicFunctors::Derivative_1<Self>	Derivative_1;
	typedef AlgebraicFunctors::Compare_1<Self>	Compare_1;

	Construct_polynomial_1 construct_polynomial_1_object () const {
		return Construct_polynomial_1 ();
	}

	Solve_1 construct_solve_1_object () const {
		return Solve_1 ();
	}

	SignAt_1 construct_signat_1_object () const {
		return SignAt_1 ();
	}

	Derivative_1 construct_derivative_1_object () const {
		return Derivative_1 ();
	}

	Compare_1 construct_compare_1_object () const {
		return Compare_1 ();
	}

};	// GBRS_algebraic_kernel

CGAL_END_NAMESPACE

#endif	// CGAL_GBRS_ALGEBRAIC_KERNEL
