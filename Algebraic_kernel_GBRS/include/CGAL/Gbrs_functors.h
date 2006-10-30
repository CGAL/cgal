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

#ifndef GBRS_FUNCTORS_H
#define GBRS_FUNCTORS_H

#include <mpfi.h>
#include <CGAL/enum.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <CGAL/Gbrs_solve_1.h>

CGAL_BEGIN_NAMESPACE

namespace AlgebraicFunctors {

template <class AK>
class Construct_polynomial_1 {
	typedef typename AK::Polynomial_1	Polynomial_1;
	public:
	template <class InputIterator>
	Polynomial_1 operator()
		(InputIterator first, InputIterator last) const {
		// count the number of elements in the container
		int elements = 0;
		InputIterator it;
		for (it = first; it != last; ++it)
			++elements;
		CGAL_assertion_msg (elements, "the container can't be empty");
		// TODO: if the last elements in the container are zero...
		int degree = elements-1;
		Polynomial_1 p (degree);
		for (it = first; it != last; ++it)
			p.set_coef (degree-(--elements), *it);
		return p;
	};

	// TODO: check this "functor" (what is a functor in C++?)

	// great! now I know what is that... now I hate C++ a lot more than
	// ten minutes ago... f_-_!
	template <class InputIterator1, class InputIterator2>
	Polynomial_1 operator()
		(InputIterator1 first_coeff, InputIterator1 last_coeff,
		InputIterator2 deg) {
		// the degree of the polynomial will be the greater degree
		unsigned int greater = 0;
		InputIterator1 c;
		InputIterator2 d = deg;
		for (c = first_coeff; c != last_coeff; ++c)
			if (d > greater)
				greater = d++;
		// now, construct the polynomial of degree d
		Polynomial_1 p (d);
		for (c = first_coeff; c != last_coeff; ++c)
			p.set_coef (deg++, *c);
		return p;
	};
};	// Construct_polynomial_1

template <class AK>
class Solve_1 {
	typedef typename AK::Algebraic_real_1	Algebraic;
	typedef typename AK::Polynomial_1	Polynomial_1;
	public:
	template <class OutputIterator>
	OutputIterator operator() (const Polynomial_1 &p,
			OutputIterator res,
			bool known_to_be_square_free) const {
		if (known_to_be_square_free)
			return res;
		mpfi_ptr *x;
		int nr;
		CGAL_assertion_msg (((nr = solve_1 (x, p)) >= 0),
				"error in resolution");
		if (nr)
			for (int i=0; i<nr; ++i) {
				// multiplicity is -1 (we didn't calculate it)
				Algebraic a (x[i], p, i, -1, CGAL_RS_DEF_PREC);
				*(res++) = a;
			}
		free (x);
		return res;
	};

	// TODO: how the hell does RS compute the multiplicity???
	template <class OutputIteratorRoots, class OutputIteratorMult>
	std::pair<OutputIteratorRoots, OutputIteratorMult>
	operator() (const Polynomial_1 &p,
			OutputIteratorRoots roots,
			OutputIteratorMult mult) const {
		CGAL_assertion_msg (false, "not implemented yet");
		return std::make_pair (roots, mult);
	};
};	// Solve_1

template <class AK>
class SignAt_1 {
	typedef typename AK::Polynomial_1	Polynomial_1;
	typedef typename AK::Algebraic_real_1	Algebraic_1;
	typedef typename AK::Coefficient	Coefficient_1;
	public:
	typedef Sign	result_type;
	inline result_type operator() (const Polynomial_1 &p,
			const Algebraic_1 &r) const {
		return sign_1 (p, r);
	};

	result_type operator() (const Polynomial_1 &p,
			const Coefficient_1 &r) const {
		Coefficient_1 c = p.eval (r);
		return ((c<0)?NEGATIVE:((c>0)?POSITIVE:ZERO));
	};
};	// SignAt_1

template <class AK>
class Derivative_1 {
	typedef typename AK::Polynomial_1	Polynomial_1;
	public:
	inline Polynomial_1 operator() (const Polynomial_1 &p) const {
		return p.derive ();
	};
};	// Derivative_1

template <class AK>
class Compare_1 {
	typedef typename AK::Algebraic_real_1	Algebraic_1;
	public:
	inline Comparison_result operator()
		(Algebraic_1 &r1, Algebraic_1 &r2) const {
			return compare_1 (r1, r2);
		};
};	// Compare_1

}	// namespace AlgebraicFunctors

CGAL_END_NAMESPACE

#endif	// GBRS_FUNCTORS_H
