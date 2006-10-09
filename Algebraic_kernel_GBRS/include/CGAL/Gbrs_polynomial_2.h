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

#ifndef CGAL_GBRS_POLYNOMIAL_2_H
#define CGAL_GBRS_POLYNOMIAL_2_H

#include <iostream>
#include <CGAL/assertions.h>
#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/MpfiInterval.h>
#include <gmp.h>

CGAL_BEGIN_NAMESPACE

// TODO: change the implementation so the type of coefficients are a parameter
// to the class (it won't be so efficient, though)

class Rational_polynomial_2 {
	// The coefficients are stored in form of a matrix, for example for the
	// polynomial of degree 3 in x and 2 in y, P(x) = 4 * x^3 * y^2 + 7 * x
	//
	// this matrix is **coef:
	//      y^0 y^1 y^2
	//     +---+---+---+
	// x^0 | 0 | 0 | 0 | <- this line is coef[0]
	//     +---+---+---+
	// x^1 | 7 | 0 | 0 | <- this line is coef[1]
	//     +---+---+---+
	// x^2 | 0 | 0 | 0 | <- this line is coef[2]
	//     +---+---+---+
	// x^3 | 0 | 0 | 4 | <- this line is coef[3]
	//     +---+---+---+
	// Clearly, the matrix has dimensions (deg_x+1) * (deg_y+1)
	private:
		int degree_x, degree_y;	// degrees in both variables
		mpq_t** coef;	// the coefficients
	public:
		// construction and destruction
		Rational_polynomial_2 ();
		Rational_polynomial_2 (unsigned int, unsigned int);
		Rational_polynomial_2 (int, int);
		Rational_polynomial_2 (const Rational_polynomial_2 &);
		~Rational_polynomial_2 ();
		// other functions
		inline int get_degree_x () const;
		inline int get_degree_y () const;
		inline mpq_t** get_coefs () const;
		void get_coef (int, int, mpq_t *) const;
		void set_coef (int, int, const mpq_t &);
		void set_coef (int, int, const mpz_t &);
		void set_coef (int, int, const CGAL::Gmpq &);
		void set_coef (int, int, const CGAL::Gmpz &);
		void set_coef (int, int, int);
		void set_coef (int, int, unsigned int);
		void set_coef (int, int, long int);
		void set_coef (int, int, unsigned long int);
		Rational_polynomial_2 operator= (const Rational_polynomial_2 &);
		Rational_polynomial_2 operator- () const;
		Rational_polynomial_2 operator+ (const Rational_polynomial_2 &) const;
		Rational_polynomial_2& operator+= (const Rational_polynomial_2 &);
		Rational_polynomial_2 operator- (const Rational_polynomial_2 &) const;
		Rational_polynomial_2& operator-= (const Rational_polynomial_2 &);
		void scale (const int);
		void scale (const mpz_t &);
		void scale (const mpq_t &);
		void scale (const CGAL::Gmpz &);
		void scale (const CGAL::Gmpq &);
		void shift (int, int);
		Rational_polynomial_2 operator* (const Rational_polynomial_2 &) const;
		Rational_polynomial_2& operator*= (const Rational_polynomial_2 &);
		bool operator== (const Rational_polynomial_2 &) const;
		inline bool operator!= (const Rational_polynomial_2 &) const;
		std::ostream& show (std::ostream &) const;
		template <class T> Rational_polynomial_2 operator* (const T &) const;
};

std::ostream& operator<< (std::ostream &, const Rational_polynomial_2 &);
template <class T> inline Rational_polynomial_2 operator* (const T &,
		Rational_polynomial_2 &);

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Gbrs_polynomial_2.C>
#endif	// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif	// CGAL_GRBS_POLYNOMIAL_2_H
