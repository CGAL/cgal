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

#include <CGAL/assertions.h>
#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <CGAL/MpfiInterval.h>
#include <iostream>
#include <gmp.h>

CGAL_BEGIN_NAMESPACE

// constructors
Rational_polynomial_1::Rational_polynomial_1 () {
	degree = 0;
	nm = 1;
	coef = (mpq_t*)malloc (sizeof(mpq_t));
	mpq_init (coef[0]);
};

Rational_polynomial_1::Rational_polynomial_1 (unsigned int d) {
	degree = (int)d;
	nm = degree+1;
	coef = (mpq_t*)malloc (sizeof(mpq_t)*nm);
	for (int i=0; i<nm; ++i)
		mpq_init (coef[i]);
};

Rational_polynomial_1::Rational_polynomial_1 (int d) {
	degree = d<0?0:d;
	nm = degree+1;
	coef = (mpq_t*)malloc (sizeof(mpq_t)*nm);
	for (int i=0; i<nm; ++i)
		mpq_init (coef[i]);
};

Rational_polynomial_1::Rational_polynomial_1 (const Rational_polynomial_1 &p) {
	degree = p.degree;
	nm = p.nm;
	mpq_t *p_coef = p.get_coefs ();
	coef = (mpq_t*)malloc (sizeof(mpq_t)*nm);
	// we have to copy the contents, not just the pointer
	for (int i=0; i<nm; ++i) {
		mpq_init (coef[i]);
		mpq_set (coef[i], p_coef[i]);
		}
};

// destructor
Rational_polynomial_1::~Rational_polynomial_1 () {
	for (int i=0; i<nm; ++i)
		mpq_clear (coef[i]);
	free (coef);
};

void Rational_polynomial_1::set_degree (int d) {	// noone should call this function
	for (int i=0; i<nm; ++i)	// free the old coefficients
		mpq_clear (coef[i]);
	free (coef);
	degree = d;	// do it again
	nm = d+1;
	coef = (mpq_t*)malloc (sizeof(mpq_t)*nm);
	for (int i=0; i<nm; ++i)
		mpq_init (coef[i]);
	return;
};

int Rational_polynomial_1::calc_index (int pow_x) const {
	// we have to calculate the index in the array to store the coefficient
	// of x^pow_x
	return degree-pow_x;
};

// functions for setting polynomial coefficients
void Rational_polynomial_1::set_coef (int pow_x, const mpq_t &q) {
	mpq_set (coef[calc_index (pow_x)], q);
};

void Rational_polynomial_1::set_coef (int pow_x, const mpz_t &z) {
	mpq_set_z (coef[calc_index (pow_x)], z);
};

void Rational_polynomial_1::set_coef (int pow_x, const CGAL::Gmpq &q) {
	mpq_set (coef[calc_index (pow_x)], q.mpq());
};

void Rational_polynomial_1::set_coef (int pow_x, const CGAL::Gmpz &z) {
	mpq_set_z (coef[calc_index (pow_x)], z.mpz());
};

void Rational_polynomial_1::set_coef (int pow_x, int c) {
	mpq_set_si (coef[calc_index (pow_x)], (long int)c, (unsigned long int)1);
};

void Rational_polynomial_1::set_coef (int pow_x, unsigned int c) {
	mpq_set_ui (coef[calc_index (pow_x)], (unsigned long int)c, (unsigned long int)1);
};

void Rational_polynomial_1::set_coef (int pow_x, long int c) {
	mpq_set_si (coef[calc_index (pow_x)], c, (unsigned long int)1);
};

void Rational_polynomial_1::set_coef (int pow_x, unsigned long int c) {
	mpq_set_ui (coef[calc_index (pow_x)], c, (unsigned long int)1);
};

void Rational_polynomial_1::set_coef_rat_li (int pow_x, long int num, long int den) {
	mpq_set_si (coef[calc_index (pow_x)], num, den);
};

// scaling
void Rational_polynomial_1::scale (const int s) {
	Gmpq rational (s);
	scale (rational);
	return;
};

void Rational_polynomial_1::scale (const CGAL::Gmpz &s) {
	Gmpq rational (s.mpz());
	scale (rational);
	return;
};

void Rational_polynomial_1::scale (const CGAL::Gmpq &s) {
	mpq_t temp;
	mpq_init (temp);
	for (int i=0; i<=degree; ++i) {
		int ind = calc_index (i);
		mpq_mul (temp, coef[ind], s.mpq());
		mpq_clear (coef[ind]);
		mpq_init (coef[ind]);
		mpq_set (coef[ind], temp);
	}
	mpq_clear(temp);
	return;
};

int Rational_polynomial_1::get_degree () const {
	return (int)degree;
};

int Rational_polynomial_1::get_number_of_monomials () const {
	return nm;
};

mpq_t* Rational_polynomial_1::get_coefs () const {
	return coef;
};

void Rational_polynomial_1::get_coef (int pow_x, mpq_t *c) const {
	mpq_set (*c, coef[calc_index (pow_x)]);
	return;
};

CGAL::MpfiInterval Rational_polynomial_1::eval (const CGAL::MpfiInterval &x) const {
	MpfiInterval result(0);
	MpfiInterval x_pow(1);
	for (int i=0; i<=degree; ++i) {
		// invariant at this point: x_pow = x^i
		result += x_pow * coef[calc_index (i)];
		x_pow *= x;
	}
	return result;
};

CGAL::Gmpq Rational_polynomial_1::eval (const CGAL::Gmpq &x) const {
	mpq_t result, x_pow, temp1, temp2;
	mpq_init (result);	// it's 0 now
	mpq_init (x_pow);
	mpq_set_si (x_pow, 1, 1);	// x^0 = 1
	mpq_init (temp1);
	mpq_init (temp2);
	for (int i=0; i<=degree; ++i) {	// invariant: x_pow = x^i
		mpq_mul (temp1, coef[calc_index (i)], x_pow);
		mpq_set (temp2, result);
		mpq_add (result, temp1, temp2);
		mpq_mul (temp1, x_pow, x.mpq());
		mpq_set (x_pow, temp1);	// to use it in the next iteration
	}
	mpq_clear (x_pow);
	mpq_clear (temp1);
	mpq_clear (temp2);
	CGAL::Gmpq ret(result);
	mpq_clear (result);
	return ret;
};

// derive polynomial
Rational_polynomial_1 Rational_polynomial_1::derive () const {
	mpq_t coef_old, coef_new, temp_x;
	mpq_init (coef_old);
	mpq_init (coef_new);
	mpq_init (temp_x);
	Rational_polynomial_1 derivative (degree-1);
	for (int x=1; x<=degree; ++x) {
		mpq_set_si (temp_x, (long int)x, (unsigned long int)1);
		get_coef (x, &coef_old);
		mpq_mul (coef_new, coef_old, temp_x);
		derivative.set_coef (x-1, coef_new);
	}
	mpq_clear (coef_old);
	mpq_clear (coef_new);
	mpq_clear (temp_x);
	return derivative;
};

Rational_polynomial_1 Rational_polynomial_1::operator= (const Rational_polynomial_1 &p) {
	// destroy the current data
	for (int i=0; i<nm; ++i)
		mpq_clear (coef[i]);
	free (coef);
	// copy data from p
	degree = p.degree;
	nm = p.nm;
	coef = (mpq_t*)malloc (sizeof(mpq_t)*nm);
	for (int i=0; i<nm; ++i) {
		mpq_init (coef[i]);
		mpq_set (coef[i], p.coef[i]);
		}
	return *this;
};

std::ostream& Rational_polynomial_1::show (std::ostream &s) const {
	bool printed = false;
	if (degree == 0) {
		s << coef[0];
		return s;
	}
	for (int i=0; i<=degree; ++i) {
		if (mpq_sgn (coef[i]) != 0) {
			if (printed && (mpq_sgn (coef[i]) == 1))
				s << "+";
			printed = true;
			bool flag = false;
			if ((mpq_cmp_si (coef[i], -1, 1) == 0) && (i != degree))
				s << "-";
			else
				if ((mpq_cmp_ui (coef[i], 1, 1) != 0) || (i == degree)) {
					flag = true;
					s << coef[i];
				}
			if (degree != i) {
				if (flag)
					s << "*";
				s << "x";
				if (degree != i+1)
					s << "^" << degree-i;
			}
		}
	}
	if (!printed)
		s << "0";
	return s;
};

// overcharging
Rational_polynomial_1 Rational_polynomial_1::operator- () const {
	Rational_polynomial_1 opposite (degree);	// the opposite has the same degree
	mpq_t temp;
	mpq_init (temp);
	for (int i=0; i<nm; ++i) {
		mpq_neg (temp, coef[calc_index (i)]);
		opposite.set_coef (i, temp);
	}
	mpq_clear (temp);
	return opposite;
};

Rational_polynomial_1 Rational_polynomial_1::operator+ (const Rational_polynomial_1 &s) const {
	int sd = s.get_degree ();
	int minord, majord;	// the minor and major of both degrees
	bool am_i_bigger;	// is *this' degree bigger than s' degree?
	mpq_t temp1, temp2;

	if (sd < degree) {
		am_i_bigger = true;
		minord = sd;
		majord = degree;
	} else {
		am_i_bigger = false;
		minord = degree;
		majord = sd;
	}
	
	Rational_polynomial_1 sum (degree>sd?degree:sd);

	mpq_init (temp1);
	mpq_init (temp2);

	for (int i=0; i<=minord; ++i) {
		s.get_coef (i, &temp1);
		mpq_add (temp2, temp1, coef[calc_index (i)]);
		sum.set_coef (i, temp2);
	}

	mpq_clear (temp1);

	if (am_i_bigger)
		for (int i=minord+1; i<=majord; ++i)
			sum.set_coef (i, coef[calc_index (i)]);
	else
		for (int i=minord+1; i<=majord; ++i) {
			s.get_coef (i, &temp2);
			sum.set_coef (i, temp2);
		}

	mpq_clear (temp2);
	return sum;
};

void Rational_polynomial_1::operator+= (const Rational_polynomial_1 &s) {
	Rational_polynomial_1 aux (*this);
	*this = aux + s;
	return;
};

Rational_polynomial_1 Rational_polynomial_1::operator- (const Rational_polynomial_1 &s) const {
	return *this + (-s);
};

void Rational_polynomial_1::operator-= (const Rational_polynomial_1 &s) {
	Rational_polynomial_1 aux (*this);
	*this = aux - s;
	return;
};

// multiplies the polynomial by x^shiftn (precondition: shiftn>=0)
void Rational_polynomial_1::shift (int shiftn) {
	CGAL_assertion (shiftn>=0);
	if (shiftn == 0)
		return;
	mpq_t *new_coef = (mpq_t*)malloc (sizeof(mpq_t)*(nm+shiftn));
	for (int i=0; i<nm; ++i) {
		mpq_init (new_coef[i]);
		mpq_set (new_coef[i], coef[i]);
		mpq_clear (coef[i]);
	}
	free (coef);
	for (int i=nm; i<nm+shiftn; ++i)
		mpq_init (new_coef[i]);
	degree += shiftn;
	nm = degree+1;
	coef = new_coef;
	return;
};

Rational_polynomial_1 Rational_polynomial_1::operator* (const Rational_polynomial_1 &f) const {
	Rational_polynomial_1 product;
	mpq_t *f_coefs;
	f_coefs = f.get_coefs ();
	int df = f.get_degree ();
	for (int i=0; i<=df; ++i) {
		Rational_polynomial_1 partial (*this);
		partial.scale (f_coefs[i]);
		partial.shift (df-i);
		product += partial;
	}
	return product;
};

void Rational_polynomial_1::operator*= (const Rational_polynomial_1 &f) {
	Rational_polynomial_1 aux (*this);
	*this = aux * f;
	return;
};

bool Rational_polynomial_1::operator== (const Rational_polynomial_1 &p) const {
	mpq_t *p_coef;
	if (degree != p.get_degree ())
		return false;
	p_coef = p.get_coefs ();
	for (int i=0; i<=degree; ++i)
		if (mpq_cmp (coef[i], p_coef[i]) != 0)
			return false;
	return true;
};

std::ostream& operator<< (std::ostream &o, const Rational_polynomial_1 &p) {
	return p.show (o);
}

CGAL_END_NAMESPACE
