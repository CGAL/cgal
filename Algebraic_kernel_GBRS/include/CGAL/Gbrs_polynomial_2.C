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
#include <CGAL/MpfiInterval.h>
#include <iostream>
#include <gmp.h>

CGAL_BEGIN_NAMESPACE

// constructors
Rational_polynomial_2::Rational_polynomial_2 () {
	degree_x = 0;
	degree_y = 0;
	coef = (mpq_t**)malloc (sizeof(mpq_t*));
	coef[0] = (mpq_t*)malloc (sizeof(mpq_t));
	mpq_init (coef[0][0]);
};

Rational_polynomial_2::Rational_polynomial_2 (unsigned int dx,
		unsigned int dy) {
	degree_x = (int)dx;
	degree_y = (int)dy;
	coef = (mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	for (int i=0; i<=degree_x; ++i) {
		coef[i] = (mpq_t*)malloc ((degree_y+1) * sizeof(mpq_t));
		for (int j=0; j<=degree_y; ++j)
			mpq_init (coef[i][j]);
	}
};

Rational_polynomial_2::Rational_polynomial_2 (int dx, int dy) {
	degree_x = dx<0?0:dx;
	degree_y = dy<0?0:dy;
	coef = (mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	for (int i=0; i<=degree_x; ++i) {
		coef[i] = (mpq_t*)malloc ((degree_y+1) * sizeof(mpq_t));
		for (int j=0; j<=degree_y; ++j)
			mpq_init (coef[i][j]);
	}
};

Rational_polynomial_2::Rational_polynomial_2 (const Rational_polynomial_2 &p) {
	degree_x = p.degree_x;
	degree_y = p.degree_y;
	mpq_t **p_coef = p.get_coefs ();
	coef = (mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	// we have to copy the contents, not just the pointer
	for (int i=0; i<=degree_x; ++i) {
		coef[i] = (mpq_t*)malloc ((degree_y+1) * sizeof(mpq_t));
		for (int j=0; j<=degree_y; ++j) {
			mpq_init (coef[i][j]);
			mpq_set (coef[i][j], p_coef[i][j]);
		}
	}
};

// destructor
Rational_polynomial_2::~Rational_polynomial_2 () {
	for (int i=0; i<=degree_x; ++i) {
		for (int j=0; j<=degree_y; ++j)
			mpq_clear (coef[i][j]);
		free (coef[i]);
	}
	free (coef);
};

// other functions
inline int Rational_polynomial_2::get_degree_x () const {
	return degree_x;
};

inline int Rational_polynomial_2::get_degree_y () const {
	return degree_y;
};

inline mpq_t** Rational_polynomial_2::get_coefs () const {
	return coef;
};

void Rational_polynomial_2::get_coef (int pow_x, int pow_y, mpq_t *c) const {
	mpq_set (*c, coef[pow_x][pow_y]);
	return;
};

// functions for setting polynomial coefficients
void Rational_polynomial_2::set_coef (int pow_x, int pow_y, const mpq_t &q) {
	mpq_set (coef[pow_x][pow_y], q);
};

void Rational_polynomial_2::set_coef (int pow_x, int pow_y, const mpz_t &z) {
	mpq_set_z (coef[pow_x][pow_y], z);
};

void Rational_polynomial_2::set_coef (int pow_x, int pow_y, const CGAL::Gmpq &q) {
	mpq_set (coef[pow_x][pow_y], q.mpq());
};

void Rational_polynomial_2::set_coef (int pow_x, int pow_y, const CGAL::Gmpz &z) {
	mpq_set_z (coef[pow_x][pow_y], z.mpz());
};

void Rational_polynomial_2::set_coef (int pow_x, int pow_y, int c) {
	mpq_set_si (coef[pow_x][pow_y], (long int)c, (unsigned long int)1);
};

void Rational_polynomial_2::set_coef (int pow_x, int pow_y, unsigned int c) {
	mpq_set_ui (coef[pow_x][pow_y], (unsigned long int)c, (unsigned long int)1);
};

void Rational_polynomial_2::set_coef (int pow_x, int pow_y, long int c) {
	mpq_set_si (coef[pow_x][pow_y], c, (unsigned long int)1);
};

void Rational_polynomial_2::set_coef (int pow_x, int pow_y, unsigned long int c) {
	mpq_set_ui (coef[pow_x][pow_y], c, (unsigned long int)1);
};

Rational_polynomial_2 Rational_polynomial_2::operator= (const Rational_polynomial_2 &p) {
	// destroy the current data
	for (int i=0; i<=degree_x; ++i) {
		for (int j=0; j<=degree_y; ++j)
			mpq_clear (coef[i][j]);
		free (coef[i]);
	}
	free (coef);
	// copy data from p
	degree_x = p.degree_x;
	degree_y = p.degree_y;
	mpq_t **p_coef = p.get_coefs ();
	coef = (mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	// we have to copy the contents, not just the pointer
	for (int i=0; i<=degree_x; ++i) {
		coef[i] = (mpq_t*)malloc ((degree_y+1) * sizeof(mpq_t));
		for (int j=0; j<=degree_y; ++j) {
			mpq_init (coef[i][j]);
			mpq_set (coef[i][j], p_coef[i][j]);
		}
	}
	return *this;
};

Rational_polynomial_2 Rational_polynomial_2::operator- () const {
	Rational_polynomial_2 opposite (degree_x, degree_y);
	mpq_t temp;
	mpq_init (temp);
	for (int i=0; i<=degree_x; ++i)
		for (int j=0; j<=degree_y; ++j) {
			mpq_neg (temp, coef[i][j]);
			opposite.set_coef (i, j, temp);
		}
	mpq_clear (temp);
	return opposite;
};

Rational_polynomial_2 Rational_polynomial_2::operator+ (const Rational_polynomial_2 &s) const {
	int sx = s.get_degree_x ();
	int sy = s.get_degree_y ();
	int minorx, majorx, minory, majory;	// the minor and major of both degrees
	bool am_i_bigger_x, am_i_bigger_y;	// is *this' degree bigger?
	mpq_t temp1, temp2, temp3;

	if (sx < degree_x) {
		am_i_bigger_x = true;
		minorx = sx;
		majorx = degree_x;
	} else {
		am_i_bigger_x = false;
		minorx = degree_x;
		majorx = sx;
	}

	if (sy < degree_y) {
		am_i_bigger_y = true;
		minory = sy;
		majory = degree_y;
	} else {
		am_i_bigger_y = false;
		minory = degree_y;
		majory = sy;
	}

	Rational_polynomial_2 sum (am_i_bigger_x?degree_x:sx,
			am_i_bigger_y>sy?degree_y:sy);

	mpq_init (temp1);
	mpq_init (temp2);
	mpq_init (temp3);
	int i, j;

	if (am_i_bigger_x)
		// copiar a "sum" los X que yo tengo y él no
		for (i=minorx+1; i<=majorx; ++i)
			for (j=0; j<=minory; ++j) {
				get_coef (i, j, &temp1);
				sum.set_coef (i, j, temp1);
			}
	else
		// copiar a "sum" los X que él tiene y yo no
		for (i=minorx+1; i<=majorx; ++i)
			for (j=0; j<=minory; ++j) {
				s.get_coef (i, j, &temp1);
				sum.set_coef (i, j, temp1);
			}

	if (am_i_bigger_y)
		// copiar a "sum" los Y que yo tengo y él no, pero no copiar de nuevo los X que ya copiamos
		for (i=0; i<=minorx; ++i)
			for (j=minory+1; j<=majory; ++j) {
				get_coef (i, j, &temp1);
				sum.set_coef (i, j, temp1);
			}
	else
		// copiar a "sum" los Y que él tiene y yo no, pero no copiar de nuevo los X que ya copiamos
		for (i=0; i<=minorx; ++i)
			for (j=minory+1; j<=majory; ++j) {
				s.get_coef (i, j, &temp1);
				sum.set_coef (i, j, temp1);
			}

	// ahora copiar las sumas de los que ambos tenemos en común
	for (i=0; i<=minorx; ++i)
		for (j=0; j<=minory; ++j) {	// XXX mmm... veremos
		//for (j=minory+1; j<=majory; ++j) {
				get_coef (i, j, &temp1);
				s.get_coef (i, j, &temp2);	// XXX acá es el segfault... existe s[i][j] ?
				mpq_add (temp3, temp2, temp1);
				sum.set_coef (i, j, temp3);
		}

	mpq_clear (temp1);
	mpq_clear (temp2);
	mpq_clear (temp3);

	return sum;
};

Rational_polynomial_2& Rational_polynomial_2::operator+= (const Rational_polynomial_2 &s) {
	Rational_polynomial_2 aux (*this);
	*this = aux + s;
	return *this;
};

Rational_polynomial_2 Rational_polynomial_2::operator- (const Rational_polynomial_2 &s) const {
	return *this + (-s);
};

Rational_polynomial_2& Rational_polynomial_2::operator-= (const Rational_polynomial_2 &s) {
	Rational_polynomial_2 aux (*this);
	*this = aux - s;
	return *this;
};

// scaling: this is the same as operator* (int)
void Rational_polynomial_2::scale (const int s) {
	Gmpq rational (s);
	scale (rational);
	return;
};

void Rational_polynomial_2::scale (const mpz_t &s) {
	Gmpq rational (s);
	scale (rational);
	return;
};

void Rational_polynomial_2::scale (const mpq_t &s) {
	Gmpq rational (s);
	scale (rational);
	return;
};

void Rational_polynomial_2::scale (const CGAL::Gmpz &s) {
	Gmpq rational (s.mpz());
	scale (rational);
	return;
};

void Rational_polynomial_2::scale (const CGAL::Gmpq &s) {
	mpq_t temp;
	mpq_init (temp);
	for (int i=0; i<=degree_x; ++i)
		for (int j=0; j<=degree_y; ++j) {
			mpq_mul (temp, coef[i][j], s.mpq());
			mpq_set (coef[i][j], temp);
		}
	mpq_clear(temp);
	return;
};

// multiplies the polynomial by x^shift_x * y^shift_y
// (preconditions: shift_[xy] >= 0)
void Rational_polynomial_2::shift (int shift_x, int shift_y) {
	CGAL_assertion ((shift_x >= 0) && (shift_y >= 0));
	if ((shift_x == 0) && (shift_y == 0))
		return;
	int i, j;
	degree_x += shift_x;
	degree_y += shift_y;
	mpq_t **new_coef =
		(mpq_t**)malloc (sizeof(mpq_t*)*(degree_x+1));
	// init the coefficients that will be zero
	for (i=0; i<shift_x; ++i) {
		new_coef[i] = (mpq_t*)malloc ((degree_y+1)*sizeof(mpq_t));
		for (j=0; j<=degree_y; ++j)
			mpq_init (new_coef[i][j]);
	}
	// copy the rest
	for (i=shift_x; i<=degree_x; ++i) {
		new_coef[i] =
			(mpq_t*)malloc ((degree_y+1)*sizeof(mpq_t));
		for (j=0; j<shift_y; ++j)
			mpq_init (new_coef[i][j]);
		for (j=shift_y; j<=degree_y; ++j) {
			mpq_init (new_coef[i][j]);
			mpq_set (new_coef[i][j], coef[i-shift_x][j-shift_y]);
			mpq_clear (coef[i-shift_x][j-shift_y]);
		}
		free (coef [i-shift_x]);	// free the column
	}
	free (coef);
	coef = new_coef;
	return;
};

// how to multiply:
// 1. create a polynomial to store the result (with zeros)
// 2. multiply *this by every monomial of &f (with shift)
// 3. sum all
Rational_polynomial_2 Rational_polynomial_2::operator* (const Rational_polynomial_2 &f) const {
	Rational_polynomial_2 product;
	mpq_t **f_coefs;
	f_coefs = f.get_coefs ();
	int xf = f.get_degree_x ();
	int yf = f.get_degree_y ();
	for (int i=0; i<=xf; ++i)
		for (int j=0; j<=yf; ++j) {
			Rational_polynomial_2 partial (*this);
			partial.scale (f_coefs[i][j]);
			partial.shift (i, j);
			product += partial;
		}
	return product;
};

// TODO: consider the case where the two polynomials have different degrees but
// are still equal
bool Rational_polynomial_2::operator== (const Rational_polynomial_2 &p) const {
	if ((degree_x != p.get_degree_x ()) || (degree_y != p.get_degree_y ()))
		return false;
	mpq_t **p_coef;
	p_coef = p.get_coefs ();
	for (int i=0; i<=degree_x; ++i)
		for (int j=0; j<=degree_y; ++j)
			if (mpq_cmp (coef[i][j], p_coef[i][j]) != 0)
				return false;
	return true;
};

inline bool Rational_polynomial_2::operator!= (const Rational_polynomial_2 &p) const {
	return !(*this == p);
};

std::ostream& Rational_polynomial_2::show (std::ostream &s) const {
	bool printed = false;
	bool zero;
	for (int i=degree_x; i>=0; --i)
		for (int j=degree_y; j>=0; --j) {
			zero = false;
			switch (mpq_sgn (coef[i][j])) {
				case 0:	// the coefficient is 0
					zero = true;	// we don't print it
					break;
				case 1:	// the coefficient is greater than 0
					if (printed)
						s << "+";
					break;
				default:	// the coefficient is lesser than 0
					break;
			}
			if (!zero) {
				if ((mpq_cmp_ui (coef[i][j], 1, 1) != 0) &&
						!((i==0)&&(j==0))) {
					s << coef[i][j];
					if (i > 0)
						s << "*";
				}
				printed = true;
				if (i > 0) {
					s << "x";
					if (i > 1)
						s << "^" << i;
				}
				if (j > 0) {
					s << "y";
					if (j > 1)
						s << "^" << j;
				}
			}
		}
	if (!printed)
		s << "0";
	return s;
};

Rational_polynomial_2& Rational_polynomial_2::operator*= (const Rational_polynomial_2 &f) {
	Rational_polynomial_2 aux (*this);
	*this = aux * f;
	return *this;
};

template <class T>
Rational_polynomial_2 Rational_polynomial_2::operator* (const T &n) const {
	Rational_polynomial_2 r (*this);
	r.scale (n);
	return r;
};

std::ostream& operator<< (std::ostream &o, const Rational_polynomial_2 &p) {
	return p.show (o);
};

template <class T>
inline Rational_polynomial_2 operator* (const T &n, Rational_polynomial_2 &p) {
		return (p*n);
};

CGAL_END_NAMESPACE
