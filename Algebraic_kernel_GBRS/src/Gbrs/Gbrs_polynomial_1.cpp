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
#include <CGAL/Gbrs_polynomial_1.h>
#include <CGAL/Gbrs_algebraic_1.h>
#include <iostream>
#include <gmp.h>

CGAL_BEGIN_NAMESPACE

// constructors
Rational_polynomial_1::Rational_polynomial_1():degree(0),solved(false){
	coef=(mpz_t*)malloc(sizeof(mpz_t));
	mpz_init(coef[0]);
	//roots.clear();
};

Rational_polynomial_1::Rational_polynomial_1(unsigned int d):solved(false){
	degree=(int)d;
	coef=(mpz_t*)malloc(sizeof(mpz_t)*(degree+1));
	for(int i=0;i<degree+1;++i)
		mpz_init(coef[i]);
	//roots.clear();
};

Rational_polynomial_1::Rational_polynomial_1(int d):solved(false){
	degree=d<0?0:d;
	coef=(mpz_t*)malloc(sizeof(mpz_t)*(degree+1));
	for(int i=0;i<degree+1;++i)
		mpz_init(coef[i]);
	//roots.clear();
};

Rational_polynomial_1::Rational_polynomial_1(const Rational_polynomial_1 &p){
	degree=p.degree;
	mpz_t *p_coef=p.get_coefs();
	coef=(mpz_t*)malloc(sizeof(mpz_t)*(degree+1));
	// we have to copy the contents, not just the pointer
	for(int i=0;i<degree+1;++i){
		mpz_init(coef[i]);
		mpz_set(coef[i],p_coef[i]);
		}
	solved=p.get_solved();
	//roots=p.get_roots();
};

// destructor
Rational_polynomial_1::~Rational_polynomial_1 () {
	for(int i=0;i<degree+1;++i)
		mpz_clear(coef[i]);
	free(coef);
	/*if(solved)
		roots.clear();*/
};

void Rational_polynomial_1::set_degree(int d){	// dangerous function!
	for(int i=0;i<degree+1;++i)	// free the old coefficients
		mpz_clear(coef[i]);
	free(coef);
	degree=d;
	coef=(mpz_t*)malloc(sizeof(mpz_t)*(degree+1));
	for(int i=0;i<degree+1;++i)
		mpz_init(coef[i]);
	solved=false;
	return;
};

// functions for setting polynomial coefficients
void Rational_polynomial_1::set_coef (int pow_x, const mpz_t &z) {
	mpz_set (coef[calc_index (pow_x)], z);
};

void Rational_polynomial_1::set_coef (int pow_x, const CGAL::Gmpz &z) {
	mpz_set (coef[calc_index (pow_x)], z.mpz());
};

void Rational_polynomial_1::set_coef (int pow_x, int c) {
	mpz_set_si (coef[calc_index (pow_x)], (long int)c);
};

void Rational_polynomial_1::set_coef (int pow_x, unsigned int c) {
	mpz_set_ui (coef[calc_index (pow_x)], (unsigned long int)c);
};

void Rational_polynomial_1::set_coef (int pow_x, long int c) {
	mpz_set_si (coef[calc_index (pow_x)], c);
};

void Rational_polynomial_1::set_coef (int pow_x, unsigned long int c) {
	mpz_set_ui (coef[calc_index (pow_x)], c);
};

// scaling
void Rational_polynomial_1::scale (const int s) {
	Gmpz rational (s);
	scale (rational);
	return;
};

void Rational_polynomial_1::scale (const mpz_t &s) {
	Gmpz rational (s);
	scale (rational);
	return;
};

void Rational_polynomial_1::scale (const CGAL::Gmpz &s) {
	for(int i=0;i<=degree;++i)
		mpz_mul(coef[i],coef[i],s.mpz());
	/*mpz_t temp;
	mpz_init (temp);
	for (int i=0; i<=degree; ++i) {
		int ind = calc_index (i);
		mpz_mul (temp, coef[ind], s.mpz());
		mpz_clear (coef[ind]);
		mpz_init (coef[ind]);
		mpz_set (coef[ind], temp);
	}
	mpz_clear(temp);
	return;*/
};

void Rational_polynomial_1::get_coef (int pow_x, mpz_t *c) const {
	mpz_set (*c, coef[calc_index (pow_x)]);
	return;
};

void Rational_polynomial_1::set_solved(){
	solved=true;
	//roots.clear();
};

void Rational_polynomial_1::clear_solved(){
	/*if(solved)
		roots.clear();*/
	solved=false;
	return;
};

/*CGAL::Algebraic_1 Rational_polynomial_1::eval (const CGAL::Algebraic_1 &x) const {
	Algebraic_1 result(0);
	Algebraic_1 x_pow(1);
	for (int i=0; i<=degree; ++i) {
		// invariant at this point: x_pow = x^i
		result += x_pow * coef[calc_index (i)];
		x_pow *= x;
	}
	return result;
};*/

CGAL::Gmpz Rational_polynomial_1::eval(const CGAL::Gmpz &x)const{
	mpz_t result,x_pow,temp;
	mpz_init(result);	// it's 0 now
	mpz_init_set_si(x_pow,1);	// x^0 = 1
	mpz_init(temp);
	for(int i=0;i<=degree;++i){	// invariant: x_pow = x^i
		mpz_mul(temp,coef[calc_index(i)],x_pow);
		mpz_add(result,temp,result);
		mpz_mul(x_pow,x_pow,x.mpz());	// for the next iteration
	}
	mpz_clear(x_pow);
	mpz_clear(temp);
	CGAL::Gmpz ret(result);
	mpz_clear(result);
	return ret;
};

void Rational_polynomial_1::eval_mpfr
(mpfr_t &result,const mpfr_t &x,mp_prec_t prec)const{
	mpfr_t x_pow,temp;
	mp_prec_t prec_r=mpfr_get_prec(result);
	mpfr_inits2(prec<prec_r?prec_r:prec,x_pow,temp,NULL);
	mpfr_set_ui(x_pow,1,GMP_RNDN);
	mpfr_set_ui(result,0,GMP_RNDN);
	for(int i=0;i<=degree;++i){
		mpfr_mul_z(temp,x_pow,coef[calc_index(i)],GMP_RNDN);
		mpfr_add(result,temp,result,GMP_RNDN);
		mpfr_mul(x_pow,x_pow,x,GMP_RNDN);
	}
	mpfr_clears(x_pow,temp,NULL);
	return;
};

// I think RS should do this
void Rational_polynomial_1::eval_mpfi(mpfi_t &result,const mpfi_t &x)const{
	mpfi_t x_pow,temp;
	mpfi_init_set_ui(x_pow,1);
	mpfi_init(temp);
	mpfi_set_ui(result,0);
	for(int i=0;i<=degree;++i){
		mpfi_mul_z(temp,x_pow,coef[calc_index(i)]);
		mpfi_add(result,temp,result);
		mpfi_mul(x_pow,x_pow,x);
	}
	return;
};

Rational_polynomial_1 Rational_polynomial_1::derive () const {
	mpz_t coef_old, coef_new, temp_x;
	mpz_init (coef_old);
	mpz_init (coef_new);
	mpz_init (temp_x);
	Rational_polynomial_1 derivative (degree-1);
	for (int x=1; x<=degree; ++x) {
		mpz_set_si (temp_x, (long int)x);
		get_coef (x, &coef_old);
		mpz_mul (coef_new, coef_old, temp_x);
		derivative.set_coef (x-1, coef_new);
	}
	mpz_clear (coef_old);
	mpz_clear (coef_new);
	mpz_clear (temp_x);
	return derivative;
};

Rational_polynomial_1& Rational_polynomial_1::operator= (const Rational_polynomial_1 &p) {
	// destroy the current data
	for (int i=0; i<degree+1; ++i)
		mpz_clear (coef[i]);
	free (coef);
	// copy data from p
	degree = p.degree;
	coef = (mpz_t*)malloc (sizeof(mpz_t)*(degree+1));
	for (int i=0; i<degree+1; ++i) {
		mpz_init (coef[i]);
		mpz_set (coef[i], p.coef[i]);
		}
	solved=p.get_solved();
	return *this;
};

std::ostream& Rational_polynomial_1::show (std::ostream &s) const {
	bool printed = false;
	if (degree == 0) {
		s << coef[0];
		return s;
	}
	for (int i=degree; i>=0; --i) {
		if (mpz_sgn (coef[i]) != 0) {
			if (printed && (mpz_sgn (coef[i]) == 1))
				s << "+";
			printed = true;
			bool flag = false;
			if ((mpz_cmp_si (coef[i], -1) == 0) && (i != 0))
				s << "-";
			else
				if ((mpz_cmp_ui (coef[i], 1) != 0) || (i == 0)) {
					flag = true;
					s << coef[i];
				}
			if (0 != i) {
				if (flag)
					s << "*";
				s << "x";
				if (0 != i-1)
					s << "^" << i;
			}
		}
	}
	if (!printed)
		s << "0";
#ifdef CGAL_RS_DEBUG
	s<<" [ d="<<degree<<" ";
	s<<"[ ";
	for (int i=0; i<degree+1; ++i)
		s<<coef[i]<<" ";
	s<<"] ]";
#endif
	return s;
};

Rational_polynomial_1 Rational_polynomial_1::operator- () const {
	Rational_polynomial_1 opposite (degree);	// the opposite has the same degree
	mpz_t temp;
	mpz_init (temp);
	for (int i=0; i<degree+1; ++i) {
		mpz_neg (temp, coef[calc_index (i)]);
		opposite.set_coef (i, temp);
	}
	mpz_clear (temp);
	return opposite;
};

Rational_polynomial_1 Rational_polynomial_1::operator+ (const Rational_polynomial_1 &s) const {
	int sd = s.get_degree ();
	int minord, majord;	// the minor and major of both degrees
	bool am_i_bigger;	// is *this' degree bigger than s' degree?
	mpz_t temp1, temp2;

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

	mpz_init (temp1);
	mpz_init (temp2);

	for (int i=0; i<=minord; ++i) {
		s.get_coef (i, &temp1);
		mpz_add (temp2, temp1, coef[calc_index (i)]);
		sum.set_coef (i, temp2);
	}

	mpz_clear (temp1);

	if (am_i_bigger)
		for (int i=minord+1; i<=majord; ++i)
			sum.set_coef (i, coef[calc_index (i)]);
	else
		for (int i=minord+1; i<=majord; ++i) {
			s.get_coef (i, &temp2);
			sum.set_coef (i, temp2);
		}

	mpz_clear (temp2);
	return sum;
};

Rational_polynomial_1& Rational_polynomial_1::operator+= (const Rational_polynomial_1 &s) {
	Rational_polynomial_1 aux (*this);
	*this = aux + s;
	solved=false;
	return *this;
};

Rational_polynomial_1 Rational_polynomial_1::operator- (const Rational_polynomial_1 &s) const {
	return (*this+(-s));
};

Rational_polynomial_1& Rational_polynomial_1::operator-= (const Rational_polynomial_1 &s) {
	Rational_polynomial_1 aux (*this);
	*this = aux - s;
	solved=false;
	return *this;
};

// multiplies the polynomial by x^shiftn (precondition: shiftn>=0)
void Rational_polynomial_1::shift (int shiftn) {
	CGAL_assertion (shiftn>=0);
	if (shiftn == 0)
		return;
	mpz_t *new_coef = (mpz_t*)malloc (sizeof(mpz_t)*(degree+1+shiftn));
	for (int i=0; i<degree+1; ++i) {
		mpz_init (new_coef[i]);
		mpz_set (new_coef[i], coef[i]);
		mpz_clear (coef[i]);
	}
	free (coef);
	for (int i=degree+1; i<degree+1+shiftn; ++i)
		mpz_init (new_coef[i]);
	degree += shiftn;
	coef = new_coef;
	solved=false;
	return;
};

Rational_polynomial_1 Rational_polynomial_1::operator* (const Rational_polynomial_1 &f) const {
	Rational_polynomial_1 product;
	mpz_t *f_coefs;
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

Rational_polynomial_1& Rational_polynomial_1::operator*= (const Rational_polynomial_1 &f) {
	Rational_polynomial_1 aux (*this);
	*this = aux * f;
	solved=false;
	return *this;
};

// TODO: this function returns false if the two polynomials have different
// degree (it may be the case where the two are equal and one of them has
// the coefficients of higher degree set to zero)
bool Rational_polynomial_1::operator== (const Rational_polynomial_1 &p) const {
	mpz_t *p_coef;
	if (degree != p.get_degree ())
		return false;
	p_coef = p.get_coefs ();
	for (int i=0; i<=degree; ++i)
		if (mpz_cmp (coef[i], p_coef[i]) != 0)
			return false;
	return true;
};

std::ostream& operator<< (std::ostream &o, const Rational_polynomial_1 &p) {
	return p.show (o);
}

CGAL_END_NAMESPACE
