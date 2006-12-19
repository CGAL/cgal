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

#ifndef CGAL_GBRS_POLYNOMIAL_1_H
#define CGAL_GBRS_POLYNOMIAL_1_H

#include <iostream>
#include <vector>
#include <CGAL/assertions.h>
#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>

CGAL_BEGIN_NAMESPACE

class Algebraic_1;

/*typedef std::vector<Algebraic_1> rootvector;*/

class Rational_polynomial_1 {
	private:
		int degree;
		mpz_t* coef;
		bool solved;
		/*rootvector roots;*/
	public:
		// copy constructor and copy assignement operator
		Rational_polynomial_1(const Rational_polynomial_1&);
		Rational_polynomial_1& operator=(const Rational_polynomial_1&);
		// other constructors and destructor
		Rational_polynomial_1 ();
		Rational_polynomial_1 (unsigned int);
		Rational_polynomial_1 (int);
		Rational_polynomial_1(const mpq_t&);
		~Rational_polynomial_1 ();
		// member functions
		void set_degree (int);
		void set_coef (int, const mpz_t &);
		void set_coef (int, const CGAL::Gmpz &);
		void set_coef (int, int);
		void set_coef (int, unsigned int);
		void set_coef (int, long int);
		void set_coef (int, unsigned long int);
		int get_degree () const;
		int get_number_of_monomials () const;
		mpz_t* get_coefs () const;
		void get_coef (int, mpz_t *) const;
		void set_solved();
		void clear_solved();
		bool get_solved()const;
		/*rootvector get_roots()const;
		void set_root(const Algebraic_1&);*/
		//CGAL::Algebraic_1 eval_alg(const CGAL::Algebraic_1&)const;
		CGAL::Gmpz eval (const CGAL::Gmpz &) const;
		void eval_mpfr(mpfr_t&,const mpfr_t&,mp_prec_t)const;
		void eval_mpfi(mpfi_ptr,mpfi_srcptr)const;
		Rational_polynomial_1 derive () const;
		std::ostream& show (std::ostream &) const;
		Rational_polynomial_1 operator- () const;
		Rational_polynomial_1 operator+ (const Rational_polynomial_1 &) const;
		Rational_polynomial_1& operator+= (const Rational_polynomial_1 &);
		Rational_polynomial_1 operator- (const Rational_polynomial_1 &) const;
		Rational_polynomial_1& operator-= (const Rational_polynomial_1 &);
		/*Rational_polynomial_1& scale_and_shift(const mpz_t&,int);*/
		Rational_polynomial_1 operator* (const Rational_polynomial_1 &) const;
		template <class T> Rational_polynomial_1 operator* (const T &) const;
		Rational_polynomial_1& operator*=(const Rational_polynomial_1 &);
		Rational_polynomial_1& operator*=(const mpz_t &);
		Rational_polynomial_1& operator*=(const CGAL::Gmpz &);
		template<class T>Rational_polynomial_1& operator*=(const T&);
		bool operator== (const Rational_polynomial_1 &) const;
};

std::ostream& operator<< (std::ostream &, const Rational_polynomial_1 &);

// /////////////////
// inline functions
inline void Rational_polynomial_1::set_coef(int pow_x,const mpz_t &z){
	mpz_set(coef[pow_x],z);};
inline void Rational_polynomial_1::set_coef(int pow_x,const CGAL::Gmpz &z){
	mpz_set(coef[pow_x],z.mpz());};
inline void Rational_polynomial_1::set_coef(int pow_x,int c){
	mpz_set_si(coef[pow_x],(long int)c);};
inline void Rational_polynomial_1::set_coef(int pow_x,unsigned int c){
	mpz_set_ui(coef[pow_x],(unsigned long int)c);};
inline void Rational_polynomial_1::set_coef(int pow_x,long int c){
	mpz_set_si(coef[pow_x],c);};
inline void Rational_polynomial_1::set_coef(int pow_x,unsigned long int c){
	mpz_set_ui(coef[pow_x],c);};
inline int Rational_polynomial_1::get_degree()const{return(int)degree;};
inline int Rational_polynomial_1::get_number_of_monomials()const{
	return(degree+1);};
inline mpz_t* Rational_polynomial_1::get_coefs()const{return coef;};
inline void Rational_polynomial_1::get_coef(int pow_x,mpz_t *c)const{
	mpz_set(*c,coef[pow_x]);};
inline bool Rational_polynomial_1::get_solved()const{return solved;};
/*inline rootvector Rational_polynomial_1::get_roots()const{return roots;};
inline void Rational_polynomial_1::set_root(const Algebraic_1 &r){
	roots.push_back(r);};*/

inline std::ostream& operator<<(std::ostream &o,const Rational_polynomial_1 &p)
{return p.show (o);}

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Gbrs_polynomial_1_impl.h>
#endif	// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif	// CGAL_GRBS_POLYNOMIAL_1_H
