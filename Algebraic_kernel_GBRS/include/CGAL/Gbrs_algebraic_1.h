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

#ifndef CGAL_MPFIINTERVAL_H
#define CGAL_MPFIINTERVAL_H

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <exception>
#include <iostream>
#include <mpfi.h>
#include <mpfi_io.h>
#include <mpfr.h>

CGAL_BEGIN_NAMESPACE

// the exception thrown when it's not clear how to handle inequality
class comparison_overlap_exn : public std::exception {
	public:
	virtual const char* what() const throw () {;
		return "Intervals overlap in comparison";
	};
};

void overlap ();

// The representation of intervals.
class Algebraic_1_rep {
public:
	mpfi_t mpfI;
	Rational_polynomial_1 *poly;
	int nr;
	int mult;
	int rsprec;

	Algebraic_1_rep():poly(NULL),nr(0),mult(0),rsprec(0){mpfi_init(mpfI);}
	~Algebraic_1_rep () { mpfi_clear (mpfI); }

private:
	// Make sure it does not get accidentally copied.
	Algebraic_1_rep (const Algebraic_1_rep &);
	Algebraic_1_rep & operator= (const Algebraic_1_rep &);
};

// The class of the MPFI intervals. It's a model of the RingNumberType concept
class Algebraic_1
: Handle_for<Algebraic_1_rep>,
	/*boost::field_operators1<Algebraic_1,*/
		boost::field_operators2<Algebraic_1, int,
		boost::field_operators2<Algebraic_1, mpz_t,
		boost::field_operators2<Algebraic_1, mpq_t,
		boost::field_operators2<Algebraic_1, Gmpz,
		boost::field_operators2<Algebraic_1, Gmpq,
		boost::field_operators2<Algebraic_1, mpfr_t > > > > > > /*>*/ {
			// XXX: the functions provided by BOOST (or supposed
			// to be) are commented

	typedef Handle_for<Algebraic_1_rep> Base;
public:
	typedef CGAL::Tag_false	Has_gcd;
	typedef CGAL::Tag_true	Has_division;
	typedef CGAL::Tag_true	Has_sqrt;

	typedef CGAL::Tag_true	Has_exact_ring_operations;
	typedef CGAL::Tag_true	Has_exact_division;
	typedef CGAL::Tag_false	Has_exact_sqrt;

	// constructors I
	Algebraic_1 ();
	Algebraic_1 (int);
	Algebraic_1 (unsigned int);
	Algebraic_1 (long int);
	Algebraic_1 (unsigned long int);
	Algebraic_1 (double);
	Algebraic_1 (const mpz_t &);
	Algebraic_1 (const mpq_t &);
	Algebraic_1 (const CGAL::Gmpq &);
	Algebraic_1 (const CGAL::Gmpz &);

	// constructors II
	Algebraic_1 (int, int);
	Algebraic_1 (unsigned int, unsigned int);
	Algebraic_1 (long int, long int);
	Algebraic_1 (unsigned long int, unsigned long int);
	Algebraic_1 (double, double);
	Algebraic_1 (const mpz_t &, const mpz_t &);
	Algebraic_1 (const mpq_t &, const mpq_t &);
	Algebraic_1 (const CGAL::Gmpq &, const CGAL::Gmpq &);
	Algebraic_1 (const CGAL::Gmpz &, const CGAL::Gmpz &);
	Algebraic_1 (const mpfi_t &);
	Algebraic_1 (const Algebraic_1 &);

	// the only interesting constructor
	Algebraic_1 (const mpfi_t &, const Rational_polynomial_1 &,
			const int, const int, const int);

	Algebraic_1& operator= (const long int);
	Algebraic_1& operator= (const mpz_t &);
	Algebraic_1& operator= (const mpq_t &);
	Algebraic_1& operator= (const CGAL::Gmpz &);
	Algebraic_1& operator= (const CGAL::Gmpq &);
	Algebraic_1& operator= (const Algebraic_1 &);

	// destructor
	/* not needed
	~Algebraic_1 ();
	*/

	// functions related to the member data
	inline const mpfi_t & mpfi () const;
	inline mpfi_t & mpfi ();
	inline const Rational_polynomial_1 & pol () const;
	inline Rational_polynomial_1 & pol ();
	inline const int nr () const;
	inline const int mult () const;
	inline const int rsprec () const;
	void clear_pol ();
	void set_pol (const Rational_polynomial_1 &);
	void set_nr (const int);
	void set_mult (const int);
	void set_rsprec (const int);
	inline void set_prec (mp_prec_t);
	inline mp_prec_t get_prec ();
	inline void get_left (mpfr_t &) const;
	inline void get_right (mpfr_t &) const;
	inline void get_endpoints (mpfr_t &, mpfr_t &) const;
	inline bool is_consistent () const;
	inline bool is_point () const;	// are endpoints equal?
	inline bool contains (const int n) const;
	inline bool contains (const mpfr_t &n) const;
	inline bool contains (const mpz_t &n) const;
	inline bool contains (const mpq_t &n) const;
	inline bool contains (const Gmpz &n) const;
	inline bool contains (const Gmpq &n) const;

	// Arithmetic functions required by RingNumberType:
	// 1. comparisons between Algebraic_1's
	// 2. comparisons with int's
	// 3. arithmetic between Algebraic_1's
	// 4. arithmetic with int's
	// 5. accesory functions
	// 6. extra arithmetic functions (not required)
	//
	// Arithmetic functions required by FieldNumberType:
	// 7. division functions
	//
	// Arithmetic functions required by SqrtFieldNumberType:
	// 8. square root
	//
	// I/O functions:
	// 9. <<
	//
	// 10. misc functions
	// 11. all functions with mpfr_t

	// 1
	// 2
	bool operator== (const int) const;
	bool operator!= (const int) const;
	bool operator< (const int) const;
	bool operator> (const int) const;
	bool operator<= (const int) const;
	bool operator>= (const int) const;
	// this template classes should work at least with Gmpz and Gmpq
	template <class T> bool operator== (const T &) const;
	template <class T> bool operator!= (const T &) const;
	bool operator< (const CGAL::Gmpz &) const;
	bool operator< (const CGAL::Gmpq &) const;
	bool operator> (const CGAL::Gmpz &) const;
	bool operator> (const CGAL::Gmpq &) const;
	template <class T> bool operator<= (const T &) const;
	template <class T> bool operator>= (const T &) const;
	// 3
	Algebraic_1 operator- () const;
	Algebraic_1 operator+ (const Algebraic_1 &) const;
	inline Algebraic_1 operator- (const Algebraic_1 &) const;
	Algebraic_1 operator* (const Algebraic_1 &) const;
	Algebraic_1& operator+= (const Algebraic_1 &);
	Algebraic_1& operator-= (const Algebraic_1 &);
	Algebraic_1& operator*= (const Algebraic_1 &);
	// 4
	//--------------------------------------------------
	// Algebraic_1 operator+ (const int) const;
	// Algebraic_1 operator- (const int) const;
	// Algebraic_1 operator* (const int) const;
	//-------------------------------------------------- 
	Algebraic_1& operator+= (const int);
	Algebraic_1& operator-= (const int);
	Algebraic_1& operator*= (const int);
	//--------------------------------------------------
	// Algebraic_1 operator+ (const CGAL::Gmpz &) const;
	// Algebraic_1 operator- (const CGAL::Gmpz &) const;
	// Algebraic_1 operator* (const CGAL::Gmpz &) const;
	//-------------------------------------------------- 
	Algebraic_1& operator+= (const CGAL::Gmpz &);
	Algebraic_1& operator-= (const CGAL::Gmpz &);
	Algebraic_1& operator*= (const CGAL::Gmpz &);
	//--------------------------------------------------
	// Algebraic_1 operator+ (const CGAL::Gmpq &) const;
	// Algebraic_1 operator- (const CGAL::Gmpq &) const;
	// Algebraic_1 operator* (const CGAL::Gmpq &) const;
	//-------------------------------------------------- 
	Algebraic_1& operator+= (const CGAL::Gmpq &);
	Algebraic_1& operator-= (const CGAL::Gmpq &);
	Algebraic_1& operator*= (const CGAL::Gmpq &);
	// 5: the required functions are outside the class 
	// 6
	bool is_valid () const;
	bool is_finite () const;
	double to_double () const;
	std::pair <double,double> to_interval () const;
	// 7
	//--------------------------------------------------
	// Algebraic_1 operator/ (const int) const;
	//-------------------------------------------------- 
	Algebraic_1 operator/ (const Algebraic_1 &) const;
	Algebraic_1& operator/= (const Algebraic_1 &);
	Algebraic_1& operator/= (const int);
	// 8
	Algebraic_1 sqrt () const;
	// 9
	std::ostream& show (std::ostream &);
	// 10
	// (the other comparison cases for mp[zq]_t should be covered by the
	// template functions)
	bool operator< (const mpz_t &) const;
	bool operator> (const mpz_t &) const;
	bool operator< (const mpq_t &) const;
	bool operator> (const mpq_t &) const;
	//--------------------------------------------------
	// Algebraic_1 operator+ (const mpz_t &) const;
	// Algebraic_1 operator- (const mpz_t &) const;
	// Algebraic_1 operator* (const mpz_t &) const;
	//-------------------------------------------------- 
	Algebraic_1& operator+= (const mpz_t &);
	Algebraic_1& operator-= (const mpz_t &);
	Algebraic_1& operator*= (const mpz_t &);
	//--------------------------------------------------
	// Algebraic_1 operator+ (const mpq_t &) const;
	// Algebraic_1 operator- (const mpq_t &) const;
	// Algebraic_1 operator* (const mpq_t &) const;
	//-------------------------------------------------- 
	Algebraic_1& operator+= (const mpq_t &);
	Algebraic_1& operator-= (const mpq_t &);
	Algebraic_1& operator*= (const mpq_t &);

	// 11
	Algebraic_1 (const mpfr_t &);	// constructor I
	Algebraic_1 (const mpfr_t &, const mpfr_t &);	// constructor II
	Algebraic_1& operator= (const mpfr_t &);	// assigning
	// comparison (previous template definitions should work with mpfr_t)
	bool operator< (const mpfr_t &) const;
	bool operator> (const mpfr_t &) const;
	// arithmetics
	//--------------------------------------------------
	// Algebraic_1 operator+ (const mpfr_t &) const;
	// Algebraic_1 operator- (const mpfr_t &) const;
	// Algebraic_1 operator* (const mpfr_t &) const;
	// Algebraic_1 operator/ (const mpfr_t &) const;
	//-------------------------------------------------- 
	Algebraic_1& operator+= (const mpfr_t &);
	Algebraic_1& operator-= (const mpfr_t &);
	Algebraic_1& operator*= (const mpfr_t &);
	Algebraic_1& operator/= (const mpfr_t &);
};

std::ostream& operator<< (std::ostream &, Algebraic_1 &);

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
// the implementation of functions inside and outside the class
#include <CGAL/Gbrs_algebraic_1.C>
#endif	// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif	// CGAL_MPFIINTERVAL_H
