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

	Algebraic_1_rep():poly(NULL),nr(-1),mult(-1),rsprec(0){mpfi_init(mpfI);}
	~Algebraic_1_rep () {}

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
	Algebraic_1 (mpfi_t &);
	Algebraic_1 (const Algebraic_1 &);

	// the only interesting constructor
	Algebraic_1 (const mpfi_ptr &, const Rational_polynomial_1 &,
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
	const mpfi_t & mpfi () const;
	mpfi_t & mpfi ();
	const Rational_polynomial_1 & pol () const;
	Rational_polynomial_1 & pol ();
	const int nr () const;
	const int mult () const;
	const int rsprec () const;
	void clear_pol ();
	void set_pol (const Rational_polynomial_1 &);
	void set_nr (const int);
	void set_mult (const int);
	void set_rsprec (const int);
	void set_prec (mp_prec_t);
	mp_prec_t get_prec ();
	void get_left (mpfr_t &) const;
	void get_right (mpfr_t &) const;
	void get_endpoints (mpfr_t &, mpfr_t &) const;
	bool is_consistent () const;
	bool is_point () const;	// are endpoints equal?
	bool contains (const int n) const;
	bool contains (const mpfr_t &n) const;
	bool contains (const mpz_t &n) const;
	bool contains (const mpq_t &n) const;
	bool contains (const Gmpz &n) const;
	bool contains (const Gmpq &n) const;

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
	Algebraic_1 operator- (const Algebraic_1 &) const;
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

// //////////////////////////
// inline functions of the class

inline const mpfi_t & Algebraic_1::mpfi() const{
	return Ptr()->mpfI;
};

inline mpfi_t & Algebraic_1::mpfi(){
	return ptr()->mpfI;
};

inline const Rational_polynomial_1 & Algebraic_1::pol() const{
	return *(Ptr()->poly);
};

inline Rational_polynomial_1 & Algebraic_1::pol(){
	return *(ptr()->poly);
};

inline const int Algebraic_1::nr() const{
	return ptr()->nr;
};

inline const int Algebraic_1::mult() const{
	return ptr()->mult;
};

inline const int Algebraic_1::rsprec() const{
	return ptr()->rsprec;
};

inline void Algebraic_1::set_prec(mp_prec_t p){
	mpfi_round_prec(mpfi(),p);
};

inline mp_prec_t Algebraic_1::get_prec() { return mpfi_get_prec(mpfi()); };

inline void Algebraic_1::get_left(mpfr_t &f) const{
	mpfi_get_left(f,mpfi());
}

inline void Algebraic_1::get_right(mpfr_t &f) const{
	mpfi_get_right(f,mpfi());
}

inline void Algebraic_1::get_endpoints(mpfr_t &l,mpfr_t &r) const{
	mpfi_get_left(l,mpfi());
	mpfi_get_right(r,mpfi());
}

inline bool Algebraic_1::is_consistent() const{
	return (ptr()->poly);
};

inline bool Algebraic_1::is_point() const{
	mpfr_t l,r;
	mpfr_inits(l,r,NULL);
	get_endpoints(l,r);
	int comp=mpfr_equal_p(l,r);
	mpfr_clears(l,r,NULL);
	return (comp!=0);
}

inline bool Algebraic_1::contains(const int n) const{
	mpfr_t end;
	mpfr_init(end);
	int comp;
	get_left(end);	// first, we compare the left end
	comp=mpfr_cmp_si(end,n);
	if (comp>0){	// n is lower than the left end
		mpfr_clear(end);
		return false;
	}
	get_right(end);	// now, the right one
	comp=mpfr_cmp_si(end,n);
	if (comp<0) {	// n is higher than the right end
		mpfr_clear(end);
		return false;
	}
	return true;
}

inline bool Algebraic_1::contains(const mpfr_t &n) const{
	mpfr_t end;
	mpfr_init(end);
	int comp;
	get_left(end);	// first, we compare the left end
	comp=mpfr_cmp(end,n);
	if (comp>0) {	// n is lower than the left end
		mpfr_clear(end);
		return false;
	}
	get_right(end);	// now, the right one
	comp=mpfr_cmp(end,n);
	if (comp<0) {	// n is higher than the right end
		mpfr_clear(end);
		return false;
	}
	return true;
}

inline bool Algebraic_1::contains(const mpz_t &n) const{
	mpfr_t end;
	mpfr_init(end);
	int comp;
	get_left(end);	// first, we compare the left end
	comp=mpfr_cmp_z(end,n);
	if (comp>0) {	// n is lower than the left end
		mpfr_clear(end);
		return false;
	}
	get_right(end);	// now, the right one
	comp=mpfr_cmp_z(end,n);
	if (comp<0) {	// n is higher than the right end
		mpfr_clear(end);
		return false;
	}
	return true;
}

inline bool Algebraic_1::contains(const mpq_t &n) const{
	mpfr_t end;
	mpfr_init(end);
	int comp;
	get_left(end);	// first, we compare the left end
	comp=mpfr_cmp_q(end,n);
	if (comp>0) {	// n is lower than the left end
		mpfr_clear(end);
		return false;
	}
	get_right(end);	// now, the right one
	comp=mpfr_cmp_q(end,n);
	if (comp<0) {	// n is higher than the right end
		mpfr_clear(end);
		return false;
	}
	return true;
}

inline bool Algebraic_1::contains (const Gmpz &n) const{
	mpfr_t end;
	mpfr_init(end);
	int comp;
	get_left(end);	// first, we compare the left end
	comp=mpfr_cmp_z(end,n.mpz());
	if (comp>0) {	// n is lower than the left end
		mpfr_clear(end);
		return false;
	}
	get_right(end);	// now, the right one
	comp=mpfr_cmp_z(end,n.mpz());
	if (comp<0) {	// n is higher than the right end
		mpfr_clear(end);
		return false;
	}
	return true;
}

inline bool Algebraic_1::contains (const Gmpq &n) const {
	mpfr_t end;
	mpfr_init(end);
	int comp;
	get_left(end);	// first, we compare the left end
	comp=mpfr_cmp_q(end,n.mpq());
	if (comp>0) {	// n is lower than the left end
		mpfr_clear(end);
		return false;
	}
	get_right(end);	// now, the right one
	comp=mpfr_cmp_q(end,n.mpq());
	if (comp<0) {	// n is higher than the right end
		mpfr_clear(end);
		return false;
	}
	return true;
}

// //////////////////////////
// other functions coded not as class members
std::ostream& operator<<(std::ostream&,Algebraic_1&);

bool operator==(const Algebraic_1&,const Algebraic_1&);
bool operator!=(const Algebraic_1&,const Algebraic_1&);
bool operator<(const Algebraic_1&,const Algebraic_1&);
bool operator>(const Algebraic_1&,const Algebraic_1 &n2);
bool operator<=(const Algebraic_1&,const Algebraic_1&);
bool operator>=(const Algebraic_1&,const Algebraic_1&);
template <class T> bool operator==(const T&,const Algebraic_1&);
template <class T> bool operator!=(const T&,const Algebraic_1&);
template <class T> bool operator<(const T&,const Algebraic_1&);
template <class T> bool operator>(const T&,const Algebraic_1&);
template <class T> bool operator<=(const T&,const Algebraic_1&);
template <class T> bool operator>=(const T&,const Algebraic_1&);

bool is_valid(const Algebraic_1&);
bool is_finite(const Algebraic_1&);
double to_double(const Algebraic_1&);
std::pair<double, double> to_interval(const Algebraic_1&);

Algebraic_1 sqrt(const Algebraic_1&);

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Gbrs_algebraic_1_impl.h>
#endif	// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif	// CGAL_MPFIINTERVAL_H
