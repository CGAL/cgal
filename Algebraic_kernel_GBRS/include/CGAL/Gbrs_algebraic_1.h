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

#ifndef CGAL_GBRS_ALGEBRAIC_1_H
#define CGAL_GBRS_ALGEBRAIC_1_H

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
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

class Rational_polynomial_1;

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
	Algebraic_1_rep (const Algebraic_1_rep &);
	Algebraic_1_rep & operator= (const Algebraic_1_rep &);
};

// The class of the MPFI intervals. It's a model of the RingNumberType concept
class Algebraic_1
: Handle_for<Algebraic_1_rep>,
	boost::field_operators2<Algebraic_1, int,
	boost::field_operators2<Algebraic_1, mpz_t,
	boost::field_operators2<Algebraic_1, mpq_t,
	boost::field_operators2<Algebraic_1, Gmpz,
	boost::field_operators2<Algebraic_1, Gmpq,
	boost::field_operators2<Algebraic_1, mpfr_t > > > > > > {
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

	// copy constructor and copy assignement operator
	Algebraic_1(const Algebraic_1&);
	Algebraic_1& operator=(const Algebraic_1&);

	// constructors I
	// these constructors create an algebraic number which is a point; they
	// will also create a polynomial of which they are the only root
	Algebraic_1 ();
	Algebraic_1 (int);
	Algebraic_1 (unsigned int);
	Algebraic_1 (long int);
	Algebraic_1 (unsigned long int);
	Algebraic_1 (double);
	Algebraic_1(mpz_srcptr);
	Algebraic_1(mpq_srcptr);
	Algebraic_1 (const CGAL::Gmpq &);
	Algebraic_1 (const CGAL::Gmpz &);

	// constructors II
	Algebraic_1 (int, int);
	Algebraic_1 (unsigned int, unsigned int);
	Algebraic_1 (long int, long int);
	Algebraic_1 (unsigned long int, unsigned long int);
	Algebraic_1 (double, double);
	Algebraic_1(mpz_srcptr,mpz_srcptr);
	Algebraic_1(mpq_srcptr,mpq_srcptr);
	Algebraic_1 (const CGAL::Gmpq &, const CGAL::Gmpq &);
	Algebraic_1 (const CGAL::Gmpz &, const CGAL::Gmpz &);
	Algebraic_1 (mpfi_srcptr);

	// the only interesting constructor
	Algebraic_1(mpfi_srcptr,Rational_polynomial_1&,
			const int, const int, const int);

	Algebraic_1& operator= (const long int);
	Algebraic_1& operator=(mpz_srcptr);
	Algebraic_1& operator=(mpq_srcptr);
	Algebraic_1& operator= (const CGAL::Gmpz &);
	Algebraic_1& operator= (const CGAL::Gmpq &);

	// destructor
	/* not needed
	~Algebraic_1 ();
	*/

	// functions related to the member data
	mpfi_srcptr mpfi()const;
	mpfi_ptr mpfi();
	const Rational_polynomial_1 & pol () const;
	//Rational_polynomial_1 & pol ();
	const int nr () const;
	const int mult () const;
	const int rsprec () const;
	void set_mpfi(mpfi_srcptr);
	void set_mpfi_ptr(mpfi_srcptr);
	void clear_pol ();
	void set_pol (const Rational_polynomial_1 &);
	void set_nr (const int);
	void set_mult (const int);
	void set_rsprec (const int);
	void set_prec (mp_prec_t);
	mp_prec_t get_prec()const;
	void get_left(mpfr_ptr)const;
	void get_right(mpfr_ptr)const;
	mpfr_srcptr left()const;
	mpfr_srcptr right()const;
	void get_endpoints(mpfr_ptr,mpfr_ptr)const;
	bool is_consistent () const;
	bool is_point () const;	// are endpoints equal?
	bool overlaps(const Algebraic_1&)const;
	bool contains (const int n) const;
	bool contains(mpfr_srcptr)const;
	bool contains(mpz_srcptr) const;
	bool contains(mpq_srcptr) const;
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
	std::ostream& show(std::ostream&,int=0)const;
	// 10
	// (the other comparison cases for mp[zq]_t should be covered by the
	// template functions)
	bool operator<(mpz_srcptr)const;
	bool operator>(mpz_srcptr)const;
	bool operator<(mpq_srcptr)const;
	bool operator>(mpq_srcptr)const;
	//--------------------------------------------------
	// Algebraic_1 operator+(mpz_srcptr)const;
	// Algebraic_1 operator-(mpz_srcptr)const;
	// Algebraic_1 operator*(mpz_srcptr)const;
	//-------------------------------------------------- 
	Algebraic_1& operator+=(mpz_srcptr);
	Algebraic_1& operator-=(mpz_srcptr);
	Algebraic_1& operator*=(mpz_srcptr);
	//--------------------------------------------------
	// Algebraic_1 operator+(mpq_srcptr)const;
	// Algebraic_1 operator-(mpq_srcptr)const;
	// Algebraic_1 operator*(mpq_srcptr)const;
	//-------------------------------------------------- 
	Algebraic_1& operator+=(mpq_srcptr);
	Algebraic_1& operator-=(mpq_srcptr);
	Algebraic_1& operator*=(mpq_srcptr);

	// 11
	Algebraic_1(mpfr_srcptr);	// constructor I
	Algebraic_1(mpfr_srcptr,mpfr_srcptr);	// constructor II
	Algebraic_1& operator=(mpfr_srcptr);	// assigning
	// comparison (previous template definitions should work with mpfr_t)
	bool operator<(mpfr_srcptr)const;
	bool operator>(mpfr_srcptr)const;
	// arithmetics
	//--------------------------------------------------
	// Algebraic_1 operator+(mpfr_srcptr)const;
	// Algebraic_1 operator-(mpfr_srcptr)const;
	// Algebraic_1 operator*(mpfr_srcptr)const;
	// Algebraic_1 operator/(mpfr_srcptr)const;
	//-------------------------------------------------- 
	Algebraic_1& operator+=(mpfr_srcptr);
	Algebraic_1& operator-=(mpfr_srcptr);
	Algebraic_1& operator*=(mpfr_srcptr);
	Algebraic_1& operator/=(mpfr_srcptr);
};

// //////////////////////////
// inline functions of the class
inline mpfi_srcptr Algebraic_1::mpfi()const{return Ptr()->mpfI;};
inline mpfi_ptr Algebraic_1::mpfi(){return ptr()->mpfI;};
inline const Rational_polynomial_1 & Algebraic_1::pol() const{
	return *(Ptr()->poly);};
//inline Rational_polynomial_1& Algebraic_1::pol(){return *(ptr()->poly);};
inline const int Algebraic_1::nr()const{return ptr()->nr;};
inline const int Algebraic_1::mult()const{return ptr()->mult;};
inline const int Algebraic_1::rsprec()const{
	int p1,p2;
	if((p1=ptr()->rsprec)>(p2=get_prec()))
		return p1;
	else
		return p2;
	};
inline void Algebraic_1::set_mpfi(mpfi_srcptr x){
	mpfi_set(mpfi(),x);};
inline void Algebraic_1::set_mpfi_ptr(mpfi_srcptr x){/**mpfi()=*x;*/mpfi_set(mpfi(),x);};
inline void Algebraic_1::clear_pol(){ptr()->poly=NULL;};
inline void Algebraic_1::set_pol(const Rational_polynomial_1 &p){
	ptr()->poly=const_cast<Rational_polynomial_1*>(&p);}; // thanks Julien!
inline void Algebraic_1::set_nr(const int n){ptr()->nr=n;};
inline void Algebraic_1::set_mult(const int m){ptr()->mult=m;};
inline void Algebraic_1::set_rsprec(const int p){ptr()->rsprec=p;};
inline void Algebraic_1::set_prec(mp_prec_t p){mpfi_round_prec(mpfi(),p);};
inline mp_prec_t Algebraic_1::get_prec()const{return mpfi_get_prec(mpfi());};
inline void Algebraic_1::get_left(mpfr_ptr f)const{mpfi_get_left(f,mpfi());};
inline void Algebraic_1::get_right(mpfr_ptr f)const{mpfi_get_right(f,mpfi());};
inline mpfr_srcptr Algebraic_1::left()const{return &(mpfi()->left);};
inline mpfr_srcptr Algebraic_1::right()const{return &(mpfi()->right);};
inline bool Algebraic_1::is_consistent()const{
	return(&pol()==NULL?false:true);};
inline bool Algebraic_1::is_point()const{
	return(mpfr_equal_p(&(mpfi()->left),&(mpfi()->right)));};
inline bool Algebraic_1::contains(const int n)const{
	return((mpfr_cmp_si(&(mpfi()->left),n)<=0)&&
			(mpfr_cmp_si(&(mpfi()->right),n)>=0));};
inline bool Algebraic_1::contains(mpfr_srcptr n)const{
	return((mpfr_lessequal_p(&(mpfi()->left),n))&&
			(mpfr_greaterequal_p(&(mpfi()->right),n)));};
inline bool Algebraic_1::contains(mpz_srcptr n)const{
	return((mpfr_cmp_z(&(mpfi()->left),n)<=0)&&
			(mpfr_cmp_z(&(mpfi()->right),n)>=0));};
inline bool Algebraic_1::contains(mpq_srcptr n)const{
	return((mpfr_cmp_q(&(mpfi()->left),n)<=0)&&
			(mpfr_cmp_q(&(mpfi()->right),n)>=0));};
inline bool Algebraic_1::contains(const Gmpz &n)const{
	return((mpfr_cmp_z(&(mpfi()->left),n.mpz())<=0)&&
			(mpfr_cmp_z(&(mpfi()->right),n.mpz())>=0));};
inline bool Algebraic_1::contains(const Gmpq &n)const{
	return((mpfr_cmp_q(&(mpfi()->left),n.mpq())<=0)&&
			(mpfr_cmp_q(&(mpfi()->right),n.mpq())>=0));};

// //////////////////////////
// other functions coded not as class members
std::ostream& operator<<(std::ostream&,Algebraic_1&);
std::ostream& operator<<(std::ostream&,const Algebraic_1&);

bool operator==(const Algebraic_1&,const Algebraic_1&);
bool operator!=(const Algebraic_1&,const Algebraic_1&);
bool operator<(const Algebraic_1&,const Algebraic_1&);
bool operator>(const Algebraic_1&,const Algebraic_1 &n2);
bool operator<=(const Algebraic_1&,const Algebraic_1&);
bool operator>=(const Algebraic_1&,const Algebraic_1&);
//template <class T> bool operator==(const T&,const Algebraic_1&);
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

// //////////////////////////
// inline functions outside the class
/*1.5*/inline bool operator!=(const Algebraic_1 &n1,const Algebraic_1 &n2){
	return !(n1==n2);}
inline bool operator>(const Algebraic_1 &n1,const Algebraic_1 &n2){
	return(n2<n1);}
inline bool operator<=(const Algebraic_1 &n1,const Algebraic_1 &n2){
	return((n1==n2)||(n1<n2));}
inline bool operator>=(const Algebraic_1 &n1,const Algebraic_1 &n2){
	return((n1==n2)||(n1>n2));}
/*5.5*/inline bool is_valid(const Algebraic_1 &n){return n.is_valid();};
inline bool is_finite(const Algebraic_1 &n){return n.is_finite();};
inline double to_double(const Algebraic_1 &n){return n.to_double();};
inline std::pair<double,double>to_interval(const Algebraic_1 &n){
	return n.to_interval();};
/*8.5*/inline Algebraic_1 sqrt(const Algebraic_1 &ntval){return ntval.sqrt();};
/*9.5*/inline std::ostream& operator<<(std::ostream &o,Algebraic_1 &n){
	return n.show(o);};
inline std::ostream& operator<<(std::ostream &o,const Algebraic_1 &n){
	return n.show(o);};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Gbrs_algebraic_1_impl.h>
#endif	// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif	// CGAL_GBRS_ALGEBRAIC_1_H
