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

// TODO:
//	-change the order in which functions are written in the file, so it
//	becomes more readable (someone, someday, will do this)
//	-think about precision propagation in arithmetic functions (MPFI is
//	supposed to do this)
//	-avoid the use of BOOST, because the operators that it provides here
//	are not correct (note that alg += non_alg will need to clear its
//	pointer to the polynomial, while alg + non_alg won't)

// NOTE:
// some functions are not coded because BOOST will provide them; they are
// commented in the code

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gbrs_polynomial_1.h>
#include <CGAL/Gbrs_algebraic_1.h>
#include <iostream>
#include <mpfr.h>
#include <mpfi.h>
#include <mpfi_io.h>

CGAL_BEGIN_NAMESPACE

// the exception object
comparison_overlap_exn exn_overlap;

// what to do when a comparison fails?
void overlap () {
	throw exn_overlap;
}

// ////////////////////
// Algebraic_1_rep member functions
// ////////////////////
// copy constructor
Algebraic_1_rep::Algebraic_1_rep(const Algebraic_1_rep &a){
	mpfi_set(mpfI,a.mpfI);
	poly=a.poly;
	nr=a.nr;
	mult=a.mult;
	rsprec=a.rsprec;
};
// copy assignement operator
Algebraic_1_rep& Algebraic_1_rep::operator=(const Algebraic_1_rep &a){
	mpfi_set(mpfI,a.mpfI);
	poly=a.poly;
	nr=a.nr;
	mult=a.mult;
	rsprec=a.rsprec;
	return *this;
};

// ////////////////////
// Algebraic_1 member functions
// ////////////////////

// XXX: are copy constructor and copy assignement operator really needed?

// copy constructor
Algebraic_1::Algebraic_1(const Algebraic_1 &i){
	//set_mpfi(mpfi(),i.mpfi());	// this copies the mpfi
	set_mpfi_ptr(i.mpfi());	// this copies the pointer to the mpfi
	set_pol(i.pol());
	set_nr(i.nr());
	set_mult(i.mult());
	set_rsprec(i.rsprec());
};

// copy assignement operator
Algebraic_1& Algebraic_1::operator=(const Algebraic_1 &i){
	//mpfi_set(mpfi(),i.mpfi());
	set_mpfi(i.mpfi());
	set_pol(i.pol());
	set_nr(i.nr());
	set_mult(i.mult());
	set_rsprec(i.rsprec());
	return *this;
};

// constructors of a "point" interval
Algebraic_1::Algebraic_1 () {};

Algebraic_1::Algebraic_1(int i){
	mpq_t temp;
	mpfi_set_si(mpfi(),(long int)i);
	mpq_init(temp);
	mpq_set_si(temp,(long int)i,1);
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	set_nr(0);
};

Algebraic_1::Algebraic_1(unsigned int i){
	mpq_t temp;
	mpfi_set_ui(mpfi(),i);
	mpq_init(temp);
	mpq_set_ui(temp,(unsigned int)i,1);
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	set_nr(0);
};

Algebraic_1::Algebraic_1(long int i){
	mpq_t temp;
	mpfi_set_si(mpfi(),i);
	mpq_init(temp);
	mpq_set_si(temp,i,1);
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	set_nr(0);
};

Algebraic_1::Algebraic_1(unsigned long int i){
	mpq_t temp;
	mpfi_set_ui(mpfi(),i);
	mpq_init(temp);
	mpq_set_ui(temp,i,1);
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	set_nr(0);
};

Algebraic_1::Algebraic_1(double d){
	mpq_t temp;
	mpfi_set_d(mpfi(),d);
	mpq_init(temp);
	mpq_set_d(temp,d);
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	set_nr(0);
};

Algebraic_1::Algebraic_1(mpz_srcptr z){
	mpq_t temp;
	mpfi_set_z(mpfi(),z);
	mpq_init(temp);
	mpq_set_z(temp,z);
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	set_nr(0);
};

Algebraic_1::Algebraic_1(mpq_srcptr q){
	mpfi_set_q(mpfi(),q);
	Rational_polynomial_1 *p=new Rational_polynomial_1(q);
	set_pol(*p);
	set_nr(0);
};

Algebraic_1::Algebraic_1(const CGAL::Gmpz &z){
	mpq_t temp;
	mpfi_set_z(mpfi(),z.mpz());
	mpq_init(temp);
	mpq_set_z(temp,z.mpz());
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	set_nr(0);
};

Algebraic_1::Algebraic_1(const CGAL::Gmpq &q){
	mpfi_set_q(mpfi(),q.mpq());
	Rational_polynomial_1 *p=new Rational_polynomial_1(q.mpq());
	set_pol(*p);
	set_nr(0);
};

// constructors of a "proper" interval (these aren't interesting at all)
Algebraic_1::Algebraic_1 (int l, int r) {
	mpfi_interv_si (mpfi (), (long int)l, (long int)r);
};

Algebraic_1::Algebraic_1 (unsigned int l, unsigned int r) {
	mpfi_interv_ui (mpfi (), (unsigned long int)l, (unsigned long int)r);
};

Algebraic_1::Algebraic_1 (long int l, long int r) {
	mpfi_interv_si (mpfi (), l, r);
};

Algebraic_1::Algebraic_1 (unsigned long int l, unsigned long int r) {
	mpfi_interv_ui (mpfi (), l, r);
};

Algebraic_1::Algebraic_1 (double l, double r) {
	mpfi_interv_d (mpfi (), l, r);
};

Algebraic_1::Algebraic_1(mpz_srcptr l,mpz_srcptr r){
	mpfi_interv_z (mpfi (), l, r);
};

Algebraic_1::Algebraic_1(mpq_srcptr l,mpq_srcptr r){
	mpfi_interv_q (mpfi (), l, r);
};

Algebraic_1::Algebraic_1 (const CGAL::Gmpz &l, const CGAL::Gmpz &r) {
	mpfi_interv_z (mpfi (), l.mpz(), r.mpz());
};

Algebraic_1::Algebraic_1 (const CGAL::Gmpq &l, const CGAL::Gmpq &r) {
	mpfi_interv_q (mpfi (), l.mpq(), r.mpq());
};

inline Algebraic_1::Algebraic_1 (mpfi_srcptr i) {set_mpfi(i);};

// interesting constructor
Algebraic_1::Algebraic_1(mpfi_srcptr i,Rational_polynomial_1 &p,
		const int n,const int m,const int rsp){
	set_mpfi_ptr(i);
	set_pol(p);
	set_nr(n);
	set_mult(m);
	set_rsprec(rsp);
	//p.set_root(*this);
};

// destructor
/* not needed
Algebraic_1::~Algebraic_1 () {};
*/

void Algebraic_1::get_endpoints(mpfr_ptr l,mpfr_ptr r)const{
	mpfi_get_left(l,mpfi());
	mpfi_get_right(r,mpfi());
};

bool Algebraic_1::overlaps(const Algebraic_1&a)const{
	if(mpfr_lessequal_p(left(),a.left()))
		return mpfr_lessequal_p(a.left(),right());
	else
		return mpfr_lessequal_p(left(),a.right());
};

// assignment
Algebraic_1& Algebraic_1::operator=(const long int i){
	mpq_t temp;
	mpfi_set_si(mpfi(),(long int)i);
	mpq_init(temp);
	mpq_set_si(temp,(long int)i,1);
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	return *this;
};

Algebraic_1& Algebraic_1::operator=(mpz_srcptr z){
	mpq_t temp;
	mpfi_set_z(mpfi(),z);
	mpq_init(temp);
	mpq_set_z(temp,z);
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	return *this;
};

Algebraic_1& Algebraic_1::operator=(mpq_srcptr q){
	mpfi_set_q(mpfi(),q);
	Rational_polynomial_1 *p=new Rational_polynomial_1(q);
	set_pol(*p);
	return *this;
};

Algebraic_1& Algebraic_1::operator=(const CGAL::Gmpz &z){
	mpq_t temp;
	mpfi_set_z(mpfi(),z.mpz());
	mpq_init(temp);
	mpq_set_z(temp,z.mpz());
	Rational_polynomial_1 *p=new Rational_polynomial_1(temp);
	mpq_clear(temp);
	set_pol(*p);
	return *this;
};

Algebraic_1& Algebraic_1::operator=(const CGAL::Gmpq &q){
	mpfi_set_q(mpfi(),q.mpq());
	Rational_polynomial_1 *p=new Rational_polynomial_1(q.mpq());
	set_pol(*p);
	return *this;
};

// 1
// 2
// comparisons with ints
bool Algebraic_1::operator== (const int n2) const {
	if (contains (n2))
		if (is_point ())
			return true;
		else
			overlap ();
	return false;
};

bool Algebraic_1::operator!= (const int n2) const {
	return !(operator== (n2));
};

bool Algebraic_1::operator< (const int n2) const {
	if (contains (n2))
		overlap();
	if(mpfr_cmp_si(right(),n2)<0)
		return true;
	return false;
};

bool Algebraic_1::operator> (const int n2) const {
	if (contains (n2))
		overlap();
	if(mpfr_cmp_si(left(),n2)>0)
		return true;
	return false;
};

bool Algebraic_1::operator<= (const int n2) const {
	return ((operator== (n2)) || (operator< (n2)));
};

bool Algebraic_1::operator>= (const int n2) const {
	return ((operator== (n2)) || (operator> (n2)));
};

// comparisons with Gmpz and Gmpq
bool Algebraic_1::operator< (const CGAL::Gmpz &n2) const {
	if (contains (n2))
		overlap();
	if(mpfr_cmp_z(right(),n2.mpz())<0)
		return true;
	return false;
};

bool Algebraic_1::operator< (const CGAL::Gmpq &n2) const {
	if (contains (n2))
		overlap();
	if(mpfr_cmp_q(right(),n2.mpq())<0)
		return true;
	return false;
};

bool Algebraic_1::operator> (const CGAL::Gmpz &n2) const {
	if (contains (n2))
		overlap();
	if(mpfr_cmp_z(left(),n2.mpz())>0)
		return true;
	return false;
};

bool Algebraic_1::operator> (const CGAL::Gmpq &n2) const {
	if (contains (n2))
		overlap();
	if(mpfr_cmp_q(left(),n2.mpq())>0)
		return true;
	return false;
};

// 3

Algebraic_1 Algebraic_1::operator- () const {
	mpfi_t n;
	mpfi_init (n);
	mpfi_neg (n, mpfi ());
	Algebraic_1 ret (n);
	return ret;
};

Algebraic_1 Algebraic_1::operator+ (const Algebraic_1 &n2) const {
	mpfi_t n;
	mpfi_init (n);
	mpfi_add (n, mpfi (), n2.mpfi());
	Algebraic_1 ret (n);
	return ret;
};

Algebraic_1 Algebraic_1::operator- (const Algebraic_1 &n2) const {
	return (*this + (-n2));
};

Algebraic_1 Algebraic_1::operator* (const Algebraic_1 &n2) const {
	mpfi_t n;
	mpfi_init (n);
	mpfi_mul (n, mpfi (), n2.mpfi());
	Algebraic_1 ret (n);
	return ret;
};

Algebraic_1& Algebraic_1::operator+= (const Algebraic_1 &n2) {
	mpfi_add (mpfi (), mpfi (), n2.mpfi ());
	clear_pol ();
	return *this;
};

Algebraic_1& Algebraic_1::operator-= (const Algebraic_1 &n2) {
	mpfi_sub (mpfi (), mpfi (), n2.mpfi ());
	clear_pol ();
	return *this;
};

Algebraic_1& Algebraic_1::operator*= (const Algebraic_1 &n2) {
	mpfi_mul (mpfi (), mpfi (), n2.mpfi ());
	clear_pol ();
	return *this;
};

// 4
// this (op) int
//--------------------------------------------------
// BOOST:
// Algebraic_1 Algebraic_1::operator+ (const int n2) const
// Algebraic_1 Algebraic_1::operator- (const int n2) const
// Algebraic_1 Algebraic_1::operator* (const int n2) const
//-------------------------------------------------- 

Algebraic_1& Algebraic_1::operator+= (const int n2) {
	mpfi_add_si (mpfi (), mpfi (), (long int)n2);
	return *this;
};

Algebraic_1& Algebraic_1::operator-= (const int n2) {
	mpfi_sub_si (mpfi (), mpfi (), (long int)n2);
	return *this;
};

Algebraic_1& Algebraic_1::operator*= (const int n2) {
	mpfi_mul_si (mpfi (), mpfi (), (long int)n2);
	return *this;
};

// this (op) Gmpz
//--------------------------------------------------
// BOOST:
// Algebraic_1 Algebraic_1::operator+ (const CGAL::Gmpz &n2) const
// Algebraic_1 Algebraic_1::operator- (const CGAL::Gmpz &n2) const
// Algebraic_1 Algebraic_1::operator* (const CGAL::Gmpz &n2) const
//-------------------------------------------------- 

Algebraic_1& Algebraic_1::operator+= (const CGAL::Gmpz &n2) {
	mpfi_add_z (mpfi (), mpfi (), n2.mpz());
	return *this;
};

Algebraic_1& Algebraic_1::operator-= (const CGAL::Gmpz &n2) {
	mpfi_sub_z (mpfi (), mpfi (), n2.mpz());
	return *this;
};

// this (op) Gmpq
//--------------------------------------------------
// BOOST:
// Algebraic_1 Algebraic_1::operator+ (const CGAL::Gmpq &n2) const
// Algebraic_1 Algebraic_1::operator- (const CGAL::Gmpq &n2) const
// Algebraic_1 Algebraic_1::operator* (const CGAL::Gmpq &n2) const
//-------------------------------------------------- 

Algebraic_1& Algebraic_1::operator+= (const CGAL::Gmpq &n2) {
	mpfi_add_q (mpfi (), mpfi (), n2.mpq());
	return *this;
};

Algebraic_1& Algebraic_1::operator-= (const CGAL::Gmpq &n2) {
	mpfi_sub_q (mpfi (), mpfi (), n2.mpq());
	return *this;
};

Algebraic_1& Algebraic_1::operator*= (const CGAL::Gmpq &n2) {
	mpfi_mul_q (mpfi (), mpfi (), n2.mpq());
	return *this;
};

// 5

// 6
std::pair <double, double> Algebraic_1::to_interval () const {
	return std::make_pair(
			mpfr_get_d(&(mpfi()->left),GMP_RNDN),
			mpfr_get_d(&(mpfi()->right),GMP_RNDN));
};

// 7
//--------------------------------------------------
// BOOST:
// Algebraic_1 Algebraic_1::operator/ (const int n2) const
//-------------------------------------------------- 

Algebraic_1 Algebraic_1::operator/ (const Algebraic_1 &n2) const {
	mpfi_t n;
	mpfi_init (n);
	mpfi_div (n, mpfi (), n2.mpfi());
	Algebraic_1 ret (n);
	return ret;
};

Algebraic_1& Algebraic_1::operator/= (const Algebraic_1 &n2) {
	mpfi_div (mpfi (), mpfi (), n2.mpfi ());
	clear_pol ();
	return *this;
};

Algebraic_1& Algebraic_1::operator/= (const int n2) {
	mpfi_div_si (mpfi (), mpfi (), (long int)n2);
	return *this;
};

// 8
Algebraic_1 Algebraic_1::sqrt () const {
	mpfi_t s;
	mpfi_init (s);
	mpfi_sqrt (s, mpfi ());
	Algebraic_1 ret (s);
	return ret;
};

// 9
std::ostream& Algebraic_1::show(std::ostream &o,int digits)const{
	if(is_point())
		return(o<<mpfi_get_d(mpfi()));
	return(o<<"["<<mpfr_get_d(left(),GMP_RNDN)<<","<<
			mpfr_get_d(right(),GMP_RNDN)<<"]");

	// this is to exactly display the interval
	/*char *str1, *str2;
	mpfr_t op1, op2;
	mp_exp_t *expptr1, *expptr2;

	expptr1 = (mp_exp_t*)malloc(sizeof(mp_exp_t));
	expptr2 = (mp_exp_t*)malloc(sizeof(mp_exp_t));

	mpfr_inits (op1, op2, NULL);

	mpfi_get_left (op1, mpfi ());
	mpfi_get_right (op2, mpfi ());

	str1=mpfr_get_str(NULL,expptr1,10,digits,op1,GMP_RNDN);
	str2=mpfr_get_str(NULL,expptr2,10,digits,op2,GMP_RNDN);

	if (str1[0] == '-')
		o << "[-." << str1+sizeof(char) << "e" << *expptr1;
	else
		o << "[." << str1 << "e" << *expptr1;
	if (str2[0] == '-')
		o << ",-." << str2+sizeof(char) << "e" << *expptr2 << "]";
	else
		o << ",." << str2 << "e" << *expptr2 << "]";

	mpfr_free_str (str1);
	mpfr_free_str (str2);
	mpfr_clears (op1, op2, NULL);

	return o;*/
};

// 10
bool Algebraic_1::operator<(mpz_srcptr n2)const{
	if (contains (n2))
		overlap();
	if(mpfr_cmp_z(right(),n2)<0)
		return true;
	return false;
};

bool Algebraic_1::operator>(mpz_srcptr n2)const{
	if (contains (n2))
		overlap();
	if(mpfr_cmp_z(left(),n2)>0)
		return true;
	return false;
};

bool Algebraic_1::operator<(mpq_srcptr n2)const{
	if (contains (n2))
		overlap();
	if(mpfr_cmp_q(right(),n2)<0)
		return true;
	return false;
};

bool Algebraic_1::operator>(mpq_srcptr n2)const{
	if (contains (n2))
		overlap();
	if(mpfr_cmp_q(left(),n2)>0)
		return true;
	return false;
};

//--------------------------------------------------
// BOOST:
// Algebraic_1 Algebraic_1::operator+(mpz_srcptr n2)const
// Algebraic_1 Algebraic_1::operator-(mpz_srcptr n2)const
// Algebraic_1 Algebraic_1::operator*(mpz_srcptr n2)const
//-------------------------------------------------- 

Algebraic_1& Algebraic_1::operator+=(mpz_srcptr n2){
	mpfi_add_z (mpfi (), mpfi (), n2);
	return *this;
};

Algebraic_1& Algebraic_1::operator-=(mpz_srcptr n2){
	mpfi_sub_z (mpfi (), mpfi (), n2);
	return *this;
};

Algebraic_1& Algebraic_1::operator*=(mpz_srcptr n2){
	mpfi_mul_z(mpfi(),mpfi(),n2);
	return *this;
};

// this (op) mpq_t
//--------------------------------------------------
// BOOST:
// Algebraic_1 Algebraic_1::operator+(mpq_srcptr n2)const
// Algebraic_1 Algebraic_1::operator-(mpq_srcptr n2)const
// Algebraic_1 Algebraic_1::operator*(mpq_srcptr n2)const
//-------------------------------------------------- 

Algebraic_1& Algebraic_1::operator+=(mpq_srcptr n2){
	mpfi_add_q (mpfi (), mpfi (), n2);
	return *this;
};

Algebraic_1& Algebraic_1::operator-=(mpq_srcptr n2){
	mpfi_sub_q (mpfi (), mpfi (), n2);
	return *this;
};

Algebraic_1& Algebraic_1::operator*=(mpq_srcptr n2){
	mpfi_mul_q (mpfi (), mpfi (), n2);
	return *this;
};

// 11. all the functions with mpfr_t that need to be inside the class

// constructor I
// we can't convert an mpfr_t to a rational easily, so calling this
// constructor isn't a good idea
Algebraic_1::Algebraic_1(mpfr_srcptr r){
	mpfi_set_fr (mpfi (), r);
};

// constructor II
Algebraic_1::Algebraic_1(mpfr_srcptr l,mpfr_srcptr r){
	mpfi_interv_fr (mpfi (), l, r);
};

// assigning: mpfi = mpfr
// big problem: the same that the constructor I above
Algebraic_1& Algebraic_1::operator=(mpfr_srcptr r){
	mpfi_set_fr (mpfi (), r);
	return *this;
};

// comparison: mpfi (op) mpfr
//	NOTE: the previous template definitions of operators =, !=, >= and <=
//	should work with mpfr_t
bool Algebraic_1::operator<(mpfr_srcptr n2)const{
	if (contains (n2))
		overlap();
	return(mpfr_less_p(right(),n2));
};

bool Algebraic_1::operator>(mpfr_srcptr n2)const{
	if (contains (n2))
		overlap();
	return(mpfr_greater_p(left(),n2));
};

// arithmetics: mpfi (op) mpfr
//--------------------------------------------------
// BOOST:
// Algebraic_1 Algebraic_1::operator+(mpfr_srcptr f)const
// Algebraic_1 Algebraic_1::operator-(mpfr_srcptr f)const
// Algebraic_1 Algebraic_1::operator*(mpfr_srcptr f)const
// Algebraic_1 Algebraic_1::operator/(mpfr_srcptr f)const
//-------------------------------------------------- 

Algebraic_1& Algebraic_1::operator+=(mpfr_srcptr f){
	mpfi_add_fr (mpfi (), mpfi (), f);
	return *this;
};

Algebraic_1& Algebraic_1::operator-=(mpfr_srcptr f){
	mpfi_sub_fr (mpfi (), mpfi (), f);
	return *this;
};

Algebraic_1& Algebraic_1::operator*=(mpfr_srcptr f){
	mpfi_mul_fr (mpfi (), mpfi (), f);
	return *this;
};

Algebraic_1& Algebraic_1::operator/=(mpfr_srcptr f){
	mpfi_div_fr (mpfi (), mpfi (), f);
	return *this;
};




// ----------------------------------------
// end of the Algebraic_1 class
// ----------------------------------------



// These functions are required, but they need to be coded outside the class:

// 1.5
bool operator== (const Algebraic_1 &n1, const Algebraic_1 &n2) {
	if(n1.is_point()&&n2.is_point()&&
			(!mpfr_cmp(n1.left(),n2.left())))
		return true;
	if((mpfr_less_p(n1.right(),n2.left()))
			||(mpfr_less_p(n2.right(),n1.left())))
		return false;
	// if intervals are not points and they have the same bounds, they are
	// not necessarily equal
	overlap();
	return false;	// never reached
};

bool operator<(const Algebraic_1 &n1,const Algebraic_1 &n2){
	if(mpfr_less_p(n1.right(),n2.left()))
		return true;
	if(mpfr_less_p(n2.right(),n1.left()))
		return false;
	overlap();
	return false;	// this never occurs
};

// 2.5
// comparison between int|mpfr_t|mp[zq]_t|Gmp[zq] and intervals
// 5.5

// 7.5
// anything / interval
//--------------------------------------------------
// BOOST:
// Algebraic_1 operator/ (const int n1, const Algebraic_1 &n2)
// Algebraic_1 operator/ (const CGAL::Gmpz &n1, const Algebraic_1 &n2)
// Algebraic_1 operator/ (const CGAL::Gmpq &n1, const Algebraic_1 &n2)
//-------------------------------------------------- 

// 8.5
// 9.5
// 11.5
// all the mpfr functions that can't be inside the class (and aren't covered
// by the template functions

// arithmetics
//--------------------------------------------------
// BOOST:
// template <class T>
// Algebraic_1 operator+ (const T &n1, const Algebraic_1 &n2)
// template <class T>
// Algebraic_1 operator- (const T &n1, const Algebraic_1 &n2)
// template <class T>
// Algebraic_1 operator* (const T &n1, const Algebraic_1 &n2)
// Algebraic_1 operator/ (mpfr_srcptr n1, const Algebraic_1 &n2)
//-------------------------------------------------- 

// not implemented: mpfr (op=) mpfi (because they must not return an interval)
// XXX: should them be implemented?
// XXX: will BOOST implement them?

CGAL_END_NAMESPACE
