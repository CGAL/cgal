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
// $URL:  $
// $Id:  $
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

// TODO:
//	-change the order in which functions are written in the file, so it
//	becames more readable
//	-think about precision propagation in arithmetic functions
//	-throw an exception when comparison between intervals is not known
//	(when they have points in common)
//	-add interfaces to more number types

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <iostream>
#include <mpfi.h>
#include <mpfi_io.h>
#include <mpfr.h>

CGAL_BEGIN_NAMESPACE

// constructors
MpfiInterval::MpfiInterval () {};

MpfiInterval::MpfiInterval (int i) {
	mpfi_set_si (mpfi (), (long int)i);
};

MpfiInterval::MpfiInterval (unsigned int i) {
	mpfi_set_ui (mpfi (), i);
};

MpfiInterval::MpfiInterval (long int i) {
	mpfi_set_si (mpfi (), i);
};

MpfiInterval::MpfiInterval (unsigned long int i) {
	mpfi_set_ui (mpfi (), i);
};

MpfiInterval::MpfiInterval (double d) {
	mpfi_set_d (mpfi (), d);
};

MpfiInterval::MpfiInterval (const CGAL::Gmpz &z) {
	mpfi_set_z (mpfi (), z.mpz());
};

MpfiInterval::MpfiInterval (const CGAL::Gmpq &q) {
	mpfi_set_q (mpfi (), q.mpq());
};

MpfiInterval::MpfiInterval (const mpfi_t &i) {
	mpfi_set (mpfi (), i);
};

MpfiInterval::MpfiInterval (const MpfiInterval &i) {
	mpfi_set (mpfi (), i.mpfi ());
};

// destructor
/* not needed
MpfiInterval::~MpfiInterval () {};
*/

inline const mpfi_t & MpfiInterval::mpfi () const { return Ptr()->mpfI; };

inline mpfi_t & MpfiInterval::mpfi () { return ptr()->mpfI; };

inline void MpfiInterval::set_prec (mp_prec_t p) { mpfi_set_prec (mpfi (), p); };

inline mp_prec_t MpfiInterval::get_prec () { return mpfi_get_prec (mpfi ()); };

inline void MpfiInterval::get_left (mpft_r &f) {
	mpfi_get_left (f, mpfi ());
}

inline void MpfiInterval::get_right (mpft_r &f) {
	mpfi_get_right (f, mpfi ());
}

// overcharge for assignment
MpfiInterval MpfiInterval::operator= (const long int i) {
	mpfi_set_si (mpfi (), i);
	return *this;
};

MpfiInterval MpfiInterval::operator= (const CGAL::Gmpz &z) {
	mpfi_set_z (mpfi (), z.mpz());
	return *this;
};

MpfiInterval MpfiInterval::operator= (const CGAL::Gmpq &q) {
	mpfi_set_q (mpfi (), q.mpq());
	return *this;
};

MpfiInterval MpfiInterval::operator= (const MpfiInterval &i) {
	mpfi_set (mpfi (), i.mpfi ());
	return *this;
};

// 1
// 2
// comparisons with ints
bool MpfiInterval::operator== (const int n2) const {
	return (mpfi_is_inside_si ((long int)n2, mpfi ()) > 0);
};

bool MpfiInterval::operator!= (const int n2) const {
	return (mpfi_is_inside_si ((long int)n2, mpfi ()) == 0);
};

bool MpfiInterval::operator< (const int n2) const {
	return (mpfi_cmp_si (mpfi (), (long int)n2) < 0);
};

bool MpfiInterval::operator> (const int n2) const {
	return (mpfi_cmp_si (mpfi (), (long int)n2) > 0);
};

bool MpfiInterval::operator<= (const int n2) const {
	return (mpfi_cmp_si (mpfi (), (long int)n2) <= 0);
};

bool MpfiInterval::operator>= (const int n2) const {
	return (mpfi_cmp_si (mpfi (), (long int)n2) >= 0);
};
// comparisons with Gmpz
bool MpfiInterval::operator== (const CGAL::Gmpz &n2) const {
	return (mpfi_is_inside_z (n2.mpz(), mpfi ()) > 0);
};

bool MpfiInterval::operator!= (const CGAL::Gmpz &n2) const {
	return (mpfi_is_inside_z (n2.mpz(), mpfi ()) == 0);
};

bool MpfiInterval::operator< (const CGAL::Gmpz &n2) const {
	return (mpfi_cmp_z (mpfi (), n2.mpz()) < 0);
};

bool MpfiInterval::operator> (const CGAL::Gmpz &n2) const {
	return (mpfi_cmp_z (mpfi (), n2.mpz()) > 0);
};

bool MpfiInterval::operator<= (const CGAL::Gmpz &n2) const {
	return (mpfi_cmp_z (mpfi (), n2.mpz()) <= 0);
};

bool MpfiInterval::operator>= (const CGAL::Gmpz &n2) const {
	return (mpfi_cmp_z (mpfi (), n2.mpz()) >= 0);
};
// comparisons with Gmpq
bool MpfiInterval::operator== (const CGAL::Gmpq &n2) const {
	return (mpfi_is_inside_q (n2.mpq(), mpfi ()) > 0);
};

bool MpfiInterval::operator!= (const CGAL::Gmpq &n2) const {
	return (mpfi_is_inside_q (n2.mpq(), mpfi ()) == 0);
};

bool MpfiInterval::operator< (const CGAL::Gmpq &n2) const {
	return (mpfi_cmp_q (mpfi (), n2.mpq()) < 0);
};

bool MpfiInterval::operator> (const CGAL::Gmpq &n2) const {
	return (mpfi_cmp_q (mpfi (), n2.mpq()) > 0);
};

bool MpfiInterval::operator<= (const CGAL::Gmpq &n2) const {
	return (mpfi_cmp_q (mpfi (), n2.mpq()) <= 0);
};

bool MpfiInterval::operator>= (const CGAL::Gmpq &n2) const {
	return (mpfi_cmp_q (mpfi (), n2.mpq()) >= 0);
};

// 3
MpfiInterval MpfiInterval::operator+ (const MpfiInterval &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add (r, mpfi (), n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator- (const MpfiInterval &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_sub (r, mpfi (), n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator* (const MpfiInterval &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul (r, mpfi (), n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator- () const {
	mpfi_t n;
	mpfi_init (n);
	mpfi_neg (n, mpfi ());
	MpfiInterval ret (n);
	mpfi_clear (n);
	return ret;
}

MpfiInterval MpfiInterval::operator+= (const MpfiInterval &n2) {
	mpfi_add (mpfi (), mpfi (), n2.mpfi ());
	return *this;
};

MpfiInterval MpfiInterval::operator-= (const MpfiInterval &n2) {
	mpfi_sub (mpfi (), mpfi (), n2.mpfi ());
	return *this;
};

MpfiInterval MpfiInterval::operator*= (const MpfiInterval &n2) {
	mpfi_mul (mpfi (), mpfi (), n2.mpfi ());
	return *this;
};

// 4
// this (op) int
MpfiInterval MpfiInterval::operator+ (const int n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_si (r, mpfi (), (long int)n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator- (const int n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_sub_si (r, mpfi (), (long int)n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator* (const int n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_si (r, mpfi (), (long int)n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator+= (const int n2) {
	mpfi_add_si (mpfi (), mpfi (), (long int)n2);
	return *this;
};

MpfiInterval MpfiInterval::operator-= (const int n2) {
	mpfi_sub_si (mpfi (), mpfi (), (long int)n2);
	return *this;
};

MpfiInterval MpfiInterval::operator*= (const int n2) {
	mpfi_mul_si (mpfi (), mpfi (), (long int)n2);
	return *this;
};

// this (op) Gmpz
MpfiInterval MpfiInterval::operator+ (const CGAL::Gmpz &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_z (r, mpfi (), n2.mpz());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator- (const CGAL::Gmpz &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_sub_z (r, mpfi (), n2.mpz());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator* (const CGAL::Gmpz &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_z (r, mpfi (), n2.mpz());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator+= (const CGAL::Gmpz &n2) {
	mpfi_add_z (mpfi (), mpfi (), n2.mpz());
	return *this;
};

MpfiInterval MpfiInterval::operator-= (const CGAL::Gmpz &n2) {
	mpfi_sub_z (mpfi (), mpfi (), n2.mpz());
	return *this;
};

// this (op) Gmpq
MpfiInterval MpfiInterval::operator+ (const CGAL::Gmpq &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_q (r, mpfi (), n2.mpq());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator- (const CGAL::Gmpq &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_sub_q (r, mpfi (), n2.mpq());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator* (const CGAL::Gmpq &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_q (r, mpfi (), n2.mpq());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator+= (const CGAL::Gmpq &n2) {
	mpfi_add_q (mpfi (), mpfi (), n2.mpq());
	return *this;
};

MpfiInterval MpfiInterval::operator-= (const CGAL::Gmpq &n2) {
	mpfi_sub_q (mpfi (), mpfi (), n2.mpq());
	return *this;
};

MpfiInterval MpfiInterval::operator*= (const CGAL::Gmpq &n2) {
	mpfi_mul_q (mpfi (), mpfi (), n2.mpq());
	return *this;
};

// 5

// 6
bool MpfiInterval::is_valid () const {
	return (mpfi_nan_p (mpfi ()) == 0);
};

bool MpfiInterval::is_finite () const {
	return (mpfi_inf_p (mpfi ()) == 0);
};

double MpfiInterval::to_double () const {
	return mpfi_get_d (mpfi ());
};

std::pair <double, double> MpfiInterval::to_interval () const {
	mpfr_t temp;
	double left, right;
	mpfr_init (temp);
	mpfi_get_left (temp, mpfi ());
	left = mpfr_get_d (temp, GMP_RNDN);
	mpfi_get_left (temp, mpfi ());
	right = mpfr_get_d (temp, GMP_RNDN);
	mpfr_clear (temp);
	return std::make_pair (left, right);
};

// 7
MpfiInterval MpfiInterval::operator/ (const MpfiInterval &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_div (r, mpfi (), n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator/ (const int n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_div_si (r, mpfi (), (long int)n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator/= (const MpfiInterval &n2) {
	mpfi_div (mpfi (), mpfi (), n2.mpfi ());
	return *this;
};

MpfiInterval MpfiInterval::operator/= (const int n2) {
	mpfi_div_si (mpfi (), mpfi (), (long int)n2);
	return *this;
};

// 8
MpfiInterval MpfiInterval::sqrt () const {
	mpfi_t s;
	mpfi_init (s);
	mpfi_sqrt (s, mpfi ());
	MpfiInterval ret (s);
	mpfi_clear (s);
	return ret;
};

// 9
// TODO: rewrite this to better show the results
std::ostream& MpfiInterval::show (std::ostream &o) {
	char *str1, *str2;
	mpfr_t op1, op2;
	mp_exp_t *expptr1, *expptr2;

	expptr1 = (mp_exp_t*)malloc(sizeof(mp_exp_t));
	expptr2 = (mp_exp_t*)malloc(sizeof(mp_exp_t));

	mpfr_init (op1);
	mpfr_init (op2);

	mpfi_get_left (op1, mpfi ());
	mpfi_get_right (op2, mpfi ());

	str1 = mpfr_get_str (NULL, expptr1, 10, 0, op1, GMP_RNDN);
	str2 = mpfr_get_str (NULL, expptr2, 10, 0, op2, GMP_RNDN);

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
	mpfr_clear (op1);
	mpfr_clear (op2);

	return o;
};

// 10
MpfiInterval MpfiInterval::operator+ (const mpz_t &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_z (r, mpfi (), n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator- (const mpz_t &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_sub_z (r, mpfi (), n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator* (const mpz_t &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_z (r, mpfi (), n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator+= (const mpz_t &n2) {
	mpfi_add_z (mpfi (), mpfi (), n2);
	return *this;
};

MpfiInterval MpfiInterval::operator-= (const mpz_t &n2) {
	mpfi_sub_z (mpfi (), mpfi (), n2);
	return *this;
};

// this (op) mpq_t
MpfiInterval MpfiInterval::operator+ (const mpq_t &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_q (r, mpfi (), n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator- (const mpq_t &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_sub_q (r, mpfi (), n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator* (const mpq_t &n2) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_q (r, mpfi (), n2);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator+= (const mpq_t &n2) {
	mpfi_add_q (mpfi (), mpfi (), n2);
	return *this;
};

MpfiInterval MpfiInterval::operator-= (const mpq_t &n2) {
	mpfi_sub_q (mpfi (), mpfi (), n2);
	return *this;
};

MpfiInterval MpfiInterval::operator*= (const mpq_t &n2) {
	mpfi_mul_q (mpfi (), mpfi (), n2);
	return *this;
};

// 11. all the functions with mpfr_t that need to be inside the class

// constructor
MpfiInterval::MpfiInterval (const mpfr_t &r) {
	mpfi_set_fr (mpfi (), r);
};

// assigning: mpfi = mpfr
MpfiInterval MpfiInterval::operator= (const mpfr_t &r) {
	mpfi_set_fr (mpfi (), r);
	return *this;
};

// comparison: mpfi (op) mpfr
bool MpfiInterval::operator== (const mpfr_t &r) const {
	return (mpfi_cmp_fr (mpfi (), r) == 0);
};

bool MpfiInterval::operator!= (const mpfr_t &r) const {
	return (mpfi_cmp_fr ((long int)n2, r) != 0);
};

bool MpfiInterval::operator< (const mpfr_t &r) const {
	return (mpfi_cmp_fr (mpfi (), r) < 0);
};

bool MpfiInterval::operator> (const mpfr_t &r) const {
	return (mpfi_cmp_fr (mpfi (), r) > 0);
};

bool MpfiInterval::operator<= (const mpfr_t &r) const {
	return (mpfi_cmp_fr (mpfi (), r) <= 0);
};

bool MpfiInterval::operator>= (const mpfr_t &r) const {
	return (mpfi_cmp_fr (mpfi (), r) >= 0);
};

// arithmetics: mpfi (op) mpfr
MpfiInterval MpfiInterval::operator+ (const mpfr_t &f) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_fr (r, mpfi (), f);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator- (const mpfr_t &f) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_sub_fr (r, mpfi (), f);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator* (const mpfr_t &f) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_fr (r, mpfi (), f);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator/ (const mpfr_t &f) const {
	mpfi_t r;
	mpfi_init (r);
	mpfi_div_fr (r, mpfi (), f);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval MpfiInterval::operator+= (const mpfr_t &f) {
	mpfi_add_fr (mpfi (), mpfi (), f);
	return *this;
};

MpfiInterval MpfiInterval::operator-= (const mpfr_t &f) {
	mpfi_sub_fr (mpfi (), mpfi (), f);
	return *this;
};

MpfiInterval MpfiInterval::operator*= (const mpfr_t &f) {
	mpfi_mul_fr (mpfi (), mpfi (), f);
	return *this;
};

MpfiInterval MpfiInterval::operator/= (const mpfr_t &f) {
	mpfi_div_fr (mpfi (), mpfi (), f);
	return *this;
};




// ----------------------------------------
// end of the MpfiInterval class
// ----------------------------------------



// These functions are required, but they need to be coded outside the class:

// 1.5
bool operator== (const MpfiInterval &n1, const MpfiInterval &n2) {
	// mpfi_cmp returns 0 iff the two intervals have some point in common
	return (mpfi_cmp (n1.mpfi (), n2.mpfi ()) == 0);
}

bool operator!= (const MpfiInterval &n1, const MpfiInterval &n2) {
	return (mpfi_cmp (n1.mpfi (), n2.mpfi ()) != 0);
}

bool operator< (const MpfiInterval &n1, const MpfiInterval &n2) {
	return (mpfi_cmp (n1.mpfi (), n2.mpfi ()) < 0);
}

bool operator> (const MpfiInterval &n1, const MpfiInterval &n2) {
	return (mpfi_cmp (n1.mpfi (), n2.mpfi ()) > 0);
}

bool operator<= (const MpfiInterval &n1, const MpfiInterval &n2) {
	return (mpfi_cmp (n1.mpfi (), n2.mpfi ()) <= 0);
}

bool operator>= (const MpfiInterval &n1, const MpfiInterval &n2) {
	return (mpfi_cmp (n1.mpfi (), n2.mpfi ()) >= 0);
}

// 2.5
// -between ints and intervals
bool operator== (const int n1, const MpfiInterval &n2) {
	return (mpfi_is_inside_si ((long int)n1, n2.mpfi ()) > 0);
}

bool operator!= (const int n1, const MpfiInterval &n2) {
	return (mpfi_is_inside_si ((long int)n1, n2.mpfi ()) == 0);
}

bool operator< (const int n1, const MpfiInterval &n2) {
	return (mpfi_cmp_si (n2.mpfi (), (long int)n1) < 0);
}

bool operator> (const int n1, const MpfiInterval &n2) {
	return (mpfi_cmp_si (n2.mpfi (), (long int)n1) > 0);
}

bool operator<= (const int n1, const MpfiInterval &n2) {
	return (mpfi_cmp_si (n2.mpfi (), (long int)n1) <= 0);
}

bool operator>= (const int n1, const MpfiInterval &n2) {
	return (mpfi_cmp_si (n2.mpfi (), (long int)n1) >= 0);
}

// -between Gmpz and intervals
bool operator== (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	return (mpfi_is_inside_z (n1.mpz(), n2.mpfi ()) > 0);
}

bool operator!= (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	return (mpfi_is_inside_z (n1.mpz(), n2.mpfi ()) == 0);
}

bool operator< (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	return (mpfi_cmp_z (n2.mpfi (), n1.mpz()) < 0);
}

bool operator> (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	return (mpfi_cmp_z (n2.mpfi (), n1.mpz()) > 0);
}

bool operator<= (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	return (mpfi_cmp_z (n2.mpfi (), n1.mpz()) <= 0);
}

bool operator>= (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	return (mpfi_cmp_z (n2.mpfi (), n1.mpz()) >= 0);
}

// -between Gmpq and intervals
bool operator== (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	return (mpfi_is_inside_q (n1.mpq(), n2.mpfi ()) > 0);
}

bool operator!= (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	return (mpfi_is_inside_q (n1.mpq(), n2.mpfi ()) == 0);
}

bool operator< (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	return (mpfi_cmp_q (n2.mpfi (), n1.mpq()) < 0);
}

bool operator> (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	return (mpfi_cmp_q (n2.mpfi (), n1.mpq()) > 0);
}

bool operator<= (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	return (mpfi_cmp_q (n2.mpfi (), n1.mpq()) <= 0);
}

bool operator>= (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	return (mpfi_cmp_q (n2.mpfi (), n1.mpq()) >= 0);
}

// 4.5
// int (op) interval
MpfiInterval operator+ (const int n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_si (r, n2.mpfi (), (long int)n1);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator- (const int n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_si_sub (r, (long int)n1, n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator* (const int n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_si (r, n2.mpfi (), (long int)n1);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

// Gmpz int (op) interval
MpfiInterval operator+ (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_z (r, n2.mpfi (), n1.mpz());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator- (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_z_sub (r, n1.mpz(), n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator* (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_z (r, n2.mpfi (), n1.mpz());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

// Gmpq (op) interval
MpfiInterval operator+ (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_q (r, n2.mpfi (), n1.mpq());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator- (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_q_sub (r, n1.mpq(), n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator* (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_q (r, n2.mpfi (), n1.mpq());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

// 5.5
bool is_valid (const MpfiInterval &n) {
	return n.is_valid ();
};

bool is_finite (const MpfiInterval &n) {
	return n.is_finite ();
};

double to_double (const MpfiInterval &n) {
	return n.to_double ();
};

std::pair<double, double> to_interval (const MpfiInterval &n) {
	return n.to_interval ();
};

// 7.5
// anything / interval
MpfiInterval operator/ (const int n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_si_div (r, (long int)n1, n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval operator/ (const CGAL::Gmpz &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_z_div (r, n1.mpz(), n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

MpfiInterval operator/ (const CGAL::Gmpq &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_q_div (r, n1.mpq(), n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
};

// 8.5
MpfiInterval sqrt (const MpfiInterval &ntval) {
	return ntval.sqrt ();
};

// 9.5
std::ostream& operator<< (std::ostream &o, MpfiInterval &n) {
	return n.show(o);
};

// 10.5
// mpz_t (op) interval
MpfiInterval operator+ (const mpz_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_z (r, n2.mpfi (), n1);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator- (const mpz_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_z_sub (r, n1, n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator* (const mpz_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_z (r, n2.mpfi (), n1);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

// mpq_t (op) interval
MpfiInterval operator+ (const mpq_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_q (r, n2.mpfi (), n1);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator- (const mpq_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_q_sub (r, n1, n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator* (const mpq_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_q (r, n2.mpfi (), n1);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

// 11.5
// all the mpfr functions that can't be inside the class

// comparison
bool operator== (const mpfr_t &r, const MpfiInterval &n2) {
	return (mpfi_is_inside_fr ((long int)n1, n2.mpfi ()) > 0);
}

bool operator!= (const mpfr_t &r, const MpfiInterval &n2) {
	return (mpfi_is_inside_fr ((long int)n1, n2.mpfi ()) == 0);
}

bool operator< (const mpfr_t &r, const MpfiInterval &n2) {
	return (mpfi_cmp_fr (n2.mpfi (), (long int)n1) < 0);
}

bool operator> (const mpfr_t &r, const MpfiInterval &n2) {
	return (mpfi_cmp_fr (n2.mpfi (), (long int)n1) > 0);
}

bool operator<= (const mpfr_t &r, const MpfiInterval &n2) {
	return (mpfi_cmp_fr (n2.mpfi (), (long int)n1) <= 0);
}

bool operator>= (const mpfr_t &r, const MpfiInterval &n2) {
	return (mpfi_cmp_fr (n2.mpfi (), (long int)n1) >= 0);
}

// arithmetics
MpfiInterval operator+ (const mpfr_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_add_fr (r, n2.mpfi (), n1);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator- (const mpfr_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_fr_sub (r, n1, n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator* (const mpfr_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_mul_fr (r, n2.mpfi (), n1);
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

MpfiInterval operator/ (const mpfr_t &n1, const MpfiInterval &n2) {
	mpfi_t r;
	mpfi_init (r);
	mpfi_fr_div (r, n1, n2.mpfi ());
	MpfiInterval ret (r);
	mpfi_clear (r);
	return ret;
}

// not implemented: mpfr (op=) mpfi (because they must not return an interval)
// XXX: should them be implemented?

CGAL_END_NAMESPACE
