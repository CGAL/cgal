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
//	becames more readable (someone, someday, will do this)
//	-think about precision propagation in arithmetic functions
//	-enhance the exception mechanism (when comparison between intervals is
//	not known)
//	-add interfaces to more CGAL number types

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <iostream>
#include <mpfi.h>
#include <mpfi_io.h>
#include <mpfr.h>

CGAL_BEGIN_NAMESPACE

// the exception object
comparison_overlap_exn exn_overlap;

// what to do when a comparison fails?
void overlap () {
	throw exn_overlap;
}

// constructors of a "point" interval
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

// constructors of a "proper" interval
MpfiInterval::MpfiInterval (int l, int r) {
	mpfi_interv_si (mpfi (), (long int)l, (long int)r);
};

MpfiInterval::MpfiInterval (unsigned int l, unsigned int r) {
	mpfi_interv_ui (mpfi (), (unsigned long int)l, (unsigned long int)r);
};

MpfiInterval::MpfiInterval (long int l, long int r) {
	mpfi_interv_si (mpfi (), l, r);
};

MpfiInterval::MpfiInterval (unsigned long int l, unsigned long int r) {
	mpfi_interv_ui (mpfi (), l, r);
};

MpfiInterval::MpfiInterval (double l, double r) {
	mpfi_interv_d (mpfi (), l, r);
};

MpfiInterval::MpfiInterval (const CGAL::Gmpz &l, const CGAL::Gmpz &r) {
	mpfi_interv_z (mpfi (), l.mpz(), r.mpz());
};

MpfiInterval::MpfiInterval (const CGAL::Gmpq &l, const CGAL::Gmpq &r) {
	mpfi_interv_q (mpfi (), l.mpq(), r.mpq());
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

inline void MpfiInterval::get_left (mpfr_t &f) const {
	mpfi_get_left (f, mpfi ());
}

inline void MpfiInterval::get_right (mpfr_t &f) const {
	mpfi_get_right (f, mpfi ());
}

inline void MpfiInterval::get_endpoints (mpfr_t &l, mpfr_t &r) const {
	mpfi_get_left (l, mpfi ());
	mpfi_get_right (r, mpfi ());
}

inline bool MpfiInterval::is_point () const {
	mpfr_t l, r;
	mpfr_inits (l, r, NULL);
	get_endpoints (l, r);
	int comp = mpfr_equal_p (l, r);
	mpfr_clears (l, r, NULL);
	return (comp != 0);
}

inline bool MpfiInterval::contains (const int n) const {
	mpfr_t end;
	mpfr_init (end);
	int comp;
	get_left (end);	// first, we compare the left end
	comp = mpfr_cmp_si (end, n);
	if (comp > 0) {	// n is lower than the left end
		mpfr_clear (end);
		return false;
	}
	get_right (end);	// now, the right one
	comp = mpfr_cmp_si (end, n);
	if (comp < 0) {	// n is higher than the right end
		mpfr_clear (end);
		return false;
	}
	return true;
}

inline bool MpfiInterval::contains (const mpfr_t &n) const {
	mpfr_t end;
	mpfr_init (end);
	int comp;
	get_left (end);	// first, we compare the left end
	comp = mpfr_cmp (end, n);
	if (comp > 0) {	// n is lower than the left end
		mpfr_clear (end);
		return false;
	}
	get_right (end);	// now, the right one
	comp = mpfr_cmp (end, n);
	if (comp < 0) {	// n is higher than the right end
		mpfr_clear (end);
		return false;
	}
	return true;
}

inline bool MpfiInterval::contains (const mpz_t &n) const {
	mpfr_t end;
	mpfr_init (end);
	int comp;
	get_left (end);	// first, we compare the left end
	comp = mpfr_cmp_z (end, n);
	if (comp > 0) {	// n is lower than the left end
		mpfr_clear (end);
		return false;
	}
	get_right (end);	// now, the right one
	comp = mpfr_cmp_z (end, n);
	if (comp < 0) {	// n is higher than the right end
		mpfr_clear (end);
		return false;
	}
	return true;
}

inline bool MpfiInterval::contains (const mpq_t &n) const {
	mpfr_t end;
	mpfr_init (end);
	int comp;
	get_left (end);	// first, we compare the left end
	comp = mpfr_cmp_q (end, n);
	if (comp > 0) {	// n is lower than the left end
		mpfr_clear (end);
		return false;
	}
	get_right (end);	// now, the right one
	comp = mpfr_cmp_q (end, n);
	if (comp < 0) {	// n is higher than the right end
		mpfr_clear (end);
		return false;
	}
	return true;
}

inline bool MpfiInterval::contains (const Gmpz &n) const {
	mpfr_t end;
	mpfr_init (end);
	int comp;
	get_left (end);	// first, we compare the left end
	comp = mpfr_cmp_z (end, n.mpz());
	if (comp > 0) {	// n is lower than the left end
		mpfr_clear (end);
		return false;
	}
	get_right (end);	// now, the right one
	comp = mpfr_cmp_z (end, n.mpz());
	if (comp < 0) {	// n is higher than the right end
		mpfr_clear (end);
		return false;
	}
	return true;
}

inline bool MpfiInterval::contains (const Gmpq &n) const {
	mpfr_t end;
	mpfr_init (end);
	int comp;
	get_left (end);	// first, we compare the left end
	comp = mpfr_cmp_q (end, n.mpq());
	if (comp > 0) {	// n is lower than the left end
		mpfr_clear (end);
		return false;
	}
	get_right (end);	// now, the right one
	comp = mpfr_cmp_q (end, n.mpq());
	if (comp < 0) {	// n is higher than the right end
		mpfr_clear (end);
		return false;
	}
	return true;
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
	if (contains (n2)) {
		if (is_point ())
			return true;
		else
			overlap ();
	}
	return false;
};

bool MpfiInterval::operator!= (const int n2) const {
	return !(operator== (n2));
};

bool MpfiInterval::operator< (const int n2) const {
	if (contains (n2))
		if (is_point ())
			return false;
		else
			overlap ();
	mpfr_t end;
	mpfr_init (end);
	get_right (end);	// end is the right endpoint
	int comp1 = mpfr_cmp_si (end, n2);
	if (comp1 < 0) {
		mpfr_clear (end);
		return true;
	}
	mpfr_clear (end);
	return false;
};

bool MpfiInterval::operator> (const int n2) const {
	if (contains (n2))
		if (is_point ())
			return false;
		else
			overlap ();
	mpfr_t end;
	mpfr_init (end);
	get_left (end);	// end is the left endpoint
	int comp1 = mpfr_cmp_si (end, n2);
	if (comp1 > 0) {
		mpfr_clear (end);
		return true;
	}
	mpfr_clear (end);
	return false;
};

bool MpfiInterval::operator<= (const int n2) const {
	return ((operator== (n2)) || (operator< (n2)));
};

bool MpfiInterval::operator>= (const int n2) const {
	return ((operator== (n2)) || (operator> (n2)));
};

// comparisons with Gmpz and Gmpq
template <class T>
bool MpfiInterval::operator== (const T &n2) const {
	if (contains (n2)) {
		if (is_point ())
			return true;
		else
			overlap ();
	}
	return false;
};

template <class T>
bool MpfiInterval::operator!= (const T &n2) const {
	return !(operator== (n2));
};

bool MpfiInterval::operator< (const CGAL::Gmpz &n2) const {
	if (contains (n2))
		if (is_point ())
			return false;
		else
			overlap ();
	mpfr_t end;
	mpfr_init (end);
	get_right (end);	// end is the right endpoint
	int comp1 = mpfr_cmp_z (end, n2.mpz());
	if (comp1 < 0) {
		mpfr_clear (end);
		return true;
	}
	mpfr_clear (end);
	return false;
};

bool MpfiInterval::operator< (const CGAL::Gmpq &n2) const {
	if (contains (n2))
		if (is_point ())
			return false;
		else
			overlap ();
	mpfr_t end;
	mpfr_init (end);
	get_right (end);	// end is the right endpoint
	int comp1 = mpfr_cmp_q (end, n2.mpq());
	if (comp1 < 0) {
		mpfr_clear (end);
		return true;
	}
	mpfr_clear (end);
	return false;
};

bool MpfiInterval::operator> (const CGAL::Gmpz &n2) const {
	if (contains (n2))
		overlap ();
	mpfr_t end;
	mpfr_init (end);
	get_left (end);	// end is the left endpoint
	int comp1 = mpfr_cmp_z (end, n2.mpz());
	if (comp1 > 0) {
		mpfr_clear (end);
		return true;
	}
	mpfr_clear (end);
	return false;
};

bool MpfiInterval::operator> (const CGAL::Gmpq &n2) const {
	if (contains (n2))
		overlap ();
	mpfr_t end;
	mpfr_init (end);
	get_left (end);	// end is the left endpoint
	int comp1 = mpfr_cmp_q (end, n2.mpq());
	if (comp1 > 0) {
		mpfr_clear (end);
		return true;
	}
	mpfr_clear (end);
	return false;
};

template <class T>
bool MpfiInterval::operator<= (const T &n2) const {
	return ((operator== (n2)) || (operator< (n2)));
};

template <class T>
bool MpfiInterval::operator>= (const T &n2) const {
	return ((operator== (n2)) || (operator> (n2)));
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

	mpfr_inits (op1, op2, NULL);

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
	mpfr_clears (op1, op2, NULL);

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

// constructor I
MpfiInterval::MpfiInterval (const mpfr_t &r) {
	mpfi_set_fr (mpfi (), r);
};

// constructor II
MpfiInterval::MpfiInterval (const mpfr_t &l, const mpfr_t &r) {
	mpfi_interv_fr (mpfi (), l, r);
};

// assigning: mpfi = mpfr
MpfiInterval MpfiInterval::operator= (const mpfr_t &r) {
	mpfi_set_fr (mpfi (), r);
	return *this;
};

// comparison: mpfi (op) mpfr
//	NOTE: the previous template definitions of operators =, !=, >= and <=
//	should work with mpfr_t
bool MpfiInterval::operator< (const mpfr_t &n2) const {
	if (contains (n2))
		if (is_point ())
			return false;
		else
			overlap ();
	mpfr_t end;
	mpfr_init (end);
	get_right (end);	// end is the right endpoint
	int comp1 = mpfr_cmp (end, n2);
	if (comp1 < 0) {
		mpfr_clear (end);
		return true;
	}
	mpfr_clear (end);
	return false;
};

bool MpfiInterval::operator> (const mpfr_t &n2) const {
	if (contains (n2))
		overlap ();
	mpfr_t end;
	mpfr_init (end);
	get_left (end);	// end is the left endpoint
	int comp1 = mpfr_cmp (end, n2);
	if (comp1 > 0) {
		mpfr_clear (end);
		return true;
	}
	mpfr_clear (end);
	return false;
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

// 1.5 TODO
bool operator== (const MpfiInterval &n1, const MpfiInterval &n2) {
	mpfr_t n1_l, n1_r, n2_l, n2_r;
	mpfr_inits (n1_l, n1_r, n2_l, n2_r, NULL);
	n1.get_endpoints (n1_l, n1_r);
	n2.get_endpoints (n2_l, n2_r);
	// we assume here that the left endpoint is lesser or equal than the
	// right one (it is responsibility of the user to assure that)
	if ((mpfr_cmp (n1_l, n2_l) == 0) && (mpfr_cmp (n1_r, n2_r) == 0)) {
		mpfr_clears (n1_l, n1_r, n2_l, n2_r, NULL);
		return true;
	}
	if ((mpfr_cmp (n1_r, n2_l) < 0) || (mpfr_cmp (n2_r, n1_l) < 0)) {
		mpfr_clears (n1_l, n1_r, n2_l, n2_r, NULL);
		return false;
	}
	mpfr_clears (n1_l, n1_r, n2_l, n2_r, NULL);
	overlap ();
	return false;	// this never occurs
}

bool operator!= (const MpfiInterval &n1, const MpfiInterval &n2) {
	return !(n1 == n2);
}

bool operator< (const MpfiInterval &n1, const MpfiInterval &n2) {
	mpfr_t n1_r, n2_l;
	mpfr_inits (n1_r, n2_l, NULL);
	n1.get_right (n1_r);
	n2.get_left (n2_l);
	if (mpfr_cmp (n1_r, n2_l) < 0) {
		mpfr_clears (n1_r, n2_l, NULL);
		return true;
	}
	mpfr_clears (n1_r, n2_l, NULL);
	mpfr_t n1_l, n2_r;
	mpfr_inits (n1_l, n2_r, NULL);
	n1.get_left (n1_l);
	n2.get_right (n2_r);
	if (mpfr_cmp (n2_r, n1_l) < 0) {
		mpfr_clears (n1_l, n2_r, NULL);
		return false;
	}
	mpfr_clears (n1_l, n2_r, NULL);
	overlap ();
	return false;	// this never occurs
}

bool operator> (const MpfiInterval &n1, const MpfiInterval &n2) {
	return (n2 < n1);
}

bool operator<= (const MpfiInterval &n1, const MpfiInterval &n2) {
	return ((n1 == n2) || (n1 < n2));
}

bool operator>= (const MpfiInterval &n1, const MpfiInterval &n2) {
	return ((n1 == n2) || (n1 > n2));
}

// 2.5
// comparison between int|mpfr_t|mp[zq]_t|Gmp[zq] and intervals
template <class T>
bool operator== (const T &n1, const MpfiInterval &n2) {
	return (n2 == n1);
}

template <class T>
bool operator!= (const T &n1, const MpfiInterval &n2) {
	return (n2 != n1);
}

template <class T>
bool operator< (const T &n1, const MpfiInterval &n2) {
	return (n2 > n1);
}

template <class T>
bool operator> (const T &n1, const MpfiInterval &n2) {
	return (n2 < n1);
}

template <class T>
bool operator<= (const T &n1, const MpfiInterval &n2) {
	return (n2 >= n1);
}

template <class T>
bool operator>= (const T &n1, const MpfiInterval &n2) {
	return (n2 <= n1);
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

// 11.5
// all the mpfr functions that can't be inside the class (and aren't covered
// by the template functions

// arithmetics
template <class T>
MpfiInterval operator+ (const T &n1, const MpfiInterval &n2) {
	return (n2 + n1);
}

template <class T>
MpfiInterval operator- (const T &n1, const MpfiInterval &n2) {
	return ((-n2) + n1);
}

template <class T>
MpfiInterval operator* (const T &n1, const MpfiInterval &n2) {
	return (n2 * n1);
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
