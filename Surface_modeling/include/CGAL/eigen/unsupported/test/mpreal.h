/*
	Multi-precision real number class. C++ wrapper fo MPFR library.
	Project homepage: http://www.holoborodko.com/pavel/
	Contact e-mail:   pavel@holoborodko.com

	Copyright (c) 2008-2010 Pavel Holoborodko

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

	Contributors:
	Brian Gladman, Helmut Jarausch, Fokko Beekhof, Ulrich Mutze, 
	Heinz van Saanen, Pere Constans, Dmitriy Gubanov
*/

#ifndef __MP_REAL_H__
#define __MP_REAL_H__

#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cfloat>
#include <cmath>

#include <mpfr.h>

// Detect compiler using signatures from http://predef.sourceforge.net/
#if defined(__GNUC__) && defined(__INTEL_COMPILER)
	#define IsInf(x) isinf(x)                                // GNU C/C++ + Intel ICC compiler

#elif defined(__GNUC__)
	#define IsInf(x) std::isinf(x)                          // GNU C/C++ 

#elif defined(_MSC_VER)		
	#define IsInf(x) (!_finite(x))				// Microsoft Visual C++

#else
	#define IsInf(x) std::isinf(x)				// C99 conformance
#endif

namespace mpfr {

class mpreal {
private:
	mpfr_t mp;

public:
	static mp_rnd_t		default_rnd;	
	static mp_prec_t	default_prec;	
	static int			default_base;
	static int			double_bits;

public:
	// Constructors && type conversion
	mpreal();
	mpreal(const mpreal& u);

	mpreal(const mpfr_t u);	
	mpreal(const mpf_t u);	

	mpreal(const mpz_t u, mp_prec_t prec = default_prec, mp_rnd_t mode = default_rnd);	
	mpreal(const mpq_t u, mp_prec_t prec = default_prec, mp_rnd_t mode = default_rnd);	
	mpreal(const double u, mp_prec_t prec = default_prec, mp_rnd_t mode = default_rnd);
	mpreal(const long double u, mp_prec_t prec = default_prec, mp_rnd_t mode = default_rnd);
	mpreal(const unsigned long int u, mp_prec_t prec = default_prec, mp_rnd_t mode = default_rnd);
	mpreal(const unsigned int u, mp_prec_t prec = default_prec, mp_rnd_t mode = default_rnd);
	mpreal(const long int u, mp_prec_t prec = default_prec, mp_rnd_t mode = default_rnd);
	mpreal(const int u, mp_prec_t prec = default_prec, mp_rnd_t mode = default_rnd);
	mpreal(const char* s, mp_prec_t prec = default_prec, int base = default_base, mp_rnd_t mode = default_rnd);

	~mpreal();                           

	// Operations
	// =
	// +, -, *, /, ++, --, <<, >> 
	// *=, +=, -=, /=,
	// <, >, ==, <=, >=

	// =
	mpreal& operator=(const mpreal& v);
	mpreal& operator=(const mpf_t v);
	mpreal& operator=(const mpz_t v);
	mpreal& operator=(const mpq_t v);
	mpreal& operator=(const long double v);
	mpreal& operator=(const double v);		
	mpreal& operator=(const unsigned long int v);
	mpreal& operator=(const unsigned int v);
	mpreal& operator=(const long int v);
	mpreal& operator=(const int v);
	mpreal& operator=(const char* s);

	// +
	mpreal& operator+=(const mpreal& v);
	mpreal& operator+=(const mpf_t v);
	mpreal& operator+=(const mpz_t v);
	mpreal& operator+=(const mpq_t v);
	mpreal& operator+=(const long double u);
	mpreal& operator+=(const double u);
	mpreal& operator+=(const unsigned long int u);
	mpreal& operator+=(const unsigned int u);
	mpreal& operator+=(const long int u);
	mpreal& operator+=(const int u);
	const mpreal operator+() const;
	mpreal& operator++ ();
	const mpreal  operator++ (int); 

	// -
	mpreal& operator-=(const mpreal& v);
	mpreal& operator-=(const mpz_t v);
	mpreal& operator-=(const mpq_t v);
	mpreal& operator-=(const long double u);
	mpreal& operator-=(const double u);
	mpreal& operator-=(const unsigned long int u);
	mpreal& operator-=(const unsigned int u);
	mpreal& operator-=(const long int u);
	mpreal& operator-=(const int u);
	const mpreal operator-() const;
	friend const mpreal operator-(const unsigned long int b, const mpreal& a);
	friend const mpreal operator-(const unsigned int b, const mpreal& a);
	friend const mpreal operator-(const long int b, const mpreal& a);
	friend const mpreal operator-(const int b, const mpreal& a);
	friend const mpreal operator-(const double b, const mpreal& a);
	mpreal& operator-- ();    
	const mpreal  operator-- (int);

	// *
	mpreal& operator*=(const mpreal& v);
	mpreal& operator*=(const mpz_t v);
	mpreal& operator*=(const mpq_t v);
	mpreal& operator*=(const long double v);
	mpreal& operator*=(const double v);
	mpreal& operator*=(const unsigned long int v);
	mpreal& operator*=(const unsigned int v);
	mpreal& operator*=(const long int v);
	mpreal& operator*=(const int v);
	
	// /
	mpreal& operator/=(const mpreal& v);
	mpreal& operator/=(const mpz_t v);
	mpreal& operator/=(const mpq_t v);
	mpreal& operator/=(const long double v);
	mpreal& operator/=(const double v);
	mpreal& operator/=(const unsigned long int v);
	mpreal& operator/=(const unsigned int v);
	mpreal& operator/=(const long int v);
	mpreal& operator/=(const int v);
	friend const mpreal operator/(const unsigned long int b, const mpreal& a);
	friend const mpreal operator/(const unsigned int b, const mpreal& a);
	friend const mpreal operator/(const long int b, const mpreal& a);
	friend const mpreal operator/(const int b, const mpreal& a);
	friend const mpreal operator/(const double b, const mpreal& a);

	//<<= Fast Multiplication by 2^u
	mpreal& operator<<=(const unsigned long int u);
	mpreal& operator<<=(const unsigned int u);
	mpreal& operator<<=(const long int u);
	mpreal& operator<<=(const int u);

	//>>= Fast Division by 2^u
	mpreal& operator>>=(const unsigned long int u);
	mpreal& operator>>=(const unsigned int u);
	mpreal& operator>>=(const long int u);
	mpreal& operator>>=(const int u);

	// Boolean Operators
	friend bool operator >  (const mpreal& a, const mpreal& b);
	friend bool operator >= (const mpreal& a, const mpreal& b);
	friend bool operator <  (const mpreal& a, const mpreal& b);
	friend bool operator <= (const mpreal& a, const mpreal& b);
	friend bool operator == (const mpreal& a, const mpreal& b);
	friend bool operator != (const mpreal& a, const mpreal& b);

	// Type Conversion operators
	inline operator long double() const;
	inline operator double() const;
	inline operator float() const;
	inline operator unsigned long() const;
	inline operator unsigned int() const;
	inline operator long() const;
	operator std::string() const;
	inline operator mpfr_ptr();

	// Math Functions
	friend const mpreal sqr(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal sqrt(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal sqrt(const unsigned long int v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal cbrt(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal root(const mpreal& v, unsigned long int k, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal pow(const mpreal& a, const mpreal& b, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal pow(const mpreal& a, const mpz_t b, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal pow(const mpreal& a, const unsigned long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal pow(const mpreal& a, const long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal pow(const unsigned long int a, const mpreal& b, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal pow(const unsigned long int a, const unsigned long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal fabs(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal abs(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal dim(const mpreal& a, const mpreal& b, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend inline const mpreal mul_2ui(const mpreal& v, unsigned long int k, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend inline const mpreal mul_2si(const mpreal& v, long int k, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend inline const mpreal div_2ui(const mpreal& v, unsigned long int k, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend inline const mpreal div_2si(const mpreal& v, long int k, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend int cmpabs(const mpreal& a,const mpreal& b);
	
	friend const mpreal log  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal log2 (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal log10(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal exp  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd); 
	friend const mpreal exp2 (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal exp10(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);

	friend const mpreal cos(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal sin(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal tan(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal sec(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal csc(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal cot(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend int sin_cos(mpreal& s, mpreal& c, const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);

	friend const mpreal acos  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal asin  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal atan  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal atan2 (const mpreal& y, const mpreal& x, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal cosh  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal sinh  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal tanh  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal sech  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal csch  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal coth  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);

	friend const mpreal acosh  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal asinh  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal atanh  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal fac_ui (unsigned long int v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal log1p  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal expm1  (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal eint   (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);

	friend const mpreal gamma (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal lngamma (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal lgamma (const mpreal& v, int *signp, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal zeta (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal erf (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal erfc (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal _j0 (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd); 
	friend const mpreal _j1 (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd); 
	friend const mpreal _jn (long n, const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal _y0 (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal _y1 (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal _yn (long n, const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd); 
	friend const mpreal fma (const mpreal& v1, const mpreal& v2, const mpreal& v3, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal fms (const mpreal& v1, const mpreal& v2, const mpreal& v3, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal agm (const mpreal& v1, const mpreal& v2, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal hypot (const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal sum (const mpreal tab[], unsigned long int n, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend int sgn(const mpreal& v); // -1 or +1

// MPFR 2.4.0 Specifics
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
	friend int sinh_cosh(mpreal& s, mpreal& c, const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal li2(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal fmod (const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal rec_sqrt(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
#endif

// MPFR 3.0.0 Specifics
#if (MPFR_VERSION >= MPFR_VERSION_NUM(3,0,0))
	friend const mpreal digamma(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
    friend const mpreal urandom (gmp_randstate_t& state,mp_rnd_t rnd_mode = mpreal::default_rnd); 	// use gmp_randinit_default() to init state, gmp_randclear() to clear
	friend bool _isregular(const mpreal& v);
#endif

	// Exponent and mantissa manipulation
	friend const mpreal frexp(const mpreal& v, mp_exp_t* exp);	
	friend const mpreal ldexp(const mpreal& v, mp_exp_t exp);

	// Splits mpreal value into fractional and integer parts.
	// Returns fractional part and stores integer part in n.
	friend const mpreal modf(const mpreal& v, mpreal& n);	

	// Constants
	// don't forget to call mpfr_free_cache() for every thread where you are using const-functions
	friend const mpreal const_log2 (mp_prec_t prec = mpreal::default_prec, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal const_pi (mp_prec_t prec = mpreal::default_prec, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal const_euler (mp_prec_t prec = mpreal::default_prec, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal const_catalan (mp_prec_t prec = mpreal::default_prec, mp_rnd_t rnd_mode = mpreal::default_rnd);
	// returns +inf iff sign>=0 otherwise -inf
	friend const mpreal const_infinity(int sign = 1, mp_prec_t prec = mpreal::default_prec, mp_rnd_t rnd_mode = mpreal::default_rnd);

	// Output/ Input
	friend std::ostream& operator<<(std::ostream& os, const mpreal& v);
    friend std::istream& operator>>(std::istream& is, mpreal& v);

	// Integer Related Functions
	friend const mpreal rint (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal ceil (const mpreal& v);
	friend const mpreal floor(const mpreal& v);
	friend const mpreal round(const mpreal& v);
	friend const mpreal trunc(const mpreal& v);
	friend const mpreal rint_ceil (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal rint_floor(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal rint_round(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal rint_trunc(const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal frac (const mpreal& v, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal remainder (const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode = mpreal::default_rnd);
	friend const mpreal remquo (long* q, const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode = mpreal::default_rnd);
	
	// Miscellaneous Functions
	friend const mpreal nexttoward (const mpreal& x, const mpreal& y);
	friend const mpreal nextabove  (const mpreal& x);
	friend const mpreal nextbelow  (const mpreal& x);

	// use gmp_randinit_default() to init state, gmp_randclear() to clear
	friend const mpreal urandomb (gmp_randstate_t& state); 

// MPFR < 2.4.2 Specifics
#if (MPFR_VERSION <= MPFR_VERSION_NUM(2,4,2))
	friend const mpreal random2 (mp_size_t size, mp_exp_t exp);
#endif

	// Instance Checkers
	friend bool _isnan(const mpreal& v);
	friend bool _isinf(const mpreal& v);
	friend bool _isnum(const mpreal& v);
	friend bool _iszero(const mpreal& v);
	friend bool _isint(const mpreal& v);

	// Set/Get instance properties
	inline mp_prec_t	get_prec() const;
	inline void		set_prec(mp_prec_t prec, mp_rnd_t rnd_mode = default_rnd);	// Change precision with rounding mode
	
	// Set mpreal to +-inf, NaN
	void      set_inf(int sign = +1);	
	void	  set_nan();

	// sign = -1 or +1
	void set_sign(int sign, mp_rnd_t rnd_mode = default_rnd);

	//Exponent
	mp_exp_t get_exp();
	int set_exp(mp_exp_t e);
	int check_range (int t, mp_rnd_t rnd_mode = default_rnd);
	int subnormalize (int t,mp_rnd_t rnd_mode = default_rnd);

	// Inexact conversion from float
	inline bool fits_in_bits(double x, int n);

	// Set/Get global properties
	static void			set_default_prec(mp_prec_t prec);
	static mp_prec_t	get_default_prec();
	static void			set_default_base(int base);
	static int			get_default_base();
	static void			set_double_bits(int dbits);
	static int			get_double_bits();
	static void			set_default_rnd(mp_rnd_t rnd_mode);
	static mp_rnd_t		get_default_rnd();
	static mp_exp_t get_emin (void);
	static mp_exp_t get_emax (void);
	static mp_exp_t get_emin_min (void);
	static mp_exp_t get_emin_max (void);
	static mp_exp_t get_emax_min (void);
	static mp_exp_t get_emax_max (void);
	static int set_emin (mp_exp_t exp);
	static int set_emax (mp_exp_t exp);

	// Get/Set conversions
	// Convert mpreal to string with n significant digits in base b
	// n = 0 -> convert with the maximum available digits 
	std::string to_string(size_t n = 0, int b = default_base, mp_rnd_t mode = default_rnd) const;
	
	// Efficient swapping of two mpreal values
	friend void swap(mpreal& x, mpreal& y);
	
	//Min Max - macros is evil. Needed for systems which defines max and min globally as macros (e.g. Windows)
	//Hope that globally defined macros use > < operations only
	#ifndef max
		friend const mpreal max(const mpreal& x, const mpreal& y);
	#endif

	#ifndef min
		friend const mpreal min(const mpreal& x, const mpreal& y);
	#endif
};

//////////////////////////////////////////////////////////////////////////
// Exceptions
class conversion_overflow : public std::exception {
public:
	std::string why() { return "inexact conversion from floating point"; }
};

//////////////////////////////////////////////////////////////////////////
// + Addition
const mpreal operator+(const mpreal& a, const mpreal& b);

// + Fast specialized addition - implemented through fast += operations
const mpreal operator+(const mpreal& a, const mpz_t b);
const mpreal operator+(const mpreal& a, const mpq_t b);
const mpreal operator+(const mpreal& a, const long double b);
const mpreal operator+(const mpreal& a, const double b);
const mpreal operator+(const mpreal& a, const unsigned long int b);
const mpreal operator+(const mpreal& a, const unsigned int b);
const mpreal operator+(const mpreal& a, const long int b);
const mpreal operator+(const mpreal& a, const int b);
const mpreal operator+(const mpreal& a, const char* b);
const mpreal operator+(const char* a, const mpreal& b);
const std::string operator+(const mpreal& a, const std::string b);
const std::string operator+(const std::string a, const mpreal& b);

const mpreal operator+(const mpz_t b, const mpreal& a);
const mpreal operator+(const mpq_t b, const mpreal& a);
const mpreal operator+(const long double b, const mpreal& a);
const mpreal operator+(const double  b, const mpreal& a);
const mpreal operator+(const unsigned long int b, const mpreal& a);
const mpreal operator+(const unsigned int b, const mpreal& a);
const mpreal operator+(const long int b, const mpreal& a);
const mpreal operator+(const int b, const mpreal& a);

//////////////////////////////////////////////////////////////////////////
// - Subtraction
const mpreal operator-(const mpreal& a, const mpreal& b);

// - Fast specialized subtraction - implemented through fast -= operations
const mpreal operator-(const mpreal& a, const mpz_t b);
const mpreal operator-(const mpreal& a, const mpq_t b);
const mpreal operator-(const mpreal& a, const long double b);
const mpreal operator-(const mpreal& a, const double b);
const mpreal operator-(const mpreal& a, const unsigned long int b);
const mpreal operator-(const mpreal& a, const unsigned int b);
const mpreal operator-(const mpreal& a, const long int b);
const mpreal operator-(const mpreal& a, const int b);
const mpreal operator-(const mpreal& a, const char* b);
const mpreal operator-(const char* a, const mpreal& b);

const mpreal operator-(const mpz_t b, const mpreal& a);
const mpreal operator-(const mpq_t b, const mpreal& a);
const mpreal operator-(const long double b, const mpreal& a);
//const mpreal operator-(const double  b, const mpreal& a);

//////////////////////////////////////////////////////////////////////////
// * Multiplication
const mpreal operator*(const mpreal& a, const mpreal& b);

// * Fast specialized multiplication - implemented through fast *= operations
const mpreal operator*(const mpreal& a, const mpz_t b);
const mpreal operator*(const mpreal& a, const mpq_t b);
const mpreal operator*(const mpreal& a, const long double b);
const mpreal operator*(const mpreal& a, const double b);
const mpreal operator*(const mpreal& a, const unsigned long int b);
const mpreal operator*(const mpreal& a, const unsigned int b);
const mpreal operator*(const mpreal& a, const long int b);
const mpreal operator*(const mpreal& a, const int b);

const mpreal operator*(const mpz_t b, const mpreal& a);
const mpreal operator*(const mpq_t b, const mpreal& a);
const mpreal operator*(const long double b, const mpreal& a);
const mpreal operator*(const double  b, const mpreal& a);
const mpreal operator*(const unsigned long int b, const mpreal& a);
const mpreal operator*(const unsigned int b, const mpreal& a);
const mpreal operator*(const long int b, const mpreal& a);
const mpreal operator*(const int b, const mpreal& a);

//////////////////////////////////////////////////////////////////////////
// / Division
const mpreal operator/(const mpreal& a, const mpreal& b);

// / Fast specialized division - implemented through fast /= operations
const mpreal operator/(const mpreal& a, const mpz_t b);
const mpreal operator/(const mpreal& a, const mpq_t b);
const mpreal operator/(const mpreal& a, const long double b);
const mpreal operator/(const mpreal& a, const double b);
const mpreal operator/(const mpreal& a, const unsigned long int b);
const mpreal operator/(const mpreal& a, const unsigned int b);
const mpreal operator/(const mpreal& a, const long int b);
const mpreal operator/(const mpreal& a, const int b);

const mpreal operator/(const long double b, const mpreal& a);

//////////////////////////////////////////////////////////////////////////
// Shifts operators - Multiplication/Division by a power of 2
const mpreal operator<<(const mpreal& v, const unsigned long int k);
const mpreal operator<<(const mpreal& v, const unsigned int k);
const mpreal operator<<(const mpreal& v, const long int k);
const mpreal operator<<(const mpreal& v, const int k);

const mpreal operator>>(const mpreal& v, const unsigned long int k);
const mpreal operator>>(const mpreal& v, const unsigned int k);
const mpreal operator>>(const mpreal& v, const long int k);
const mpreal operator>>(const mpreal& v, const int k);

//////////////////////////////////////////////////////////////////////////
// Boolean operators
bool operator <  (const mpreal& a, const unsigned long int b);
bool operator <  (const mpreal& a, const unsigned int b);
bool operator <  (const mpreal& a, const long int b);
bool operator <  (const mpreal& a, const int b);
bool operator <  (const mpreal& a, const long double b);
bool operator <  (const mpreal& a, const double b);

bool operator <  (const unsigned long int a,const mpreal& b);
bool operator <  (const unsigned int a,		const mpreal& b);
bool operator <  (const long int a,			const mpreal& b);
bool operator <  (const int a,				const mpreal& b);
bool operator <  (const long double a,		const mpreal& b);
bool operator <  (const double a,			const mpreal& b);

bool operator >  (const mpreal& a, const unsigned long int b);
bool operator >  (const mpreal& a, const unsigned int b);
bool operator >  (const mpreal& a, const long int b);
bool operator >  (const mpreal& a, const int b);
bool operator >  (const mpreal& a, const long double b);
bool operator >  (const mpreal& a, const double b);

bool operator >  (const unsigned long int a,const mpreal& b);
bool operator >  (const unsigned int a,		const mpreal& b);
bool operator >  (const long int a,			const mpreal& b);
bool operator >  (const int a,				const mpreal& b);
bool operator >  (const long double a,		const mpreal& b);
bool operator >  (const double a,			const mpreal& b);

bool operator >=  (const mpreal& a, const unsigned long int b);
bool operator >=  (const mpreal& a, const unsigned int b);
bool operator >=  (const mpreal& a, const long int b);
bool operator >=  (const mpreal& a, const int b);
bool operator >=  (const mpreal& a, const long double b);
bool operator >=  (const mpreal& a, const double b);

bool operator >=  (const unsigned long int a,const mpreal& b);
bool operator >=  (const unsigned int a,		const mpreal& b);
bool operator >=  (const long int a,			const mpreal& b);
bool operator >=  (const int a,				const mpreal& b);
bool operator >=  (const long double a,		const mpreal& b);
bool operator >=  (const double a,			const mpreal& b);

bool operator <=  (const mpreal& a, const unsigned long int b);
bool operator <=  (const mpreal& a, const unsigned int b);
bool operator <=  (const mpreal& a, const long int b);
bool operator <=  (const mpreal& a, const int b);
bool operator <=  (const mpreal& a, const long double b);
bool operator <=  (const mpreal& a, const double b);

bool operator <=  (const unsigned long int a,const mpreal& b);
bool operator <=  (const unsigned int a,		const mpreal& b);
bool operator <=  (const long int a,			const mpreal& b);
bool operator <=  (const int a,				const mpreal& b);
bool operator <=  (const long double a,		const mpreal& b);
bool operator <=  (const double a,			const mpreal& b);

bool operator ==  (const mpreal& a, const unsigned long int b);
bool operator ==  (const mpreal& a, const unsigned int b);
bool operator ==  (const mpreal& a, const long int b);
bool operator ==  (const mpreal& a, const int b);
bool operator ==  (const mpreal& a, const long double b);
bool operator ==  (const mpreal& a, const double b);

bool operator ==  (const unsigned long int a,const mpreal& b);
bool operator ==  (const unsigned int a,		const mpreal& b);
bool operator ==  (const long int a,			const mpreal& b);
bool operator ==  (const int a,				const mpreal& b);
bool operator ==  (const long double a,		const mpreal& b);
bool operator ==  (const double a,			const mpreal& b);

bool operator !=  (const mpreal& a, const unsigned long int b);
bool operator !=  (const mpreal& a, const unsigned int b);
bool operator !=  (const mpreal& a, const long int b);
bool operator !=  (const mpreal& a, const int b);
bool operator !=  (const mpreal& a, const long double b);
bool operator !=  (const mpreal& a, const double b);

bool operator !=  (const unsigned long int a,const mpreal& b);
bool operator !=  (const unsigned int a,		const mpreal& b);
bool operator !=  (const long int a,			const mpreal& b);
bool operator !=  (const int a,				const mpreal& b);
bool operator !=  (const long double a,		const mpreal& b);
bool operator !=  (const double a,			const mpreal& b);

//////////////////////////////////////////////////////////////////////////
// sqrt
const mpreal sqrt(const unsigned int v, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal sqrt(const long int v, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal sqrt(const int v, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal sqrt(const long double v, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal sqrt(const double v, mp_rnd_t rnd_mode = mpreal::default_rnd);

//////////////////////////////////////////////////////////////////////////
// pow
const mpreal pow(const mpreal& a, const unsigned int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const mpreal& a, const int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const mpreal& a, const long double b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const mpreal& a, const double b, mp_rnd_t rnd_mode = mpreal::default_rnd);

const mpreal pow(const unsigned int a, const mpreal& b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long int a, const mpreal& b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const int a, const mpreal& b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long double a, const mpreal& b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const double a, const mpreal& b, mp_rnd_t rnd_mode = mpreal::default_rnd);

const mpreal pow(const unsigned long int a, const unsigned int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned long int a, const long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned long int a, const int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned long int a, const long double b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned long int a, const double b, mp_rnd_t rnd_mode = mpreal::default_rnd);

const mpreal pow(const unsigned int a, const unsigned long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned int a, const unsigned int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned int a, const long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned int a, const int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned int a, const long double b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const unsigned int a, const double b, mp_rnd_t rnd_mode = mpreal::default_rnd);

const mpreal pow(const long int a, const unsigned long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long int a, const unsigned int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long int a, const long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long int a, const int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long int a, const long double b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long int a, const double b, mp_rnd_t rnd_mode = mpreal::default_rnd);

const mpreal pow(const int a, const unsigned long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const int a, const unsigned int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const int a, const long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const int a, const int b, mp_rnd_t rnd_mode = mpreal::default_rnd); 
const mpreal pow(const int a, const long double b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const int a, const double b, mp_rnd_t rnd_mode = mpreal::default_rnd); 

const mpreal pow(const long double a, const long double b, mp_rnd_t rnd_mode = mpreal::default_rnd);	
const mpreal pow(const long double a, const unsigned long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long double a, const unsigned int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long double a, const long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const long double a, const int b, mp_rnd_t rnd_mode = mpreal::default_rnd);

const mpreal pow(const double a, const double b, mp_rnd_t rnd_mode = mpreal::default_rnd);	
const mpreal pow(const double a, const unsigned long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const double a, const unsigned int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const double a, const long int b, mp_rnd_t rnd_mode = mpreal::default_rnd);
const mpreal pow(const double a, const int b, mp_rnd_t rnd_mode = mpreal::default_rnd);

//////////////////////////////////////////////////////////////////////////
// Estimate machine epsilon for the given precision
inline const mpreal machine_epsilon(mp_prec_t prec);
inline const mpreal mpreal_min(mp_prec_t prec);
inline const mpreal mpreal_max(mp_prec_t prec);

//////////////////////////////////////////////////////////////////////////
// Implementation of inline functions
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// Operators - Assignment
inline mpreal& mpreal::operator=(const mpreal& v)
{
	if (this!= &v)	mpfr_set(mp,v.mp,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator=(const mpf_t v)
{
	mpfr_set_f(mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator=(const mpz_t v)
{
	mpfr_set_z(mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator=(const mpq_t v)
{
	mpfr_set_q(mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator=(const long double v)		
{	
    mpfr_set_ld(mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator=(const double v)				
{	
    if(double_bits == -1 || fits_in_bits(v, double_bits))
    {
    	mpfr_set_d(mp,v,default_rnd);
    }
    else
        throw conversion_overflow();

	return *this;
}

inline mpreal& mpreal::operator=(const unsigned long int v)	
{	
	mpfr_set_ui(mp,v,default_rnd);	
	return *this;
}

inline mpreal& mpreal::operator=(const unsigned int v)		
{	
	mpfr_set_ui(mp,v,default_rnd);	
	return *this;
}

inline mpreal& mpreal::operator=(const long int v)			
{	
	mpfr_set_si(mp,v,default_rnd);	
	return *this;
}

inline mpreal& mpreal::operator=(const int v)
{	
	mpfr_set_si(mp,v,default_rnd);	
	return *this;
}

//////////////////////////////////////////////////////////////////////////
// + Addition
inline mpreal& mpreal::operator+=(const mpreal& v)
{
	mpfr_add(mp,mp,v.mp,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator+=(const mpf_t u)
{
	*this += mpreal(u);
	return *this;
}

inline mpreal& mpreal::operator+=(const mpz_t u)
{
	mpfr_add_z(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator+=(const mpq_t u)
{
	mpfr_add_q(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator+= (const long double u)
{
	return *this += mpreal(u);	
}

inline mpreal& mpreal::operator+= (const double u)
{
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
	mpfr_add_d(mp,mp,u,default_rnd);
	return *this;
#else
	return *this += mpreal(u);	
#endif
}

inline mpreal& mpreal::operator+=(const unsigned long int u)
{
	mpfr_add_ui(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator+=(const unsigned int u)
{
	mpfr_add_ui(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator+=(const long int u)
{
	mpfr_add_si(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator+=(const int u)
{
	mpfr_add_si(mp,mp,u,default_rnd);
	return *this;
}

inline const mpreal mpreal::operator+()const
{
	return mpreal(*this);
}

inline const mpreal operator+(const mpreal& a, const mpreal& b)
{
	// prec(a+b) = max(prec(a),prec(b))
	if(a.get_prec()>b.get_prec()) return mpreal(a) += b;
	else						  return mpreal(b) += a;
}

inline const std::string operator+(const mpreal& a, const std::string b)
{
	return (std::string)a+b;
}

inline const std::string operator+(const std::string a, const mpreal& b)
{
	return a+(std::string)b;
}

inline const mpreal operator+(const mpreal& a, const mpz_t b)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpreal& a, const char* b)
{
	return a+mpreal(b);
}

inline const mpreal operator+(const char* a, const mpreal& b)
{
	return mpreal(a)+b;

}

inline const mpreal operator+(const mpreal& a, const mpq_t b)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpreal& a, const long double b)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpreal& a, const double b)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpreal& a, const unsigned long int b)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpreal& a, const unsigned int b)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpreal& a, const long int b)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpreal& a, const int b)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpz_t b, const mpreal& a)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const mpq_t b, const mpreal& a)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const long double b, const mpreal& a)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const double  b, const mpreal& a)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const unsigned long int b, const mpreal& a)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const unsigned int b, const mpreal& a)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const long int b, const mpreal& a)
{
	return mpreal(a) += b;
}

inline const mpreal operator+(const int b, const mpreal& a)
{
	return mpreal(a) += b;
}

inline mpreal& mpreal::operator++() 
{
	*this += 1;
	return *this;
}

inline const mpreal mpreal::operator++ (int)
{
	mpreal x(*this);
	*this += 1;
	return x;
}

inline mpreal& mpreal::operator--() 
{
	*this -= 1;
	return *this;
}

inline const mpreal mpreal::operator-- (int)
{
	mpreal x(*this);
	*this -= 1;
	return x;
}

//////////////////////////////////////////////////////////////////////////
// - Subtraction
inline mpreal& mpreal::operator-= (const mpreal& v)
{
	mpfr_sub(mp,mp,v.mp,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator-=(const mpz_t v)
{
	mpfr_sub_z(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator-=(const mpq_t v)
{
	mpfr_sub_q(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator-=(const long double v)
{
	return *this -= mpreal(v);	
}

inline mpreal& mpreal::operator-=(const double v)
{
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
	mpfr_sub_d(mp,mp,v,default_rnd);
	return *this;
#else
	return *this -= mpreal(v);	
#endif
}

inline mpreal& mpreal::operator-=(const unsigned long int v)
{
	mpfr_sub_ui(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator-=(const unsigned int v)
{
	mpfr_sub_ui(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator-=(const long int v)
{
	mpfr_sub_si(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator-=(const int v)
{
	mpfr_sub_si(mp,mp,v,default_rnd);
	return *this;
}

inline const mpreal mpreal::operator-()const
{
	mpreal u(*this);
	mpfr_neg(u.mp,u.mp,default_rnd);
	return u;
}

inline const mpreal operator-(const mpreal& a, const mpreal& b)
{
	// prec(a-b) = max(prec(a),prec(b))
	if(a.get_prec()>b.get_prec())	return   mpreal(a) -= b;
	else							return -(mpreal(b) -= a);		
}

inline const mpreal operator-(const mpreal& a, const mpz_t b)
{
	return mpreal(a) -= b;
}

inline const mpreal operator-(const mpreal& a, const mpq_t b)
{
	return mpreal(a) -= b;
}

inline const mpreal operator-(const mpreal& a, const long double b)
{
	return mpreal(a) -= b;
}

inline const mpreal operator-(const mpreal& a, const double b)
{
	return mpreal(a) -= b;
}

inline const mpreal operator-(const mpreal& a, const unsigned long int b)
{
	return mpreal(a) -= b;
}

inline const mpreal operator-(const mpreal& a, const unsigned int b)
{
	return mpreal(a) -= b;
}

inline const mpreal operator-(const mpreal& a, const long int b)
{
	return mpreal(a) -= b;
}

inline const mpreal operator-(const mpreal& a, const int b)
{
	return mpreal(a) -= b;
}

inline const mpreal operator-(const mpz_t b, const mpreal& a)
{
	return -(mpreal(a) -= b);
}

inline const mpreal operator-(const mpq_t b, const mpreal& a)
{
	return -(mpreal(a) -= b);
}

inline const mpreal operator-(const long double b, const mpreal& a)
{
	return -(mpreal(a) -= b);
}

inline const mpreal operator-(const double  b, const mpreal& a)
{
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
	mpreal x(a);
	mpfr_d_sub(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
#else
	return -(mpreal(a) -= b);
#endif
}

inline const mpreal operator-(const unsigned long int b, const mpreal& a)
{
	mpreal x(a);
	mpfr_ui_sub(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
}

inline const mpreal operator-(const unsigned int b, const mpreal& a)
{
	mpreal x(a);
	mpfr_ui_sub(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
}

inline const mpreal operator-(const long int b, const mpreal& a)
{
	mpreal x(a);
	mpfr_si_sub(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
}

inline const mpreal operator-(const int b, const mpreal& a)
{
	mpreal x(a);
	mpfr_si_sub(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
}

inline const mpreal operator-(const mpreal& a, const char* b)
{
	return a-mpreal(b);
}

inline const mpreal operator-(const char* a, const mpreal& b)
{
	return mpreal(a)-b;
}

//////////////////////////////////////////////////////////////////////////
// * Multiplication
inline mpreal& mpreal::operator*= (const mpreal& v)
{
	mpfr_mul(mp,mp,v.mp,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator*=(const mpz_t v)
{
	mpfr_mul_z(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator*=(const mpq_t v)
{
	mpfr_mul_q(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator*=(const long double v)
{
	return *this *= mpreal(v);	
}

inline mpreal& mpreal::operator*=(const double v)
{
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
	mpfr_mul_d(mp,mp,v,default_rnd);
	return *this;
#else
	return *this *= mpreal(v);	
#endif
}

inline mpreal& mpreal::operator*=(const unsigned long int v)
{
	mpfr_mul_ui(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator*=(const unsigned int v)
{
	mpfr_mul_ui(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator*=(const long int v)
{
	mpfr_mul_si(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator*=(const int v)
{
	mpfr_mul_si(mp,mp,v,default_rnd);
	return *this;
}

inline const mpreal operator*(const mpreal& a, const mpreal& b)
{
	// prec(a*b) = max(prec(a),prec(b))
	if(a.get_prec()>b.get_prec())	return   mpreal(a) *= b;
	else							return   mpreal(b) *= a;		
}

inline const mpreal operator*(const mpreal& a, const mpz_t b)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpreal& a, const mpq_t b)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpreal& a, const long double b)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpreal& a, const double b)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpreal& a, const unsigned long int b)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpreal& a, const unsigned int b)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpreal& a, const long int b)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpreal& a, const int b)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpz_t b, const mpreal& a)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const mpq_t b, const mpreal& a)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const long double b, const mpreal& a)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const double  b, const mpreal& a)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const unsigned long int b, const mpreal& a)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const unsigned int b, const mpreal& a)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const long int b, const mpreal& a)
{
	return mpreal(a) *= b;
}

inline const mpreal operator*(const int b, const mpreal& a)
{
	return mpreal(a) *= b;
}

//////////////////////////////////////////////////////////////////////////
// / Division
inline mpreal& mpreal::operator/=(const mpreal& v)
{
	mpfr_div(mp,mp,v.mp,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator/=(const mpz_t v)
{
	mpfr_div_z(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator/=(const mpq_t v)
{
	mpfr_div_q(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator/=(const long double v)
{
	return *this /= mpreal(v);	
}

inline mpreal& mpreal::operator/=(const double v)
{
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
	mpfr_div_d(mp,mp,v,default_rnd);
	return *this;
#else
	return *this /= mpreal(v);	
#endif
}

inline mpreal& mpreal::operator/=(const unsigned long int v)
{
	mpfr_div_ui(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator/=(const unsigned int v)
{
	mpfr_div_ui(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator/=(const long int v)
{
	mpfr_div_si(mp,mp,v,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator/=(const int v)
{
	mpfr_div_si(mp,mp,v,default_rnd);
	return *this;
}

inline const mpreal operator/(const mpreal& a, const mpreal& b)
{
	mpreal x(a);
	mp_prec_t pb;
	mp_prec_t pa;

	// prec(a/b) = max(prec(a),prec(b))
	pa = a.get_prec();
	pb = b.get_prec();
	if(pb>pa) x.set_prec(pb);

	return   x /= b;
}

inline const mpreal operator/(const mpreal& a, const mpz_t b)
{
	return mpreal(a) /= b;
}

inline const mpreal operator/(const mpreal& a, const mpq_t b)
{
	return mpreal(a) /= b;
}

inline const mpreal operator/(const mpreal& a, const long double b)
{
	return mpreal(a) /= b;
}

inline const mpreal operator/(const mpreal& a, const double b)
{
	return mpreal(a) /= b;
}

inline const mpreal operator/(const mpreal& a, const unsigned long int b)
{
	return mpreal(a) /= b;
}

inline const mpreal operator/(const mpreal& a, const unsigned int b)
{
	return mpreal(a) /= b;
}

inline const mpreal operator/(const mpreal& a, const long int b)
{
	return mpreal(a) /= b;
}

inline const mpreal operator/(const mpreal& a, const int b)
{
	return mpreal(a) /= b;
}

inline const mpreal operator/(const unsigned long int b, const mpreal& a)
{
	mpreal x(a);
	mpfr_ui_div(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
}

inline const mpreal operator/(const unsigned int b, const mpreal& a)
{
	mpreal x(a);
	mpfr_ui_div(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
}

inline const mpreal operator/(const long int b, const mpreal& a)
{
	mpreal x(a);
	mpfr_si_div(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
}

inline const mpreal operator/(const int b, const mpreal& a)
{
	mpreal x(a);
	mpfr_si_div(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
}

inline const mpreal operator/(const long double b, const mpreal& a)
{
	mpreal x(b);
	return x/a;
}

inline const mpreal operator/(const double  b, const mpreal& a)
{
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
	mpreal x(a);
	mpfr_d_div(x.mp,b,a.mp,mpreal::default_rnd);
	return x;
#else
	mpreal x(b);
	return x/a;
#endif
}

//////////////////////////////////////////////////////////////////////////
// Shifts operators - Multiplication/Division by power of 2
inline mpreal& mpreal::operator<<=(const unsigned long int u)
{
	mpfr_mul_2ui(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator<<=(const unsigned int u)
{
	mpfr_mul_2ui(mp,mp,static_cast<unsigned long int>(u),default_rnd);
	return *this;
}

inline mpreal& mpreal::operator<<=(const long int u)
{
	mpfr_mul_2si(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator<<=(const int u)
{
	mpfr_mul_2si(mp,mp,static_cast<long int>(u),default_rnd);
	return *this;
}

inline mpreal& mpreal::operator>>=(const unsigned long int u)
{
	mpfr_div_2ui(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator>>=(const unsigned int u)
{
	mpfr_div_2ui(mp,mp,static_cast<unsigned long int>(u),default_rnd);
	return *this;
}

inline mpreal& mpreal::operator>>=(const long int u)
{
	mpfr_div_2si(mp,mp,u,default_rnd);
	return *this;
}

inline mpreal& mpreal::operator>>=(const int u)
{
	mpfr_div_2si(mp,mp,static_cast<long int>(u),default_rnd);
	return *this;
}

inline const mpreal operator<<(const mpreal& v, const unsigned long int k)
{
	return mul_2ui(v,k);
}

inline const mpreal operator<<(const mpreal& v, const unsigned int k)
{
	return mul_2ui(v,static_cast<unsigned long int>(k));
}

inline const mpreal operator<<(const mpreal& v, const long int k)
{
	return mul_2si(v,k);
}

inline const mpreal operator<<(const mpreal& v, const int k)
{
	return mul_2si(v,static_cast<long int>(k));
}

inline const mpreal operator>>(const mpreal& v, const unsigned long int k)
{
	return div_2ui(v,k);
}

inline const mpreal operator>>(const mpreal& v, const long int k)
{
	return div_2si(v,k);
}

inline const mpreal operator>>(const mpreal& v, const unsigned int k)
{
	return div_2ui(v,static_cast<unsigned long int>(k));
}

inline const mpreal operator>>(const mpreal& v, const int k)
{
	return div_2si(v,static_cast<long int>(k));
}

// mul_2ui
inline const mpreal mul_2ui(const mpreal& v, unsigned long int k, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_mul_2ui(x.mp,v.mp,k,rnd_mode);
	return x;
}

// mul_2si
inline const mpreal mul_2si(const mpreal& v, long int k, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_mul_2si(x.mp,v.mp,k,rnd_mode);
	return x;
}

inline const mpreal div_2ui(const mpreal& v, unsigned long int k, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_div_2ui(x.mp,v.mp,k,rnd_mode);
	return x;
}

inline const mpreal div_2si(const mpreal& v, long int k, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_div_2si(x.mp,v.mp,k,rnd_mode);
	return x;
}

//////////////////////////////////////////////////////////////////////////
//Boolean operators
inline bool operator > (const mpreal& a, const mpreal& b)
{
	return (mpfr_greater_p(a.mp,b.mp)!=0);
}

inline bool operator >  (const mpreal& a, const unsigned long int b)
{
	return a>mpreal(b);
}

inline bool operator >  (const mpreal& a, const unsigned int b)
{
	return a>mpreal(b);
}

inline bool operator >  (const mpreal& a, const long int b)
{
	return a>mpreal(b);
}

inline bool operator >  (const mpreal& a, const int b)
{
	return a>mpreal(b);
}

inline bool operator >  (const mpreal& a, const long double b)
{
	return a>mpreal(b);
}

inline bool operator >  (const mpreal& a, const double b)
{
	return a>mpreal(b);
}

inline bool operator >  (const unsigned long int a,	const mpreal& b)
{
	return mpreal(a)>b;
}

inline bool operator >  (const unsigned int a,		const mpreal& b)
{
	return mpreal(a)>b;
}

inline bool operator >  (const long int a,			const mpreal& b)
{
	return mpreal(a)>b;
}

inline bool operator >  (const int a,				const mpreal& b)
{
	return mpreal(a)>b;
}

inline bool operator >  (const long double a,		const mpreal& b)
{
	return mpreal(a)>b;
}

inline bool operator >  (const double a,			const mpreal& b)
{
	return mpreal(a)>b;
}

inline bool operator >= (const mpreal& a, const mpreal& b)
{
	return (mpfr_greaterequal_p(a.mp,b.mp)!=0);
}

inline bool operator >=  (const mpreal& a, const unsigned long int b)
{
	return a>=mpreal(b);
}

inline bool operator >=  (const mpreal& a, const unsigned int b)
{
	return a>=mpreal(b);
}

inline bool operator >=  (const mpreal& a, const long int b)
{
	return a>=mpreal(b);
}

inline bool operator >=  (const mpreal& a, const int b)
{
	return a>=mpreal(b);
}

inline bool operator >=  (const mpreal& a, const long double b)
{
	return a>=mpreal(b);
}

inline bool operator >=  (const mpreal& a, const double b)
{
	return a>=mpreal(b);
}

inline bool operator >=  (const unsigned long int a,const mpreal& b)
{
	return mpreal(a)>=b;
}

inline bool operator >=  (const unsigned int a,		const mpreal& b)
{
	return mpreal(a)>=b;
}

inline bool operator >=  (const long int a,			const mpreal& b)
{
	return mpreal(a)>=b;
}

inline bool operator >=  (const int a,				const mpreal& b)
{
	return mpreal(a)>=b;
}

inline bool operator >=  (const long double a,		const mpreal& b)
{
	return mpreal(a)>=b;
}

inline bool operator >=  (const double a,			const mpreal& b)
{
	return mpreal(a)>=b;
}

inline bool operator <  (const mpreal& a, const mpreal& b)
{
	return (mpfr_less_p(a.mp,b.mp)!=0);
}

inline bool operator <  (const mpreal& a, const unsigned long int b)
{
	return a<mpreal(b);
}

inline bool operator <  (const mpreal& a, const unsigned int b)
{
	return a<mpreal(b);
}

inline bool operator <  (const mpreal& a, const long int b)
{
	return a<mpreal(b);
}

inline bool operator <  (const mpreal& a, const int b)
{
	return a<mpreal(b);
}

inline bool operator <  (const mpreal& a, const long double b)
{
	return a<mpreal(b);
}

inline bool operator <  (const mpreal& a, const double b)
{
	return a<mpreal(b);
}

inline bool operator <  (const unsigned long int a,	const mpreal& b)
{
	return mpreal(a)<b;
}

inline bool operator <  (const unsigned int a,const mpreal& b)
{
	return mpreal(a)<b;
}

inline bool operator <  (const long int a,const mpreal& b)
{
	return mpreal(a)<b;
}

inline bool operator <  (const int a,const mpreal& b)
{
	return mpreal(a)<b;
}

inline bool operator <  (const long double a,const mpreal& b)
{
	return mpreal(a)<b;
}

inline bool operator <  (const double a,const mpreal& b)
{
	return mpreal(a)<b;
}

inline bool operator <= (const mpreal& a, const mpreal& b)
{
	return (mpfr_lessequal_p(a.mp,b.mp)!=0);
}

inline bool operator <=  (const mpreal& a, const unsigned long int b)
{
	return a<=mpreal(b);
}

inline bool operator <=  (const mpreal& a, const unsigned int b)
{
	return a<=mpreal(b);
}

inline bool operator <=  (const mpreal& a, const long int b)
{
	return a<=mpreal(b);
}

inline bool operator <=  (const mpreal& a, const int b)
{
	return a<=mpreal(b);
}

inline bool operator <=  (const mpreal& a, const long double b)
{
	return a<=mpreal(b);
}

inline bool operator <=  (const mpreal& a, const double b)
{
	return a<=mpreal(b);
}

inline bool operator <=  (const unsigned long int a,const mpreal& b)
{
	return mpreal(a)<=b;
}

inline bool operator <=  (const unsigned int a,		const mpreal& b)
{
	return mpreal(a)<=b;
}

inline bool operator <=  (const long int a,			const mpreal& b)
{
	return mpreal(a)<=b;
}

inline bool operator <=  (const int a,				const mpreal& b)
{
	return mpreal(a)<=b;
}

inline bool operator <=  (const long double a,		const mpreal& b)
{
	return mpreal(a)<=b;
}

inline bool operator <=  (const double a,			const mpreal& b)
{
	return mpreal(a)<=b;
}

inline bool operator == (const mpreal& a, const mpreal& b)
{
	return (mpfr_equal_p(a.mp,b.mp)!=0);
}

inline bool operator ==  (const mpreal& a, const unsigned long int b)
{
	return a==mpreal(b);
}

inline bool operator ==  (const mpreal& a, const unsigned int b)
{
	return a==mpreal(b);
}

inline bool operator ==  (const mpreal& a, const long int b)
{
	return a==mpreal(b);
}

inline bool operator ==  (const mpreal& a, const int b)
{
	return a==mpreal(b);
}

inline bool operator ==  (const mpreal& a, const long double b)
{
	return a==mpreal(b);
}

inline bool operator ==  (const mpreal& a, const double b)
{
	return a==mpreal(b);
}

inline bool operator ==  (const unsigned long int a,const mpreal& b)
{
	return mpreal(a)==b;
}

inline bool operator ==  (const unsigned int a,		const mpreal& b)
{
	return mpreal(a)==b;
}

inline bool operator ==  (const long int a,			const mpreal& b)
{
	return mpreal(a)==b;
}

inline bool operator ==  (const int a,				const mpreal& b)
{
	return mpreal(a)==b;
}

inline bool operator ==  (const long double a,		const mpreal& b)
{
	return mpreal(a)==b;
}

inline bool operator ==  (const double a,			const mpreal& b)
{
	return mpreal(a)==b;
}

inline bool operator != (const mpreal& a, const mpreal& b)
{
	return (mpfr_lessgreater_p(a.mp,b.mp)!=0);
}

inline bool operator !=  (const mpreal& a, const unsigned long int b)
{
	return a!=mpreal(b);
}

inline bool operator !=  (const mpreal& a, const unsigned int b)
{
	return a!=mpreal(b);
}

inline bool operator !=  (const mpreal& a, const long int b)
{
	return a!=mpreal(b);
}

inline bool operator !=  (const mpreal& a, const int b)
{
	return a!=mpreal(b);
}

inline bool operator !=  (const mpreal& a, const long double b)
{
	return a!=mpreal(b);
}

inline bool operator !=  (const mpreal& a, const double b)
{
	return a!=mpreal(b);
}

inline bool operator !=  (const unsigned long int a,const mpreal& b)
{
	return mpreal(a)!=b;
}

inline bool operator !=  (const unsigned int a,		const mpreal& b)
{
	return mpreal(a)!=b;
}

inline bool operator !=  (const long int a,			const mpreal& b)
{
	return mpreal(a)!=b;
}

inline bool operator !=  (const int a,				const mpreal& b)
{
	return mpreal(a)!=b;
}

inline bool operator !=  (const long double a,		const mpreal& b)
{
	return mpreal(a)!=b;
}

inline bool operator !=  (const double a,			const mpreal& b)
{
	return mpreal(a)!=b;
}

inline bool _isnan(const mpreal& v)
{
	return (mpfr_nan_p(v.mp)!=0);
}

inline bool _isinf(const mpreal& v)
{
	return (mpfr_inf_p(v.mp)!=0);
}

inline bool _isnum(const mpreal& v)
{
	return (mpfr_number_p(v.mp)!=0);
}

inline bool _iszero(const mpreal& v)
{
	return (mpfr_zero_p(v.mp)!=0);
}

inline bool _isint(const mpreal& v)
{
	return (mpfr_integer_p(v.mp)!=0);
}

#if (MPFR_VERSION >= MPFR_VERSION_NUM(3,0,0))
inline bool _isregular(const mpreal& v)
{
	return (mpfr_regular_p(v.mp));
}
#endif // MPFR 3.0.0 Specifics

//////////////////////////////////////////////////////////////////////////
// Type Converters
inline mpreal::operator double() const
{
	return mpfr_get_d(mp,default_rnd);
}

inline mpreal::operator float() const
{
	return (float)mpfr_get_d(mp,default_rnd);
}

inline mpreal::operator long double() const
{
	return mpfr_get_ld(mp,default_rnd);
}

inline mpreal::operator unsigned long() const
{
	return mpfr_get_ui(mp,default_rnd);	
}

inline mpreal::operator unsigned int() const
{
	return static_cast<unsigned int>(mpfr_get_ui(mp,default_rnd));	
}

inline mpreal::operator long() const
{
	return mpfr_get_si(mp,default_rnd);	
}

inline mpreal::operator mpfr_ptr()
{
	return mp;
}

//////////////////////////////////////////////////////////////////////////
// Set/Get number properties
inline int sgn(const mpreal& v)
{
	int r = mpfr_signbit(v.mp);
	return (r>0?-1:1);
}

inline void mpreal::set_sign(int sign, mp_rnd_t rnd_mode)
{
	mpfr_setsign(mp,mp,(sign<0?1:0),rnd_mode);
}

inline mp_prec_t mpreal::get_prec() const
{
	return mpfr_get_prec(mp);
}

inline void mpreal::set_prec(mp_prec_t prec, mp_rnd_t rnd_mode)
{
	mpfr_prec_round(mp,prec,rnd_mode);
}

inline void mpreal::set_inf(int sign) 
{ 
	mpfr_set_inf(mp,sign);
}	

inline void mpreal::set_nan() 
{
	mpfr_set_nan(mp);
}

inline mp_exp_t mpreal::get_exp ()
{
	return mpfr_get_exp(mp);
}

inline int mpreal::set_exp (mp_exp_t e)
{
	return mpfr_set_exp(mp,e);
}

inline const mpreal frexp(const mpreal& v, mp_exp_t* exp)
{
	mpreal x(v);
	*exp = x.get_exp();
	x.set_exp(0);
	return x;
}

inline const mpreal ldexp(const mpreal& v, mp_exp_t exp)
{
	mpreal x(v);

	// rounding is not important since we just increasing the exponent
	mpfr_mul_2si(x.mp,x.mp,exp,mpreal::default_rnd); 
	return x;
}

inline const mpreal machine_epsilon(mp_prec_t prec)
{
	// smallest eps such that 1.0+eps != 1.0
	// depends (of cause) on the precision
	mpreal x(1,prec); 
	return nextabove(x)-x;
}

inline const mpreal mpreal_min(mp_prec_t prec)
{
	// min = 1/2*2^emin = 2^(emin-1)
	
	mpreal x(1,prec);
	return x <<= mpreal::get_emin()-1;
}

inline const mpreal mpreal_max(mp_prec_t prec)
{
	// max = (1-eps)*2^emax, assume eps = 0?, 
	// and use emax-1 to prevent value to be +inf
	// max = 2^(emax-1)

	mpreal x(1,prec);
	return x <<= mpreal::get_emax()-1;
}

inline const mpreal modf(const mpreal& v, mpreal& n)
{
	mpreal frac(v);

	// rounding is not important since we are using the same number
	mpfr_frac(frac.mp,frac.mp,mpreal::default_rnd);	
	mpfr_trunc(n.mp,v.mp);
	return frac;
}

inline int mpreal::check_range (int t, mp_rnd_t rnd_mode)
{
	return mpfr_check_range(mp,t,rnd_mode);
}

inline int mpreal::subnormalize (int t,mp_rnd_t rnd_mode)
{
	return mpfr_subnormalize(mp,t,rnd_mode);
}

inline mp_exp_t mpreal::get_emin (void)
{
	return mpfr_get_emin();
}

inline int mpreal::set_emin (mp_exp_t exp)
{
	return mpfr_set_emin(exp);
}

inline mp_exp_t mpreal::get_emax (void)
{
	return mpfr_get_emax();
}

inline int mpreal::set_emax (mp_exp_t exp)
{
	return mpfr_set_emax(exp);
}

inline mp_exp_t mpreal::get_emin_min (void)
{
	return mpfr_get_emin_min();
}

inline mp_exp_t mpreal::get_emin_max (void)
{
	return mpfr_get_emin_max();
}

inline mp_exp_t mpreal::get_emax_min (void)
{
	return mpfr_get_emax_min();
}

inline mp_exp_t mpreal::get_emax_max (void)
{
	return mpfr_get_emax_max();
}

//////////////////////////////////////////////////////////////////////////
// Mathematical Functions
//////////////////////////////////////////////////////////////////////////
inline const mpreal sqr(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_sqr(x.mp,x.mp,rnd_mode);
	return x;
}

inline const mpreal sqrt(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_sqrt(x.mp,x.mp,rnd_mode);
	return x;
}

inline const mpreal sqrt(const unsigned long int v, mp_rnd_t rnd_mode)
{
	mpreal x;
	mpfr_sqrt_ui(x.mp,v,rnd_mode);
	return x;
}

inline const mpreal sqrt(const unsigned int v, mp_rnd_t rnd_mode)
{
	return sqrt(static_cast<unsigned long int>(v),rnd_mode);
}

inline const mpreal sqrt(const long int v, mp_rnd_t rnd_mode)
{
	if (v>=0)	return sqrt(static_cast<unsigned long int>(v),rnd_mode);
	else		return mpreal(); // NaN  
}

inline const mpreal sqrt(const int v, mp_rnd_t rnd_mode)
{
	if (v>=0)	return sqrt(static_cast<unsigned long int>(v),rnd_mode);
	else		return mpreal(); // NaN
}

inline const mpreal sqrt(const long double v, mp_rnd_t rnd_mode)
{
	return sqrt(mpreal(v),rnd_mode);
}

inline const mpreal sqrt(const double v, mp_rnd_t rnd_mode)
{
	return sqrt(mpreal(v),rnd_mode);
}

inline const mpreal cbrt(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_cbrt(x.mp,x.mp,rnd_mode);
	return x;
}

inline const mpreal root(const mpreal& v, unsigned long int k, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_root(x.mp,x.mp,k,rnd_mode);
	return x;
}

inline const mpreal fabs(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_abs(x.mp,x.mp,rnd_mode);
	return x;
}

inline const mpreal abs(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_abs(x.mp,x.mp,rnd_mode);
	return x;
}

inline const mpreal dim(const mpreal& a, const mpreal& b, mp_rnd_t rnd_mode)
{
	mpreal x(a);
	mpfr_dim(x.mp,a.mp,b.mp,rnd_mode);
	return x;
}

inline int cmpabs(const mpreal& a,const mpreal& b)
{
	return mpfr_cmpabs(a.mp,b.mp);
}

inline const mpreal log  (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_log(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal log2(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_log2(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal log10(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_log10(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal exp(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_exp(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal exp2(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_exp2(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal exp10(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_exp10(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal cos(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_cos(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal sin(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_sin(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal tan(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_tan(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal sec(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_sec(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal csc(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_csc(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal cot(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_cot(x.mp,v.mp,rnd_mode);
	return x;
}

inline int sin_cos(mpreal& s, mpreal& c, const mpreal& v, mp_rnd_t rnd_mode)
{
	return mpfr_sin_cos(s.mp,c.mp,v.mp,rnd_mode);
}

inline const mpreal acos (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_acos(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal asin (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_asin(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal atan (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_atan(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal atan2 (const mpreal& y, const mpreal& x, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t yp, xp;

	yp = y.get_prec(); 
	xp = x.get_prec(); 

	a.set_prec(yp>xp?yp:xp);

	mpfr_atan2(a.mp, y.mp, x.mp, rnd_mode);

	return a;
}

inline const mpreal cosh (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_cosh(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal sinh (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_sinh(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal tanh (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_tanh(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal sech (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_sech(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal csch (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_csch(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal coth (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_coth(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal acosh  (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_acosh(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal asinh  (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_asinh(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal atanh  (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_atanh(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal fac_ui (unsigned long int v, mp_rnd_t rnd_mode)
{
	mpreal x;
	mpfr_fac_ui(x.mp,v,rnd_mode);
	return x;
}

inline const mpreal log1p  (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_log1p(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal expm1  (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_expm1(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal eint   (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_eint(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal gamma (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_gamma(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal lngamma (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_lngamma(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal lgamma (const mpreal& v, int *signp, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_lgamma(x.mp,signp,v.mp,rnd_mode);
	return x;
}

inline const mpreal zeta (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_zeta(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal erf (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_erf(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal erfc (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_erfc(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal _j0 (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_j0(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal _j1 (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_j1(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal _jn (long n, const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_jn(x.mp,n,v.mp,rnd_mode);
	return x;
}

inline const mpreal _y0 (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_y0(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal _y1 (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_y1(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal _yn (long n, const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_yn(x.mp,n,v.mp,rnd_mode);
	return x;
}

//////////////////////////////////////////////////////////////////////////
// MPFR 2.4.0 Specifics
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))

inline int sinh_cosh(mpreal& s, mpreal& c, const mpreal& v, mp_rnd_t rnd_mode)
{
	return mpfr_sinh_cosh(s.mp,c.mp,v.mp,rnd_mode);
}

inline const mpreal li2(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_li2(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal fmod (const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t yp, xp;

	yp = y.get_prec(); 
	xp = x.get_prec(); 

	a.set_prec(yp>xp?yp:xp);

	mpfr_fmod(a.mp, x.mp, y.mp, rnd_mode);

	return a;
}

inline const mpreal rec_sqrt(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_rec_sqrt(x.mp,v.mp,rnd_mode);
	return x;
}
#endif //  MPFR 2.4.0 Specifics

//////////////////////////////////////////////////////////////////////////
// MPFR 3.0.0 Specifics
#if (MPFR_VERSION >= MPFR_VERSION_NUM(3,0,0))
inline const mpreal digamma(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_digamma(x.mp,v.mp,rnd_mode);
	return x;
}
#endif // MPFR 3.0.0 Specifics

//////////////////////////////////////////////////////////////////////////
// Constants
inline const mpreal const_log2 (mp_prec_t prec, mp_rnd_t rnd_mode)
{
	mpreal x;
	x.set_prec(prec);
	mpfr_const_log2(x.mp,rnd_mode);
	return x;
}

inline const mpreal const_pi (mp_prec_t prec, mp_rnd_t rnd_mode)
{
	mpreal x;
	x.set_prec(prec);
	mpfr_const_pi(x.mp,rnd_mode);
	return x;
}

inline const mpreal const_euler (mp_prec_t prec, mp_rnd_t rnd_mode)
{
	mpreal x;
	x.set_prec(prec);
	mpfr_const_euler(x.mp,rnd_mode);
	return x;
}

inline const mpreal const_catalan (mp_prec_t prec, mp_rnd_t rnd_mode)
{
	mpreal x;
	x.set_prec(prec);
	mpfr_const_catalan(x.mp,rnd_mode);
	return x;
}

inline const mpreal const_infinity (int sign, mp_prec_t prec, mp_rnd_t rnd_mode)
{
	mpreal x;
	x.set_prec(prec,rnd_mode);
	mpfr_set_inf(x.mp, sign);
	return x;
}

//////////////////////////////////////////////////////////////////////////
// Integer Related Functions
inline const mpreal rint(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_rint(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal ceil(const mpreal& v)
{
	mpreal x(v);
	mpfr_ceil(x.mp,v.mp);
	return x;

}

inline const mpreal floor(const mpreal& v)
{
	mpreal x(v);
	mpfr_floor(x.mp,v.mp);
	return x;
}

inline const mpreal round(const mpreal& v)
{
	mpreal x(v);
	mpfr_round(x.mp,v.mp);
	return x;
}

inline const mpreal trunc(const mpreal& v)
{
	mpreal x(v);
	mpfr_trunc(x.mp,v.mp);
	return x;
}

inline const mpreal rint_ceil (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_rint_ceil(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal rint_floor(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_rint_floor(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal rint_round(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_rint_round(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal rint_trunc(const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_rint_trunc(x.mp,v.mp,rnd_mode);
	return x;
}

inline const mpreal frac (const mpreal& v, mp_rnd_t rnd_mode)
{
	mpreal x(v);
	mpfr_frac(x.mp,v.mp,rnd_mode);
	return x;
}

//////////////////////////////////////////////////////////////////////////
// Miscellaneous Functions
inline void swap(mpreal& a, mpreal& b) 
{
	mpfr_swap(a.mp,b.mp);
}

#ifndef max
inline const mpreal max(const mpreal& x, const mpreal& y)
{
	return (x>y?x:y);
}
#endif

#ifndef min
inline const mpreal min(const mpreal& x, const mpreal& y)
{
	return (x<y?x:y);
}
#endif

inline const mpreal nexttoward (const mpreal& x, const mpreal& y)
{
	mpreal a(x);
	mpfr_nexttoward(a.mp,y.mp);
	return a;
}

inline const mpreal nextabove  (const mpreal& x)
{
	mpreal a(x);
	mpfr_nextabove(a.mp);
	return a;
}

inline const mpreal nextbelow  (const mpreal& x)
{
	mpreal a(x);
	mpfr_nextbelow(a.mp);
	return a;
}

inline const mpreal urandomb (gmp_randstate_t& state)
{
	mpreal x;
	mpfr_urandomb(x.mp,state);
	return x;
}

#if (MPFR_VERSION >= MPFR_VERSION_NUM(3,0,0))
// use gmp_randinit_default() to init state, gmp_randclear() to clear
inline const mpreal urandom (gmp_randstate_t& state,mp_rnd_t rnd_mode)
{
	mpreal x;
	mpfr_urandom(x.mp,state,rnd_mode);
	return x;
}
#endif 

#if (MPFR_VERSION <= MPFR_VERSION_NUM(2,4,2))
inline const mpreal random2 (mp_size_t size, mp_exp_t exp)
{
	mpreal x;
	mpfr_random2(x.mp,size,exp);
	return x;
}
#endif

//////////////////////////////////////////////////////////////////////////
// Set/Get global properties
inline void mpreal::set_default_prec(mp_prec_t prec)
{ 
	default_prec = prec;
	mpfr_set_default_prec(prec); 
}

inline mp_prec_t mpreal::get_default_prec()
{ 
	return mpfr_get_default_prec();
}

inline void mpreal::set_default_base(int base)
{ 
	default_base = base;
}

inline int mpreal::get_default_base()
{ 
	return default_base;
}

inline void mpreal::set_default_rnd(mp_rnd_t rnd_mode)
{ 
	default_rnd =  rnd_mode;
	mpfr_set_default_rounding_mode(rnd_mode); 
}

inline mp_rnd_t mpreal::get_default_rnd()
{ 
	return mpfr_get_default_rounding_mode();
}

inline void mpreal::set_double_bits(int dbits)
{ 
	double_bits = dbits;
}

inline int mpreal::get_double_bits()
{ 
	return double_bits;
}

inline bool mpreal::fits_in_bits(double x, int n)
{   
	int i;
	double t;
	return IsInf(x) || (std::modf ( std::ldexp ( std::frexp ( x, &i ), n ), &t ) == 0.0);
}

inline const mpreal pow(const mpreal& a, const mpreal& b, mp_rnd_t rnd_mode)
{
	mpreal x(a);
	mpfr_pow(x.mp,x.mp,b.mp,rnd_mode);
	return x;
}

inline const mpreal pow(const mpreal& a, const mpz_t b, mp_rnd_t rnd_mode)
{
	mpreal x(a);
	mpfr_pow_z(x.mp,x.mp,b,rnd_mode);
	return x;
}

inline const mpreal pow(const mpreal& a, const unsigned long int b, mp_rnd_t rnd_mode)
{
	mpreal x(a);
	mpfr_pow_ui(x.mp,x.mp,b,rnd_mode);
	return x;
}

inline const mpreal pow(const mpreal& a, const unsigned int b, mp_rnd_t rnd_mode)
{
	return pow(a,static_cast<unsigned long int>(b),rnd_mode);
}

inline const mpreal pow(const mpreal& a, const long int b, mp_rnd_t rnd_mode)
{
	mpreal x(a);
	mpfr_pow_si(x.mp,x.mp,b,rnd_mode);
	return x;
}

inline const mpreal pow(const mpreal& a, const int b, mp_rnd_t rnd_mode)
{
	return pow(a,static_cast<long int>(b),rnd_mode);
}

inline const mpreal pow(const mpreal& a, const long double b, mp_rnd_t rnd_mode)
{
	return pow(a,mpreal(b),rnd_mode);
}

inline const mpreal pow(const mpreal& a, const double b, mp_rnd_t rnd_mode)
{
	return pow(a,mpreal(b),rnd_mode);
}

inline const mpreal pow(const unsigned long int a, const mpreal& b, mp_rnd_t rnd_mode)
{
	mpreal x(a);
	mpfr_ui_pow(x.mp,a,b.mp,rnd_mode);
	return x;
}

inline const mpreal pow(const unsigned int a, const mpreal& b, mp_rnd_t rnd_mode)
{
	return pow(static_cast<unsigned long int>(a),b,rnd_mode);
}

inline const mpreal pow(const long int a, const mpreal& b, mp_rnd_t rnd_mode)
{
	if (a>=0) 	return pow(static_cast<unsigned long int>(a),b,rnd_mode);
	else		return pow(mpreal(a),b,rnd_mode);
}

inline const mpreal pow(const int a, const mpreal& b, mp_rnd_t rnd_mode)
{
	if (a>=0) 	return pow(static_cast<unsigned long int>(a),b,rnd_mode);
	else		return pow(mpreal(a),b,rnd_mode);
}

inline const mpreal pow(const long double a, const mpreal& b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),b,rnd_mode);
}

inline const mpreal pow(const double a, const mpreal& b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),b,rnd_mode);
}

// pow unsigned long int
inline const mpreal pow(const unsigned long int a, const unsigned long int b, mp_rnd_t rnd_mode)
{
	mpreal x(a);
	mpfr_ui_pow_ui(x.mp,a,b,rnd_mode);
	return x;
}

inline const mpreal pow(const unsigned long int a, const unsigned int b, mp_rnd_t rnd_mode)
{
	return pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
}

inline const mpreal pow(const unsigned long int a, const long int b, mp_rnd_t rnd_mode)
{
	if(b>0)	return pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
	else	return pow(a,mpreal(b),rnd_mode); //mpfr_ui_pow
}

inline const mpreal pow(const unsigned long int a, const int b, mp_rnd_t rnd_mode)
{
	if(b>0)	return pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
	else	return pow(a,mpreal(b),rnd_mode); //mpfr_ui_pow
}

inline const mpreal pow(const unsigned long int a, const long double b, mp_rnd_t rnd_mode)
{
	return pow(a,mpreal(b),rnd_mode); //mpfr_ui_pow
}

inline const mpreal pow(const unsigned long int a, const double b, mp_rnd_t rnd_mode)
{
	return pow(a,mpreal(b),rnd_mode); //mpfr_ui_pow
}

// pow unsigned int
inline const mpreal pow(const unsigned int a, const unsigned long int b, mp_rnd_t rnd_mode)
{
	return pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
}

inline const mpreal pow(const unsigned int a, const unsigned int b, mp_rnd_t rnd_mode)
{
	return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
}

inline const mpreal pow(const unsigned int a, const long int b, mp_rnd_t rnd_mode)
{
	if(b>0)	return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
	else	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
}

inline const mpreal pow(const unsigned int a, const int b, mp_rnd_t rnd_mode)
{
	if(b>0)	return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
	else	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
}

inline const mpreal pow(const unsigned int a, const long double b, mp_rnd_t rnd_mode)
{
	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
}

inline const mpreal pow(const unsigned int a, const double b, mp_rnd_t rnd_mode)
{
	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
}

// pow long int
inline const mpreal pow(const long int a, const unsigned long int b, mp_rnd_t rnd_mode)
{
	if (a>0) return pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
	else	 return pow(mpreal(a),b,rnd_mode); //mpfr_pow_ui
}

inline const mpreal pow(const long int a, const unsigned int b, mp_rnd_t rnd_mode)
{
	if (a>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode);  //mpfr_ui_pow_ui
	else	 return pow(mpreal(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const mpreal pow(const long int a, const long int b, mp_rnd_t rnd_mode)
{
	if (a>0)
	{
		if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
		else	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
	}else{
		return pow(mpreal(a),b,rnd_mode); // mpfr_pow_si
	}
}

inline const mpreal pow(const long int a, const int b, mp_rnd_t rnd_mode)
{
	if (a>0)
	{
		if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
		else	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
	}else{
		return pow(mpreal(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
	}
}

inline const mpreal pow(const long int a, const long double b, mp_rnd_t rnd_mode)
{
	if (a>=0) 	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
	else		return pow(mpreal(a),mpreal(b),rnd_mode); //mpfr_pow
}

inline const mpreal pow(const long int a, const double b, mp_rnd_t rnd_mode)
{
	if (a>=0) 	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
	else		return pow(mpreal(a),mpreal(b),rnd_mode); //mpfr_pow
}

// pow int
inline const mpreal pow(const int a, const unsigned long int b, mp_rnd_t rnd_mode)
{
	if (a>0) return pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
	else	 return pow(mpreal(a),b,rnd_mode); //mpfr_pow_ui
}

inline const mpreal pow(const int a, const unsigned int b, mp_rnd_t rnd_mode)
{
	if (a>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode);  //mpfr_ui_pow_ui
	else	 return pow(mpreal(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const mpreal pow(const int a, const long int b, mp_rnd_t rnd_mode)
{
	if (a>0)
	{
		if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
		else	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
	}else{
		return pow(mpreal(a),b,rnd_mode); // mpfr_pow_si
	}
}

inline const mpreal pow(const int a, const int b, mp_rnd_t rnd_mode)
{
	if (a>0)
	{
		if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
		else	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
	}else{
		return pow(mpreal(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
	}
}

inline const mpreal pow(const int a, const long double b, mp_rnd_t rnd_mode)
{
	if (a>=0) 	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
	else		return pow(mpreal(a),mpreal(b),rnd_mode); //mpfr_pow
}

inline const mpreal pow(const int a, const double b, mp_rnd_t rnd_mode)
{
	if (a>=0) 	return pow(static_cast<unsigned long int>(a),mpreal(b),rnd_mode); //mpfr_ui_pow
	else		return pow(mpreal(a),mpreal(b),rnd_mode); //mpfr_pow
}

// pow long double 
inline const mpreal pow(const long double a, const long double b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),mpreal(b),rnd_mode);
}

inline const mpreal pow(const long double a, const unsigned long int b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),b,rnd_mode); //mpfr_pow_ui
}

inline const mpreal pow(const long double a, const unsigned int b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const mpreal pow(const long double a, const long int b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),b,rnd_mode); // mpfr_pow_si
}

inline const mpreal pow(const long double a, const int b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
}

inline const mpreal pow(const double a, const double b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),mpreal(b),rnd_mode);
}

inline const mpreal pow(const double a, const unsigned long int b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),b,rnd_mode); // mpfr_pow_ui
}

inline const mpreal pow(const double a, const unsigned int b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),static_cast<unsigned long int>(b),rnd_mode); // mpfr_pow_ui
}

inline const mpreal pow(const double a, const long int b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),b,rnd_mode); // mpfr_pow_si
}

inline const mpreal pow(const double a, const int b, mp_rnd_t rnd_mode)
{
	return pow(mpreal(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
}
}

// Explicit specialization of std::swap for mpreal numbers
// Thus standard algorithms will use efficient version of swap (due to Koenig lookup)
// Non-throwing swap C++ idiom: http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Non-throwing_swap
namespace std
{
	template <>
	inline void swap(mpfr::mpreal& x, mpfr::mpreal& y) 
	{ 
		return mpfr::swap(x, y); 
	}
}

#endif /* __MP_REAL_H__ */
