/******************************************************************
 * Core Library, Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: BigInt.inl
 * $Id$ 
******************************************************************/
#ifdef CORE_INLINE

// some helper functions.

CORE_INLINE void lidia_error_handler(const char *f, const char *m)
{
  std::cout << "\n error_handler";
  std::cout << "::" << f;
  std::cout << "::" << m;
  std::cout << "\n";
  std::cout.flush();
  abort();
}

CORE_INLINE void memory_handler(char *t, const char *f, const char *m)
{
  if (t == NULL) {
    std::cout << "\n memory_handler";
    std::cout << "::" << f;
    std::cout << "::" << m;
    std::cout << "memory exhausted\n";
    std::cout.flush();
    abort();
  }
}

/**
** constructors and destructor; we could leave out some of these
**/

CORE_INLINE bigint::bigint()
{ mpz_init(I); } 

CORE_INLINE bigint::bigint(int i)
{ mpz_init_set_si(I, i); } 

CORE_INLINE bigint::bigint(unsigned int ui)
{ mpz_init_set_ui(I, ui); } 

CORE_INLINE bigint::bigint(long l)
{ mpz_init_set_si(I, l); }

CORE_INLINE bigint::bigint(unsigned long ul)
{ mpz_init_set_ui(I, ul); }

CORE_INLINE void bigint::assign(double d)
{ mpz_set_d(I, d); }

// new constructor
CORE_INLINE bigint::bigint(const char *str) 
{ mpz_init_set_str(I, str, 0); }
 
CORE_INLINE bigint::bigint(const mpz_t & a)
{ mpz_init_set(I, a); }

CORE_INLINE bigint::bigint(const bigint & a)
{ mpz_init_set(I, a.I); }

CORE_INLINE bigint::~bigint()
{ mpz_clear(I); }

/**
** CORE_INLINE member functions
**/

CORE_INLINE int bigint::bit(unsigned int i) const
{ return mpz_tstbit(I, i); }

CORE_INLINE int bigint::length() const
{ return mpz_size(I); } 

CORE_INLINE int bigint::bit_length() const
{ return mpz_sizeinbase(I, 2); } 

CORE_INLINE int bigint::sign() const
{ return mpz_sgn(I); } 

CORE_INLINE bool bigint::is_positive() const
{ return mpz_sgn(I) > 0; } 

CORE_INLINE bool bigint::is_negative() const
{ return mpz_sgn(I) < 0; }

CORE_INLINE bool bigint::is_even() const
{ return mpz_even_p(I); }

CORE_INLINE bool bigint::is_odd() const
{ return mpz_odd_p(I); }

CORE_INLINE bool bigint::is_zero() const
{ return mpz_sgn(I) == 0; }

CORE_INLINE bool bigint::is_gt_zero() const
{ return mpz_sgn(I) > 0; }

CORE_INLINE bool bigint::is_ge_zero() const
{ return mpz_sgn(I) >= 0; }

CORE_INLINE bool bigint::is_lt_zero() const
{ return mpz_sgn(I) < 0; }

CORE_INLINE bool bigint::is_le_zero() const
{ return mpz_sgn(I) <= 0; }

CORE_INLINE bool bigint::is_one() const
{ return mpz_cmp_ui(I, 1) == 0; }


CORE_INLINE bool is_positive (const bigint & a)
{ return mpz_sgn(a.I) > 0; }

CORE_INLINE bool is_negative (const bigint & a)
{ return mpz_sgn(a.I) < 0; }

CORE_INLINE bool is_even (const bigint & a)
{ return mpz_even_p(a.I); }

CORE_INLINE bool is_odd (const bigint & a)
{ return mpz_odd_p(a.I); }

CORE_INLINE bool is_zero (const bigint & a)
{ return mpz_sgn(a.I) == 0; }

CORE_INLINE bool is_one (const bigint & a)
{ return mpz_cmp_ui(a.I, 1) == 0; }



CORE_INLINE bool bigint::intify(int & i) const
{ 
  if (!mpz_fits_sint_p(I))
    return 1;
  i = (int) mpz_get_si(I);
  return 0;
}

CORE_INLINE bool bigint::longify(long & i) const
{ 
  if (!mpz_fits_slong_p(I))
    return 1;
  i = mpz_get_si(I);
  return 0;
}

CORE_INLINE int bigint::abs_compare(const bigint & a) const
{
  int cmp = mpz_cmpabs(I, a.I);
  if (cmp < 0) return -1;
  return (cmp > 0);
}

CORE_INLINE int bigint::compare(const bigint & a) const
{
  int cmp = mpz_cmp(I, a.I);
  if (cmp < 0) return -1;
  return (cmp > 0);
}

CORE_INLINE unsigned long bigint::most_significant_digit() const
{ 
  long l = mpz_size(I);
  if (l == 0)
    return 0;
  return I[0]._mp_d[l - 1];
}

CORE_INLINE unsigned long bigint::least_significant_digit() const
{
  if (mpz_sgn(I) == 0)
    return 0;
  return I[0]._mp_d[0];
}

CORE_INLINE const double bigint::radix()
{ return ldexp(1.0, mp_bits_per_limb); }

CORE_INLINE const int bigint::bits_per_digit()
{ return mp_bits_per_limb; }

CORE_INLINE void bigint::absolute_value()
{ mpz_abs(I, I); }

CORE_INLINE void bigint::abs()
{ mpz_abs(I, I); }

CORE_INLINE void bigint::negate()
{ mpz_neg(I, I); }

CORE_INLINE void bigint::assign_zero()
{  mpz_set_ui(I, 0); }

CORE_INLINE void bigint::assign_one()
{ mpz_set_ui(I, 1); }

CORE_INLINE void bigint::assign(int i)
{ mpz_set_si(I, i); }

CORE_INLINE void bigint::assign(long i)
{ mpz_set_si(I, i); }

CORE_INLINE void bigint::assign(unsigned long ui)
{ mpz_set_ui(I, ui); }

CORE_INLINE void bigint::assign(const bigint & a) 
{ mpz_set(I, a.I); }

CORE_INLINE void bigint::multiply_by_2() 
{ mpz_mul_2exp(I, I, 1); }

CORE_INLINE void bigint::divide_by_2()
{ mpz_tdiv_q_2exp(I, I, 1); }

/**
** Type checking
**/

CORE_INLINE bool is_char(const bigint & a)
{ return a.bit_length() <= (BITS_PER_CHAR - 1); }

CORE_INLINE bool is_uchar(const bigint & a)
{ return ( a.bit_length() <= (BITS_PER_CHAR) && !a.is_negative() ); }

CORE_INLINE bool is_short(const bigint & a)
{ return mpz_fits_sshort_p(a.I); }

CORE_INLINE bool is_ushort(const bigint & a)
{ return mpz_fits_ushort_p(a.I); }

CORE_INLINE bool is_int(const bigint & a)
{ return mpz_fits_sint_p(a.I); }

CORE_INLINE bool is_uint(const bigint & a)
{ return mpz_fits_uint_p(a.I); }

CORE_INLINE bool is_long(const bigint & a)
{ return mpz_fits_slong_p(a.I); }

CORE_INLINE bool is_ulong(const bigint & a)
{ return mpz_fits_ulong_p(a.I); }

/**
** assignments
**/


CORE_INLINE int bigint::operator = (int i)
{ mpz_set_si(I, i); return i; }

CORE_INLINE long bigint::operator = (long l)
{ mpz_set_si(I, l); return l; }

CORE_INLINE unsigned long bigint::operator = (unsigned long ul)
{ mpz_set_ui(I, ul); return ul; }

CORE_INLINE double bigint::operator = (double d)
{
  mpz_set_d(I, d);
  return d;
} 

CORE_INLINE bigint & bigint::operator = (const bigint & a)
{ 
  if (this == &a)       // Bug fixed: avoid self-assignment
    return *this;       // Zilin Du, 07/11/01
  mpz_set(I, a.I); return *this; }

/**
** comparisons
**/

CORE_INLINE bool operator == (const bigint & a, const bigint & b)
{ return mpz_cmp(a.I, b.I) == 0;}
 
CORE_INLINE bool operator != (const bigint & a, const bigint & b)
{ return mpz_cmp(a.I, b.I) != 0;}
 
CORE_INLINE bool operator > (const bigint & a, const bigint & b)
{ return mpz_cmp(a.I, b.I) > 0;}
 
CORE_INLINE bool operator >= (const bigint & a, const bigint & b)
{ return mpz_cmp(a.I, b.I) >= 0;}
 
CORE_INLINE bool operator < (const bigint & a, const bigint & b)
{ return mpz_cmp(a.I, b.I) < 0;}
 
CORE_INLINE bool operator <= (const bigint & a, const bigint & b)
{ return mpz_cmp(a.I, b.I) <= 0;}

/**
** operator overloading
**/

CORE_INLINE bigint operator - (const bigint & a)
{ bigint c; mpz_neg(c.I, a.I); return c; }

CORE_INLINE bigint operator + (const bigint & a, const bigint & b)
{ bigint c; mpz_add(c.I, a.I, b.I); return c; }

CORE_INLINE bigint operator - (const bigint & a, const bigint & b)
{ bigint c; mpz_sub(c.I, a.I, b.I); return c; }

CORE_INLINE bigint operator * (const bigint & a, const bigint & b)
{ bigint c; mpz_mul(c.I, a.I, b.I); return c; }

CORE_INLINE bigint operator / (const bigint & a, const bigint & b)
{ bigint c; mpz_tdiv_q(c.I, a.I, b.I); return c; }

CORE_INLINE bigint operator % (const bigint & a, const bigint & b)
{ 
  bigint c;
  mpz_mod(c.I, a.I, b.I); // returns non-negative rest c

  if (a < 0 && c != 0)
   {
     // a is negative and rest != 0 -> create negative rest
     if (b < 0)
        mpz_add(c.I, c.I, b.I);
     else
        mpz_sub(c.I, c.I, b.I);
   }
  return c;
}

CORE_INLINE bigint operator << (const bigint & a, long ui)
{ 
  bigint c;
  if (ui < 0) {
#ifdef DEBUG
    lidia_error_handler("bigint", "operator<<::index is negative.");
#endif
    return (a >> (-ui));
  }
  mpz_mul_2exp(c.I, a.I, ui);
  return c; 
}

CORE_INLINE bigint operator >> (const bigint & a, long ui)
{ 
  bigint c; 
  if (ui < 0) {
#ifdef DEBUG
    lidia_error_handler("bigint", "operator>>::index is negative.");
#endif
    return (a << (-ui));
  }
  mpz_tdiv_q_2exp(c.I, a.I, ui); 
  return c; 
}

CORE_INLINE bigint operator & (const bigint & a, const bigint & b)
{ bigint c; mpz_and(c.I, a.I, b.I); return c; }

CORE_INLINE bigint operator | (const bigint & a, const bigint & b)
{ bigint c; mpz_ior(c.I, a.I, b.I); return c; }

CORE_INLINE bigint operator ^ (const bigint & a, const bigint & b)
{ bigint c = (a & ~b) | (~a & b); return c; }

CORE_INLINE bigint & bigint::operator += (const bigint & a)
{ mpz_add(I, I, a.I); return *this; }

CORE_INLINE bigint & bigint::operator -= (const bigint & a)
{ mpz_sub(I, I, a.I); return *this; }

CORE_INLINE bigint & bigint::operator *= (const bigint & a)
{ mpz_mul(I, I, a.I); return *this; }

CORE_INLINE bigint & bigint::operator /= (const bigint & a)
{ mpz_tdiv_q(I, I, a.I); return *this; }


CORE_INLINE bigint & bigint::operator %= (const bigint & a)
{
  if (&a == this)
    assign_zero();
  else
  {
    bool is_neg = (*this < 0);
    mpz_mod(I, I, a.I);      // *this = *this % |a|
    if (*this != 0)
    {
      if (a < 0)  // create negative rest
        mpz_add(I, I, a.I);
      else if (is_neg)
        mpz_sub(I, I, a.I);
    }
  }
  return *this;
}

CORE_INLINE bigint & bigint::operator <<= (long ui)
{ 
  if (ui < 0) 
    lidia_error_handler("bigint", "operator<<=::index is negative.");
  mpz_mul_2exp(I, I, ui); 
  return *this; 
}

CORE_INLINE bigint & bigint::operator >>= (long ui)
{ 
  if (ui < 0) 
    lidia_error_handler("bigint", "operator>>=::index is negative.");
  mpz_tdiv_q_2exp(I, I, ui);
  return *this;
}

CORE_INLINE bigint & bigint::operator &= (const bigint & a)
{ mpz_and(I, I, a.I); return *this; }

CORE_INLINE bigint & bigint::operator |= (const bigint & a)
{ mpz_ior(I, I, a.I); return *this; }

CORE_INLINE bigint & bigint::operator ^= (const bigint & a)
{ mpz_xor(I, I, a.I); return *this; }

CORE_INLINE bigint & bigint::operator++ ()
{ mpz_add_ui(I, I, 1); return *this; }

CORE_INLINE bigint & bigint::operator-- ()
{ mpz_sub_ui(I, I, 1); return *this; }

CORE_INLINE bigint  bigint::operator++ (int)
{ bigint a = *this; mpz_add_ui(I, I, 1); return a; }

CORE_INLINE bigint  bigint::operator-- (int)
{ bigint a =*this; mpz_sub_ui(I, I, 1); return a; }

CORE_INLINE int bigint::operator ! () const
{ return mpz_size(I) == 0; }

CORE_INLINE bigint bigint::operator ~ () const
{ bigint c; mpz_com(c.I, I); return c; }

/**
** Procedural versions
**/

//CORE_INLINE void negate(bigint & a, const bigint & b)
//{ mpz_neg(a.I, b.I); }

CORE_INLINE void add(bigint & c, const bigint & a, const bigint & b)
{ mpz_add(c.I, a.I, b.I); }

CORE_INLINE void subtract(bigint & c, const bigint & a, const bigint & b)
{ mpz_sub(c.I, a.I, b.I); }

CORE_INLINE void multiply(bigint & c, const bigint & a, const bigint & b)
{ mpz_mul(c.I, a.I, b.I); }

CORE_INLINE void divide(bigint & c, const bigint & a, const bigint & b)
{ mpz_tdiv_q(c.I, a.I, b.I); }

CORE_INLINE void div_rem(bigint & q, bigint & r, const bigint & a, const bigint & b)
{ mpz_tdiv_qr(q.I, r.I, a.I, b.I); }

CORE_INLINE void invert(bigint & a, const bigint & b)
{
  if (mpz_cmpabs_ui(b.I, 1) == 0) // b == 1 or b == -1.
    mpz_set(a.I, b.I);
  else
    lidia_error_handler("bigint", "invert::inverting of a non-unit.");
}

CORE_INLINE void shift_left(bigint & c, const bigint & a, long ui)
{ 
  if (ui < 0) 
    lidia_error_handler("bigint", "shift_left()::index is negative.");
  mpz_mul_2exp(c.I, a.I, ui);
}

CORE_INLINE void shift_right(bigint & c, const bigint & a, long ui)
{ 
  if (ui < 0) 
    lidia_error_handler("bigint", "shift_right()::index is negative.");
  mpz_tdiv_q_2exp(c.I, a.I, ui);
}

/*CORE_INLINE void and(bigint & c, const bigint & a, const bigint & b)
{ mpz_and(c.I, a.I, b.I); }

CORE_INLINE void or(bigint & c, const bigint & a, const bigint & b)
{ mpz_ior(c.I, a.I, b.I); }

CORE_INLINE void xor(bigint & c, const bigint & a, const bigint & b)
{ c = (~a & b) | (~b & a); }

CORE_INLINE void not(bigint & b, const bigint & a)
{ mpz_com(b.I, a.I); }
*/
CORE_INLINE void inc(bigint & c)
{ mpz_add_ui(c.I, c.I, 1); }

CORE_INLINE void dec(bigint & c)
{ mpz_sub_ui(c.I, c.I, 1); }

CORE_INLINE void add(bigint & c, const bigint & a, long i)
{
  if (i >= 0)
    mpz_add_ui(c.I, a.I, i);
  else
    mpz_sub_ui(c.I, a.I, -i);
}

CORE_INLINE void subtract(bigint & c, const bigint & a, long i)
{
  if (i >= 0)
    mpz_sub_ui(c.I, a.I, i);
  else
    mpz_add_ui(c.I, a.I, -i);
}

CORE_INLINE void multiply(bigint & c, const bigint & a, long i)
{
  if (i >= 0)
    mpz_mul_ui(c.I, a.I, i);
  else
  {
    mpz_mul_ui(c.I, a.I, -i);
    mpz_neg(c.I, c.I);
  }
}

CORE_INLINE void divide(bigint & q, const bigint & a, long i)
{
  if (i > 0)
    mpz_tdiv_q_ui(q.I, a.I, (unsigned long)i);
  else
  {
    mpz_tdiv_q_ui(q.I, a.I, (unsigned long)(-i));
    mpz_neg(q.I, q.I);
  }
}
  

CORE_INLINE void remainder(long &r, const bigint & a, long i)
{
  if (i > 0)
    r = (long)mpz_tdiv_ui(a.I, i);
  else
    r = (long)mpz_tdiv_ui(a.I, -i);
  
  if (r != 0 && mpz_sgn(a.I) < 0) // r must be negative
    r = -r;
}

CORE_INLINE void remainder(unsigned long &r, const bigint & a, unsigned long i)
{ r = mpz_tdiv_ui(a.I, i); }

CORE_INLINE long remainder(const bigint & a, long i)
{
  long r;
  if (i > 0)
    r = (long)mpz_tdiv_ui(a.I, i);
  else
    r = (long)mpz_tdiv_ui(a.I, -i);
  
  if (r != 0 && mpz_sgn(a.I) < 0) // r must be negative
    r = -r;
  return r;
}

CORE_INLINE void div_rem(bigint & q, long &r, const bigint & a, long i)
{
  if (i > 0)
    r = (long) mpz_tdiv_q_ui(q.I, a.I, (unsigned long)i);
  else
  {
    r = (long) mpz_tdiv_q_ui(q.I, a.I, (unsigned long)(-i));
    mpz_neg(q.I, q.I);
  }

  if (r != 0 && mpz_sgn(a.I) < 0) // r must be negative
    r = -r;
}

/**
** gcd's
**/

CORE_INLINE bigint gcd(const bigint & a, const bigint & b)
{ bigint c; mpz_gcd(c.I, a.I, b.I); return c; }

CORE_INLINE bigint bgcd(const bigint & a, const bigint & b)
{ bigint c; mpz_gcd(c.I, a.I, b.I); return c; }

CORE_INLINE bigint dgcd(const bigint & a, const bigint & b)
{ bigint c; mpz_gcd(c.I, a.I, b.I); return c; }

CORE_INLINE bigint xgcd(bigint & u, bigint & v, const bigint & a, const bigint & b)
{ bigint c; mpz_gcdext(c.I, u.I, v.I, a.I, b.I); return c; }

CORE_INLINE bigint xgcd_left(bigint & u, const bigint & a, const bigint & b)
{ bigint c; mpz_gcdext(c.I, u.I, 0, a.I, b.I); return c;}

CORE_INLINE bigint xgcd_right(bigint & v, const bigint & a, const bigint & b)
{ bigint c; mpz_gcdext(c.I, v.I, 0, b.I, a.I); return c; }

/**
** functions
**/

CORE_INLINE bigint abs(const bigint & a)
{ bigint c; mpz_abs(c.I, a.I); return c; }    

CORE_INLINE static gmp_randstate_t * get_randstate()
{
    static gmp_randstate_t rstate;
    static bool initialized = false;
    if (!initialized)
    {
        gmp_randinit(rstate, GMP_RAND_ALG_DEFAULT, 32L);
        initialized = true;
    }
    return &rstate;
}

CORE_INLINE void seed(const bigint & a)
{
    mpz_t tmp;
    mpz_init_set(tmp, a.I);
    gmp_randseed(*get_randstate(), tmp);
    mpz_clear(tmp);
}

CORE_INLINE void bigint::randomize(const bigint & a)
{
  if (a.is_zero())
     lidia_error_handler("bigint", "Bound must not be equal to zero.");
  else
   {
     bigint tmp;
     mpz_urandomb(tmp.I, *get_randstate(), mpz_sizeinbase(a.I, 2));
  
     if (a.is_lt_zero())
     {
         if (tmp <= a)
             remainder(tmp, tmp, a);
     }
     else if (tmp >= a)
         remainder(tmp, tmp, a);
     swap(*this, tmp);
   }
}

CORE_INLINE bigint randomize(const bigint & a)
{ 
  bigint c;
  c.randomize(a);
  return c;
}


#ifdef USE_OLD_DBL

CORE_INLINE double dbl(const bigint & a)
{
  char *s;
  int t = (int) mpz_sizeinbase(&a.I, 10) + 2;
  s = new char [t];
  mpz_get_str(s, 10, &a.I);
  double d;
  d = atof(s);
  delete[] s;
  return d;
}

/*xdouble xdbl(const bigint & a)
{
  char *s;
  int t = (int) mpz_sizeinbase(&a.I, 10) + 2;
  s = new char [t];
  mpz_get_str(s, 10, &a.I);
  xdouble d;
  d = string_to_xdouble(s);
  delete[] s;
  return d;
}*/

#else

CORE_INLINE double dbl(const bigint & a)
{ return mpz_get_d(a.I); }

/*xdouble xdbl(const bigint & a)
{
  long l = a.I._mp_size;
  l=((l<0)?-l:l);
  xdouble d = 0.0;
  if (l)
  {
    int i = 1;
    d = (double )a.I._mp_d[0]; 
    xdouble base = bigint::radix();
    xdouble power = base;
    while (i < l)
    {
     d += a.I._mp_d[i] * power;
     power *= base;
     i++;
    }
    if (a.I._mp_size < 0)
      d = -d;
  }
  return d;
}*/

#endif

CORE_INLINE void sqrt(bigint & a, const bigint & b)
{ 
   if (b.is_lt_zero())
      lidia_error_handler("bigint","sqrt(bigint&a,const bigint&b):: b < 0");
   else
      mpz_sqrt(a.I, b.I);
}
      
CORE_INLINE void square(bigint & a, const bigint & b)
{ mpz_mul(a.I, b.I, b.I); }

CORE_INLINE void swap(bigint & a, bigint & b)
{ mpz_swap(a.I, b.I); }



/**
 ** input / output
 **/

CORE_INLINE std::istream & operator >> (std::istream & in, bigint & a)
{     
  a.scan (in);
  return (in);
}

CORE_INLINE std::ostream & operator << (std::ostream & out, const bigint & a)
{
  int l = (int) mpz_sizeinbase(a.I, 10) + 10;
  char *s = new char [l];
  mpz_get_str(s, 10, a.I);

  a.print(out,s);

  delete[] s;
  return out;
}

CORE_INLINE int string_to_bigint(const char *s, bigint & a, int base)
{ 
  if (!mpz_set_str(a.I, s, base))
    return (strlen(s));
  else
    return 0;
}

CORE_INLINE int bigint_to_string(const bigint & a, char *s, int base)
{ mpz_get_str(s, base, a.I); return strlen(s); }

/**
 ** input / output
 **/


// number of characters in one line
// when printing a bigint

CORE_INLINE void bigint
::set_chars_per_line (int cpl)
 {
   bigint::chars_per_line = cpl;
 }

CORE_INLINE int bigint
::get_chars_per_line ()
 {
   return bigint::chars_per_line;
 }

//
// characteristic
//

CORE_INLINE bigint bigint
::characteristic() const
 { 
   return 0;
 }

/*****************************************************************
 * Core Additions to LiDIA's definitions above
 *****************************************************************/

// This function is dangerous: if bigint a does fit into
//      machine long, then tmp is set to zero.

CORE_INLINE long bigIntToLong(const bigint & a) {
  long tmp;
  if (a.longify(tmp))   // a.longify(tmp) returns 1 if
                        // a does not fit into a long, and tmp is unassigned.
        tmp = 0;        
                        // Otherwise a.longify(tmp) returns 0.
  return tmp;
}

// Probably dangerous too:
CORE_INLINE double bigIntToDouble(const bigint & a) {
  return dbl(a);
}

// lg(a) returns the bitlength of the argument a
//
// N.B. The name "lg" suggest that this is log to base 2.
//      But it is a misnomer!  The function "lg(a)" really computes
//      the bitlength of |a|, i.e., "ceiling(log_2(1 + |a| ))"
//      This function is deprecated.

CORE_INLINE int lg(const bigint & a){        // lg(a) = ceil(log_2(1+ |a|))
        return a.bit_length();          // Note: sign is ignored.
}

//  This is a more self-descriptive name for the previous "lg()"
//  USE THIS NAME in place of lg().

CORE_INLINE int bitLength(const bigint & a){
			 		// bitLength(a) = ceil(log_2(1+ |a|))
        return a.bit_length();      
}

CORE_INLINE int ceilLg(const bigint & a) {   // this is ceiling of log_2(a):
	// CONVENTION: log_2(0) = -1  (slight improvement, Zilin 8/5/02)
	if (a.is_zero()) return -1;
        unsigned int len = a.bit_length(); 
	return (mpz_scan0(a.I, 0) == len-1) ? (len-1) : len;
}

CORE_INLINE int floorLg(const bigint & a) {  // this is floor of log_2(a)
	// CONVENTION: log_2(0) = -1  (slight improvement, Zilin 8/5/02)
	if (a.is_zero()) return -1;  
        int len = a.bit_length(); 
        return (len - 1);
}
        
CORE_INLINE int compare(const bigint &a, const bigint &b){ return a.compare(b);}

CORE_INLINE int sign(const bigint &a) { return a.sign(); }

#endif

