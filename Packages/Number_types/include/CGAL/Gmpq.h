// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/Gmpq.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_GMPQ_H
#define CGAL_GMPQ_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Gmpz.h>
#include <cassert>

#include <gmp.h>


CGAL_BEGIN_NAMESPACE

class Gmpq_rep
{
public:

  mpq_t  mpQ;

  Gmpq_rep()
  { 
    mpq_init(mpQ); 
  }

  Gmpq_rep(const mpq_t z)
  { 
    mpq_init(mpQ); 
    mpq_set(mpQ, z);
  }

  Gmpq_rep(const Gmpq_rep & g)
  { 
    mpq_init(mpQ); 
    mpq_set(mpQ, g.mpQ);
  }

  Gmpq_rep & operator= (const Gmpq_rep & g)
  {
      if (&g != this) {
	  mpq_clear(mpQ);
	  mpq_set(mpQ, g.mpQ);
      }
      return *this;
  }

  Gmpq_rep(int si)
  { 
    mpq_init(mpQ); 
    mpq_set_si(mpQ, si, 1); 
  }

  Gmpq_rep(long si)
  { 
    mpq_init(mpQ); 
    mpq_set_si(mpQ, si, 1); 
  }


  Gmpq_rep(unsigned long ui)
  { 
    mpq_init(mpQ); 
    mpq_set_ui(mpQ, ui, 1); 
  }

  Gmpq_rep(const Gmpz& z)
  { 
    mpq_init(mpQ); 
    mpq_set_z(mpQ, z.mpz()); 
  }

  Gmpq_rep(unsigned long int ui1, unsigned long int ui2)
  { 
    mpq_init(mpQ); 
    mpq_set_ui(mpQ, ui1, ui2); 
    mpq_canonicalize(mpQ);
  }
  
  Gmpq_rep(signed long int si, unsigned long int ui)
  { 
    mpq_init(mpQ); 
    mpq_set_si(mpQ, si, ui);
    mpq_canonicalize(mpQ);
  }
  
  Gmpq_rep(const Gmpz& n, const Gmpz& d)
  { 
    mpq_init(mpQ); 
    mpz_set(mpq_numref(mpQ), n.mpz());
    mpz_set(mpq_denref(mpQ), d.mpz());
    
    mpq_canonicalize(mpQ);
  }
  

  Gmpq_rep(double d)
  { 
    mpq_init(mpQ);
    mpq_set_d(mpQ, d); 
}

  Gmpq_rep(const char * const str)
  { 
    mpq_init(mpQ);
    mpq_set_str(mpQ, str, 10);
    mpq_canonicalize(mpQ); 
  }

  Gmpq_rep(const char * const str, int base)
  { 
    mpq_init(mpQ);
    mpq_set_str(mpQ, str, base);
    mpq_canonicalize(mpQ);
  }

  ~Gmpq_rep()
  { mpq_clear(mpQ); }
};

//__________________________________________________
class Gmpq
  : public Handle_for<Gmpq_rep>
{
  typedef Handle_for<Gmpq_rep> Base;
public:
  typedef Tag_false  Has_gcd;
  typedef Tag_true Has_division;
  typedef Tag_false  Has_sqrt;


  Gmpq()
    : Base(Gmpq_rep()) {}

  Gmpq(const mpq_t z)
    : Base(Gmpq_rep(z)) {}


  Gmpq(int n)
    : Base(Gmpq_rep(n)) {}

  Gmpq(long n)
    : Base(Gmpq_rep(n)) {}

  Gmpq(unsigned long n)
    : Base(Gmpq_rep(n)) {}

  Gmpq(const Gmpz& n)
    : Base(Gmpq_rep(n)) {}

  Gmpq(signed long n, unsigned long d)
    : Base(Gmpq_rep(n, d)) {}

  Gmpq(unsigned long n, unsigned long d)
    : Base(Gmpq_rep(n, d)) {}

  Gmpq(const Gmpz& n, const Gmpz& d)
    : Base(Gmpq_rep(n,d)) {}

  Gmpq(double d)
    : Base(Gmpq_rep(d)) {}

  Gmpq(const char* const str)
    : Base(Gmpq_rep(str)) {}

  Gmpq(const char* const str, int base)
    : Base(Gmpq_rep(str, base)) {}


  Gmpz numerator() const
  {
    return Gmpz(mpq_numref(mpq()));
  }

  Gmpz denominator() const
  {
    return Gmpz(mpq_denref(mpq()));
    
  }
  
  bool operator==(const Gmpq &z) const;

  bool operator!=(const Gmpq &z) const;

  bool operator<(const Gmpq &z) const;

  bool operator<=(const Gmpq &z) const;

  bool operator>(const Gmpq &z) const;

  bool operator>=(const Gmpq &z) const;

  Gmpq operator-() const;

  Gmpq operator+(const Gmpq &z) const;

  Gmpq operator-(const Gmpq &z) const;

  Gmpq operator*(const Gmpq &z) const;

  Gmpq operator/(const Gmpq &z) const;

  Gmpq& operator+=(const Gmpq &z);

  Gmpq& operator-=(const Gmpq &z);

  Gmpq& operator*=(const Gmpq &z);

  Gmpq& operator/=(const Gmpq &z);


  double to_double() const;
  Sign sign() const;

  const mpq_t & mpq() const { return Ptr()->mpQ; }
  mpq_t & mpq() { return ptr()->mpQ; }
};


inline
bool
Gmpq::operator==(const Gmpq &z) const
{ return mpq_equal(mpq(), z.mpq()); }

inline
bool
Gmpq::operator<(const Gmpq &z) const
{ return mpq_cmp(mpq(), z.mpq()) < 0; }



inline
bool
Gmpq::operator<=(const Gmpq &z) const
{ return mpq_cmp(mpq(), z.mpq()) <= 0; }


inline
bool
Gmpq::operator>(const Gmpq &z) const
{ return mpq_cmp(mpq(), z.mpq()) > 0; }


inline
bool
Gmpq::operator>=(const Gmpq &z) const
{ return mpq_cmp(mpq(), z.mpq()) >= 0; }


inline
bool
Gmpq::operator!=(const Gmpq &z) const
{ return ! (*this == z); }


inline
Gmpq
Gmpq::operator-() const
{
    Gmpq Res;
    mpq_neg(Res.mpq(), mpq());
    return Res;
}

inline
Gmpq
Gmpq::operator+(const Gmpq &z) const
{
    Gmpq Res;
    mpq_add(Res.mpq(), mpq(), z.mpq());
    return Res;
}


inline
Gmpq&
Gmpq::operator+=(const Gmpq &z)
{
    *this = *this + z;
    return *this;
}


inline
Gmpq
Gmpq::operator-(const Gmpq &z) const
{
    Gmpq Res;
    mpq_sub(Res.mpq(), mpq(), z.mpq());
    return Res;
}


inline
Gmpq&
Gmpq::operator-=(const Gmpq &z)
{
    *this = *this - z;
    return *this;
}

inline
Gmpq
Gmpq::operator*(const Gmpq &z) const
{
    Gmpq Res;
    mpq_mul(Res.mpq(), mpq(), z.mpq());
    return Res;
}

inline
Gmpq&
Gmpq::operator*=(const Gmpq &z)
{
    *this = *this * z;
    return *this;
}


inline
Gmpq
Gmpq::operator/(const Gmpq &z) const
{
    Gmpq Res;
    mpq_div(Res.mpq(), mpq(), z.mpq());
    return Res;
}




inline
Gmpq&
Gmpq::operator/=(const Gmpq &z)
{
    *this = *this / z;
    return *this;
}


inline
double
Gmpq::to_double() const
{ return mpq_get_d(mpq()); }

inline
io_Operator
io_tag(const Gmpq&)
{ return io_Operator(); }

inline
Sign
Gmpq::sign() const
{ return static_cast<Sign>(mpq_sgn(mpq())); }


inline
double
to_double(const Gmpq &z)
{ return z.to_double(); }

inline
Sign
sign(const Gmpq &z)
{ return z.sign(); }

inline
bool
is_valid(const Gmpq &)
{ return true; }

inline
bool
is_finite(const Gmpq &)
{ return true; }



inline
std::ostream&
operator<<(std::ostream& os, const Gmpq &z)
{
  os << z.numerator() << "/" << z.denominator();
  return os;
}

inline
std::istream&
operator>>(std::istream& is, Gmpq &z)
{
  char c;
  Gmpz n, d;
  is >> n;
  is >> c;
  assert(c == '/');
  is >> d;
  z = Gmpq(n,d);
  
  return is;
}

inline
Interval_base
to_interval (const Gmpq& z)
{
  return Interval_nt<>(CGAL::to_interval(z.numerator())) /
         Interval_nt<>(CGAL::to_interval(z.denominator()));
}




template <>
struct Rational_traits<Gmpq> {
  typedef Gmpz RT;
 RT  numerator   (const Gmpq & r) const { return r.numerator(); }
 RT  denominator (const Gmpq & r) const { return r.denominator(); }
 Gmpq make_rational(const RT & n, const RT & d) const
 { return Gmpq(n, d); } 
};

CGAL_END_NAMESPACE

#endif // CGAL_GMPQ_H
