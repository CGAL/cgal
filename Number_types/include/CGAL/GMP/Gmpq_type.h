// Copyright (c) 2002,2003  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Andreas Fabri, Sylvain Pion


#ifndef CGAL_GMPQ_TYPE_H
#define CGAL_GMPQ_TYPE_H

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/IO/io.h>

#include <CGAL/GMP/Gmpz_type.h>
#include <CGAL/GMP/Gmpfr_type.h>

#include <CGAL/gmp.h>
#include <mpfr.h>
#include <utility>
#include <string>

#include <boost/operators.hpp>
#include <CGAL/Handle_for.h>
#include <CGAL/Profile_counter.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4146)
     // warning on - applied on unsigned number
#endif

namespace CGAL {

// Wrapper around mpq_t to get the destructor call mpq_clear.
// Contrary to mpz_t, there are no mpq_init_set_* functions,
// so we simply call mpq_init() here.
struct Gmpq_rep
{
  mpq_t mpQ;

  Gmpq_rep()  { mpq_init(mpQ); }
  ~Gmpq_rep() { mpq_clear(mpQ); }

private:
  // Make sure it does not get accidentally copied.
  Gmpq_rep(const Gmpq_rep &);
  Gmpq_rep & operator= (const Gmpq_rep &);
};


class Gmpq
  : Handle_for<Gmpq_rep>,
    boost::totally_ordered1< Gmpq
  , boost::ordered_field_operators2< Gmpq, int
  , boost::ordered_field_operators2< Gmpq, long
  , boost::ordered_field_operators2< Gmpq, long long
  , boost::ordered_field_operators2< Gmpq, double
  , boost::ordered_field_operators2< Gmpq, Gmpz
  , boost::ordered_field_operators2< Gmpq, Gmpfr
    > > > > > > >
{
  typedef Handle_for<Gmpq_rep> Base;
public:
  typedef Tag_false  Has_gcd;
  typedef Tag_true   Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_true   Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  Gmpq() {}

  Gmpq(const mpq_t q)
  { mpq_set(mpq(), q); }

  Gmpq(int n)
  { mpq_set_si(mpq(), n, 1); }

  Gmpq(unsigned int n)
  { mpq_set_ui(mpq(), n, 1); }

  Gmpq(long n)
  { mpq_set_si(mpq(), n, 1); }

  Gmpq(unsigned long n)
  { mpq_set_ui(mpq(), n, 1); }

private:
  void init_ull(unsigned long long n){
      CGAL_assertion(sizeof(long)==4 && sizeof(long long)==8);
      mpq_set_ui(mpq(), (unsigned long)(n>>32), 1);
      mpz_ptr z = mpq_numref(mpq());
      mpz_mul_2exp (z, z, 32);
      mpz_add_ui (z, z, (unsigned long)n);
  }
public:
  Gmpq(unsigned long long n)
  {
    if (n <= std::numeric_limits<unsigned long>::max BOOST_PREVENT_MACRO_SUBSTITUTION ())
      mpq_set_ui(mpq(), (unsigned long)n, 1);
    else
      init_ull(n);
  }

  Gmpq(long long n)
  {
    if (sizeof(long)==sizeof(long long))
      mpq_set_si(mpq(), (long)n, 1);
    else if (n>=0)
      init_ull(n);
    else {
      init_ull(-(unsigned long long)n);
      mpq_neg(mpq(), mpq());
    }
  }

  Gmpq(const Gmpz& n)
  { mpq_set_z(mpq(), n.mpz()); }

  Gmpq(int n, int d)
  {
    if (d < 0) {
      n = -n;
      d = -d;
    }
    mpq_set_si(mpq(), n, d);
    mpq_canonicalize(mpq());
  }

  Gmpq(signed long n, unsigned long d)
  {
    mpq_set_si(mpq(), n, d);
    mpq_canonicalize(mpq());
  }

  Gmpq(unsigned long n, unsigned long d)
  {
    mpq_set_ui(mpq(), n, d);
    mpq_canonicalize(mpq());
  }

  Gmpq(const Gmpz& n, const Gmpz& d)
  {
    mpz_set(mpq_numref(mpq()), n.mpz());
    mpz_set(mpq_denref(mpq()), d.mpz());
    mpq_canonicalize(mpq());
  }

  Gmpq(double d)
  {
    CGAL_assertion(is_finite(d));
    mpq_set_d(mpq(), d);
  }

  Gmpq(const Gmpfr &f)
  {
    std::pair<Gmpz,long> intexp=f.to_integer_exp();
    if(intexp.second<0){
            mpz_set(mpq_numref(mpq()),intexp.first.mpz());
            mpz_ui_pow_ui(mpq_denref(mpq()),2,-intexp.second);
    }else{
            mpz_mul_2exp(mpq_numref(mpq()),
                         intexp.first.mpz(),
                         (unsigned long)intexp.second);
            mpz_set_ui(mpq_denref(mpq()),1);
    }
    // mpq_canonicalize is needed only when the numerator is odd and not zero
    if(mpz_tstbit(intexp.first.mpz(),0)==0 && mpz_sgn(intexp.first.mpz())!=0)
        mpq_canonicalize(mpq());
    CGAL_assertion_msg(mpfr_cmp_q(f.fr(),mpq())==0,
                       "error in conversion Gmpfr->Gmpq");
  }

  Gmpq(const std::string& str, int base = 10)
  {
    mpq_set_str(mpq(), str.c_str(), base);
    mpq_canonicalize(mpq());
  }

#ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
  int tam() const { return 0; }  // put here a code
                                 // measuring the number of digits
                                 // of the Gmpq
// a possible code is:
//  int tam() const { return std::max(numerator().tam(),
//                                      denominator().tam()); }
// the same as Quotient<MP_Float>
#endif

  // Gives the memory size in bytes. (not documented yet)
  std::size_t size() const
  {
    std::size_t s_num = mpz_size(mpq_numref(mpq())) * (mp_bits_per_limb/8);
    std::size_t s_den = mpz_size(mpq_denref(mpq())) * (mp_bits_per_limb/8);
    return s_num + s_den;
  }

  Gmpz numerator() const
  { return Gmpz(mpq_numref(mpq())); }

  Gmpz denominator() const
  { return Gmpz(mpq_denref(mpq())); }

  Gmpq operator+() const;
  Gmpq operator-() const;

  Gmpq& operator+=(const Gmpq &q);
  Gmpq& operator-=(const Gmpq &q);
  Gmpq& operator*=(const Gmpq &q);
  Gmpq& operator/=(const Gmpq &q);

  bool operator==(const Gmpq &q) const { return mpq_equal(this->mpq(), q.mpq()) != 0;}
  bool operator< (const Gmpq &q) const { return mpq_cmp(this->mpq(), q.mpq()) < 0; }

  double to_double() const;
  Sign sign() const;

  const mpq_t & mpq() const { return Ptr()->mpQ; }
  mpq_t & mpq() { return ptr()->mpQ; }

  ~Gmpq()
  {
     CGAL_HISTOGRAM_PROFILER("[Gmpq sizes in log2 scale]",
                             (unsigned) ( ::log(double(size())) / ::log(double(2)) )  );
  }

  // Interoperability with int
  Gmpq& operator+=(int z){return (*this)+= Gmpq(z);}
  Gmpq& operator-=(int z){return (*this)-= Gmpq(z);}
  Gmpq& operator*=(int z){return (*this)*= Gmpq(z);}
  Gmpq& operator/=(int z){return (*this)/= Gmpq(z);}
  bool  operator==(int z) const {return mpq_cmp_si(mpq(),z,1)==0;}
  bool  operator< (int z) const {return mpq_cmp_si(mpq(),z,1)<0;}
  bool  operator> (int z) const {return mpq_cmp_si(mpq(),z,1)>0;}

  // Interoperability with long
  Gmpq& operator+=(long z){return (*this)+= Gmpq(z);}
  Gmpq& operator-=(long z){return (*this)-= Gmpq(z);}
  Gmpq& operator*=(long z){return (*this)*= Gmpq(z);}
  Gmpq& operator/=(long z){return (*this)/= Gmpq(z);}
  bool  operator==(long z) const {return mpq_cmp_si(mpq(),z,1)==0;}
  bool  operator< (long z) const {return mpq_cmp_si(mpq(),z,1)<0;}
  bool  operator> (long z) const {return mpq_cmp_si(mpq(),z,1)>0;}

  // Interoperability with long long
  Gmpq& operator+=(long long z){return (*this)+= Gmpq(z);}
  Gmpq& operator-=(long long z){return (*this)-= Gmpq(z);}
  Gmpq& operator*=(long long z){return (*this)*= Gmpq(z);}
  Gmpq& operator/=(long long z){return (*this)/= Gmpq(z);}
  bool  operator==(long long z) const {return (*this)== Gmpq(z);}
  bool  operator< (long long z) const {return (*this)<  Gmpq(z);}
  bool  operator> (long long z) const {return (*this)>  Gmpq(z);}

  // Interoperability with double
  Gmpq& operator+=(double d){return (*this)+= Gmpq(d);}
  Gmpq& operator-=(double d){return (*this)-= Gmpq(d);}
  Gmpq& operator*=(double d){return (*this)*= Gmpq(d);}
  Gmpq& operator/=(double d){return (*this)/= Gmpq(d);}
  bool  operator==(double d) const {return (*this)== Gmpq(d);}
  bool  operator< (double d) const {return (*this)<  Gmpq(d);}
  bool  operator> (double d) const {return (*this)>  Gmpq(d);}

  // Interoperability with Gmpz
  Gmpq& operator+=(const Gmpz&);
  Gmpq& operator-=(const Gmpz&);
  Gmpq& operator*=(const Gmpz&);
  Gmpq& operator/=(const Gmpz&);
  bool  operator==(const Gmpz &z) const {return (*this)== Gmpq(z);}
  bool  operator< (const Gmpz &z) const {return (*this)<  Gmpq(z);}
  bool  operator> (const Gmpz &z) const {return (*this)>  Gmpq(z);}

  // Interoperability with Gmpfr
  Gmpq& operator+=(const Gmpfr &f){return (*this)+= Gmpq(f);}
  Gmpq& operator-=(const Gmpfr &f){return (*this)-= Gmpq(f);}
  Gmpq& operator*=(const Gmpfr &f){return (*this)*= Gmpq(f);}
  Gmpq& operator/=(const Gmpfr &f){return (*this)/= Gmpq(f);}
  bool  operator==(const Gmpfr &f) const {return mpfr_cmp_q(f.fr(),mpq())==0;}
  bool  operator< (const Gmpfr &f) const {return mpfr_cmp_q(f.fr(),mpq())>0;}
  bool  operator> (const Gmpfr &f) const {return mpfr_cmp_q(f.fr(),mpq())<0;}
};


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
Gmpq::operator+() const
{
  return Gmpq(mpq());
}

inline
Gmpq
operator+(const Gmpq &x, const Gmpq &y)
{
    Gmpq Res;
    mpq_add(Res.mpq(), x.mpq(), y.mpq());
    return Res;
}

inline
Gmpq&
Gmpq::operator+=(const Gmpq &z)
{
    (*this + z).swap(*this);
    return *this;
}

inline
Gmpq
operator-(const Gmpq &x, const Gmpq &y)
{
    Gmpq Res;
    mpq_sub(Res.mpq(), x.mpq(), y.mpq());
    return Res;
}

inline
Gmpq&
Gmpq::operator-=(const Gmpq &z)
{
    (*this - z).swap(*this);
    return *this;
}

inline
Gmpq
operator*(const Gmpq &x, const Gmpq &y)
{
    Gmpq Res;
    mpq_mul(Res.mpq(), x.mpq(), y.mpq());
    return Res;
}

inline
Gmpq&
Gmpq::operator*=(const Gmpq &z)
{
    (*this * z).swap(*this);
    return *this;
}

inline
Gmpq
operator/(const Gmpq &x, const Gmpq &y)
{
    CGAL_precondition(y != 0);
    Gmpq Res;
    mpq_div(Res.mpq(), x.mpq(), y.mpq());
    return Res;
}

inline
Gmpq&
Gmpq::operator/=(const Gmpq &z)
{
    (*this / z).swap(*this);
    return *this;
}

inline
Gmpq& Gmpq::operator+=(const Gmpz &z){
  if(unique()){
    mpz_addmul(mpq_numref(mpq()),mpq_denref(mpq()),z.mpz());
  }else{
    Gmpq result;
    mpz_mul(mpq_numref(result.mpq()),
            mpq_denref(mpq()),
            z.mpz());
    mpz_add(mpq_numref(result.mpq()),
            mpq_numref(mpq()),
            mpq_numref(result.mpq()));
    mpz_set(mpq_denref(result.mpq()),mpq_denref(mpq()));
    swap(result);
  }
  return *this;
}

inline
Gmpq& Gmpq::operator-=(const Gmpz &z){
  if(unique()){
    mpz_submul(mpq_numref(mpq()),mpq_denref(mpq()),z.mpz());
  }else{
    Gmpq result;
    mpz_mul(mpq_numref(result.mpq()),
            mpq_denref(mpq()),
            z.mpz());
    mpz_sub(mpq_numref(result.mpq()),
            mpq_numref(mpq()),
            mpq_numref(result.mpq()));
    mpz_set(mpq_denref(result.mpq()),mpq_denref(mpq()));
    swap(result);
  }
  return *this;
}

inline
Gmpq& Gmpq::operator*=(const Gmpz &z){
  if(unique()){
    mpz_mul(mpq_numref(mpq()),mpq_numref(mpq()),z.mpz());
    mpq_canonicalize(mpq());
  }else{
    Gmpq result;
    mpz_mul(mpq_numref(result.mpq()),mpq_numref(mpq()),z.mpz());
    mpz_set(mpq_denref(result.mpq()),mpq_denref(mpq()));
    mpq_canonicalize(result.mpq());
    swap(result);
  }
  return *this;
}

inline
Gmpq& Gmpq::operator/=(const Gmpz &z){
  if(unique()){
    mpz_mul(mpq_denref(mpq()),mpq_denref(mpq()),z.mpz());
    mpq_canonicalize(mpq());
  }else{
    Gmpq result;
    mpz_mul(mpq_denref(result.mpq()),mpq_denref(mpq()),z.mpz());
    mpz_set(mpq_numref(result.mpq()),mpq_numref(mpq()));
    mpq_canonicalize(result.mpq());
    swap(result);
  }
  return *this;
}

inline
double
Gmpq::to_double() const
{ return mpq_get_d(mpq()); }

inline
Sign
Gmpq::sign() const
{ return static_cast<Sign>(mpq_sgn(mpq())); }

inline
std::ostream&
operator<<(std::ostream& os, const Gmpq &z)
{
  os << z.numerator() << "/" << z.denominator();
  return os;
}

// inline
// std::istream&
// operator>>(std::istream& is, Gmpq &z)
// {
//   char c;
//   Gmpz n, d;
//   is >> n;
//   is >> c;
//   //CGAL_assertion(!is || c == '/');
//   if (c != '/'){
//     is.setstate(std::ios_base::failbit);
//     return is;
//   }
//   is >> d;
//   if (!is.fail()) {
//     z = Gmpq(n,d);
//   }
//   return is;
// }



inline
std::istream&
operator>>(std::istream& is, Gmpq &z)
{
  internal::read_float_or_quotient<Gmpz,Gmpq>(is, z);
  return is;
}


inline Gmpq min BOOST_PREVENT_MACRO_SUBSTITUTION(const Gmpq& x,const Gmpq& y){
  return (x<=y)?x:y;
}
inline Gmpq max BOOST_PREVENT_MACRO_SUBSTITUTION(const Gmpq& x,const Gmpq& y){
  return (x>=y)?x:y;
}

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#include <CGAL/enable_warnings.h>

#endif // CGAL_GMPQ_TYPE_H
