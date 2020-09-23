// Copyright (c) 1999,2003,2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Stefan Schirra, Sylvain Pion, Michael Hemmer


#ifndef CGAL_GMPZ_TYPE_H
#define CGAL_GMPZ_TYPE_H

#include <CGAL/basic.h>
#include <CGAL/IO/io.h>

#include <CGAL/gmp.h>
#include <mpfr.h>

#include <boost/operators.hpp>
#include <CGAL/Handle_for.h>

#include <string>
#include <locale>

namespace CGAL {

// TODO : benchmark without ref-counting, and maybe give the possibility
// to select ref-counting or not, then... => template class.

// Wrapper around mpz_t to get the destructor call mpz_clear.
struct Gmpz_rep
{
// FIXME : bug if ~() is called before an mpz_init*() is called.
// not a problem in practice, but not nice.
// maybe the mpz_init_set* functions should move back to Gmpz_rep.
// But then we should use the Storage_traits::construct/get...

  mpz_t mpZ;

  Gmpz_rep() {}
  ~Gmpz_rep() { mpz_clear(mpZ); }

private:
  // Make sure it does not get accidentally copied.
  Gmpz_rep(const Gmpz_rep &);
  Gmpz_rep & operator= (const Gmpz_rep &);
};


class Gmpz
  : Handle_for<Gmpz_rep>,
    boost::ordered_euclidian_ring_operators1< Gmpz
  , boost::ordered_euclidian_ring_operators2< Gmpz, int
  , boost::ordered_euclidian_ring_operators2< Gmpz, long
  , boost::ordered_euclidian_ring_operators2< Gmpz, unsigned long
  , boost::shiftable< Gmpz , long
  , boost::unit_steppable<Gmpz
  , boost::bitwise<Gmpz
> > > > > > >
{
  typedef Handle_for<Gmpz_rep> Base;
public:
  typedef Tag_true  Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_true  Has_exact_division;
  typedef Tag_false Has_exact_sqrt;

  Gmpz()
  { mpz_init(mpz()); }

  Gmpz(const mpz_t z)
  { mpz_init_set(mpz(), z); }

  Gmpz(int i)
  { mpz_init_set_si(mpz(), i); }

  Gmpz(long l)
  { mpz_init_set_si(mpz(), l); }

  Gmpz(unsigned long l)
  { mpz_init_set_ui(mpz(), l); }

  Gmpz(double d)
  {
     CGAL_warning_msg(is_integer(d), "Gmpz constructed from non-integer double value");
     CGAL_assertion(is_finite(d));
     mpz_init_set_d(mpz(), d);
   }

  Gmpz(const std::string& str, int base = 10)
  { mpz_init_set_str(mpz(), str.c_str(), base); }

  // returns the number of bits used to represent this number
  size_t bit_size() const { return mpz_sizeinbase(mpz(),2); }

  // returns the memory size in bytes
  size_t size() const { return mpz_size(mpz()) / (mp_bits_per_limb/8); }

  // returns the number of decimal digits needed to represent this number
  size_t approximate_decimal_length() const { return mpz_sizeinbase(mpz(),10); }

  double to_double() const {return mpz_get_d(mpz());}
  Sign sign() const { return static_cast<Sign>(mpz_sgn(mpz()));}

  const mpz_t & mpz() const { return Ptr()->mpZ; }
  mpz_t & mpz() { return ptr()->mpZ; }

  #ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
  int tam() const { return 0; }  // put here a code
                                 // measuring the number of digits
                                 // of the Gmpz
#endif

#define CGAL_GMPZ_OBJECT_OPERATOR(_op,_class,_fun)    \
  Gmpz& _op(const _class& z){                        \
    Gmpz Res;                                         \
    _fun(Res.mpz(), mpz(), z.mpz());                  \
    swap(Res);                                        \
    return *this;                                     \
  }

  CGAL_GMPZ_OBJECT_OPERATOR(operator+=,Gmpz,mpz_add);
  CGAL_GMPZ_OBJECT_OPERATOR(operator-=,Gmpz,mpz_sub);
  CGAL_GMPZ_OBJECT_OPERATOR(operator*=,Gmpz,mpz_mul);
  CGAL_GMPZ_OBJECT_OPERATOR(operator/=,Gmpz,mpz_tdiv_q);
  CGAL_GMPZ_OBJECT_OPERATOR(operator%=,Gmpz,mpz_tdiv_r);
  CGAL_GMPZ_OBJECT_OPERATOR(operator&=,Gmpz,mpz_and);
  CGAL_GMPZ_OBJECT_OPERATOR(operator|=,Gmpz,mpz_ior);
  CGAL_GMPZ_OBJECT_OPERATOR(operator^=,Gmpz,mpz_xor);
#undef CGAL_GMPZ_OBJECT_OPERATOR

  bool operator<(const Gmpz &b) const
  { return mpz_cmp(this->mpz(), b.mpz()) < 0; }
  bool operator==(const Gmpz &b) const
  { return mpz_cmp(this->mpz(), b.mpz()) == 0; }


  Gmpz operator+() const {return Gmpz( mpz() );}
  Gmpz operator-() const {
    Gmpz Res;
    mpz_neg(Res.mpz(), mpz());
    return Res;
  }

  Gmpz& operator <<= (const unsigned long& i){
    Gmpz Res;
    mpz_mul_2exp(Res.mpz(),this->mpz(), i);
    swap(Res);
    return *this;
  }
  Gmpz& operator >>= (const unsigned long& i){
    Gmpz Res;
    mpz_tdiv_q_2exp(Res.mpz(),this->mpz(), i);
    swap(Res);
    return *this;
  }

  Gmpz& operator++(){return *this+=1;}
  Gmpz& operator--(){return *this-=1;}


  // interoperability with int
  Gmpz& operator+=(int i);
  Gmpz& operator-=(int i);
  Gmpz& operator*=(int i);
  Gmpz& operator/=(int i);
  bool  operator==(int i) const {return mpz_cmp_si(this->mpz(), i) == 0;};
  bool  operator< (int i) const {return mpz_cmp_si(this->mpz(), i) < 0;};
  bool  operator> (int i) const {return mpz_cmp_si(this->mpz(), i) > 0;};

  // interoperability with long
  Gmpz& operator+=(long i);
  Gmpz& operator-=(long i);
  Gmpz& operator*=(long i);
  Gmpz& operator/=(long i);
  bool  operator==(long i) const {return mpz_cmp_si(this->mpz(), i) == 0;};
  bool  operator< (long i) const {return mpz_cmp_si(this->mpz(), i) < 0;};
  bool  operator> (long i) const {return mpz_cmp_si(this->mpz(), i) > 0;};

  // interoperability with unsigned long
  Gmpz& operator+=(unsigned long i);
  Gmpz& operator-=(unsigned long i);
  Gmpz& operator*=(unsigned long i);
  Gmpz& operator/=(unsigned long i);
  bool  operator==(unsigned long i) const {return mpz_cmp_ui(this->mpz(), i) == 0;};
  bool  operator< (unsigned long i) const {return mpz_cmp_ui(this->mpz(), i) < 0;};
  bool  operator> (unsigned long i) const {return mpz_cmp_ui(this->mpz(), i) > 0;};
};



#define CGAL_GMPZ_SCALAR_OPERATOR(_op,_type,_fun)   \
  inline Gmpz& Gmpz::_op(_type z) {                 \
    Gmpz Res;                                       \
    _fun(Res.mpz(), mpz(), z);                      \
    swap(Res);                                      \
    return *this;                                   \
  }

CGAL_GMPZ_SCALAR_OPERATOR(operator*=,int,mpz_mul_si)
CGAL_GMPZ_SCALAR_OPERATOR(operator*=,long,mpz_mul_si)

CGAL_GMPZ_SCALAR_OPERATOR(operator+=,unsigned long,mpz_add_ui)
CGAL_GMPZ_SCALAR_OPERATOR(operator-=,unsigned long,mpz_sub_ui)
CGAL_GMPZ_SCALAR_OPERATOR(operator*=,unsigned long,mpz_mul_ui)
CGAL_GMPZ_SCALAR_OPERATOR(operator/=,unsigned long,mpz_tdiv_q_ui)
#undef CGAL_GMPZ_SCALAR_OPERATOR


inline Gmpz& Gmpz::operator+=(int i)
{
  Gmpz Res;
  if (i >= 0)
    mpz_add_ui(Res.mpz(), mpz(), i);
  else
    mpz_sub_ui(Res.mpz(), mpz(), -i);
  swap(Res);
  return *this;
}

inline Gmpz& Gmpz::operator+=(long i)
{
  Gmpz Res;
  if (i >= 0)
    mpz_add_ui(Res.mpz(), mpz(), i);
  else
    mpz_sub_ui(Res.mpz(), mpz(), -i);
  swap(Res);
  return *this;
}



inline Gmpz& Gmpz::operator-=(int  i){return *this+=-i;}
inline Gmpz& Gmpz::operator-=(long i){return *this+=-i;}

inline Gmpz& Gmpz::operator/=(int b) {
  if (b>0) {
    Gmpz Res;
    mpz_tdiv_q_ui(Res.mpz(), mpz(), b);
    swap(Res);
    return *this;
  }
  return *this /= Gmpz(b);
}

inline Gmpz& Gmpz::operator/=(long b) {
  if (b>0) {
    Gmpz Res;
    mpz_tdiv_q_ui(Res.mpz(), mpz(), b);
    swap(Res);
    return *this;
  }
  return *this /= Gmpz(b);
}

inline
std::ostream&
operator<<(std::ostream& os, const Gmpz &z)
{
  char *str = new char [mpz_sizeinbase(z.mpz(),10) + 2];
  str = mpz_get_str(str, 10, z.mpz());
  os << str ;
  delete[] str;
  return os;
}


inline
std::istream &
gmpz_new_read(std::istream &is, Gmpz &z)
{
  bool negative = false;
  const std::istream::char_type zero = '0';
  std::istream::int_type c;
  Gmpz r;
  std::ios::fmtflags old_flags = is.flags();

  is.unsetf(std::ios::skipws);
  internal::eat_white_space(is);

  c=is.peek();
  if (c=='-' || c=='+'){
      is.get();
      CGAL_assertion(!is.fail());
      negative=(c=='-');
      internal::eat_white_space(is);
      c=is.peek();
  }

  std::istream::char_type cc= static_cast<std::istream::char_type>(c);

  if (c== std::istream::traits_type::eof() ||
      !std::isdigit(cc, std::locale::classic() ) ){
    is.setstate(std::ios_base::failbit);
  } else {
    CGAL_assertion(cc==c);
    r= cc-zero;
    is.get();
    CGAL_assertion(!is.fail());

    // The following loop was supposed to be an infinite loop with:
    //   while (true)
    // where the break condition is that is.peek() returns and EOF or a
    // non-digit character.
    //
    // Unfortunately, the wording of the C++03 and C++11 standard was not
    // well understood by the authors of libc++ (the STL of clang++) and,
    // in the version of libc++ shipped with Apple-clang-3.2,
    // istream::peek() set the flag eofbit when it reads the last character
    // of the stream *instead* of setting it only when it *tries to read
    // past the last character*. For that reason, to avoid that the next
    // peek() sets also the failbit, one has to check for EOL twice.
    //
    // See the LWG C++ Issue 2036, classified as Not-A-Defect:
    //   http://lwg.github.com/issues/lwg-closed.html#2036
    // and a StackOverflow related question:
    //   http://stackoverflow.com/a/9020292/1728537
    // --
    // Laurent Rineau, 2013/10/10
    while (!is.eof()) {
      c=is.peek();
      if (c== std::istream::traits_type::eof()) {
        break;
      }
      cc=static_cast<std::istream::char_type>(c);
      if  ( !std::isdigit(cc, std::locale::classic() )) {
        break;
      }
      is.get();
      CGAL_assertion(!is.fail());
      CGAL_assertion(cc==c);
      r= r*10+(cc-zero);
    }
  }

  is.flags(old_flags);
  if (!is.fail()) {
    if (negative) {
      z=-r;
    } else {
      z=r;
    }
  }
  return is;
}

/*inline
std::istream&
read_gmpz(std::istream& is, Gmpz &z) {
  bool negative = false;
  bool good = false;
  const int null = '0';
  char c;
  Gmpz tmp;
  std::ios::fmtflags old_flags = is.flags();

  is.unsetf(std::ios::skipws);
  while (is.get(c) && std::isspace(c, std::locale::classic() ))
  {}

  if (c == '-')
  {
        negative = true;
        while (is.get(c) && std::isspace(c, std::locale::classic() ))
        {}
  }
  if (std::isdigit(c, std::locale::classic() ))
  {
        good = true;
        tmp = c - null;
        while (is.get(c) && std::isdigit(c, std::locale::classic() ))
        {
            tmp = 10*tmp + (c-null);
        }
  }
  if (is)
        is.putback(c);
  if (sign(tmp) != ZERO && negative)
      tmp = -tmp;
  if (good){
      z = tmp;
      }
   else
    is.clear(is.rdstate() | std::ios::failbit);

  is.flags(old_flags);
  return is;
  }*/

inline
std::istream&
operator>>(std::istream& is, Gmpz &z)
{
  return gmpz_new_read(is, z);
}

template <>
struct Split_double<Gmpz>
{
  void operator()(double d, Gmpz &num, Gmpz &den) const
  {
    std::pair<double, double> p = split_numerator_denominator(d);
    num = Gmpz(p.first);
    den = Gmpz(p.second);
  }
};

inline Gmpz min BOOST_PREVENT_MACRO_SUBSTITUTION(const Gmpz& x,const Gmpz& y){
  return (x<=y)?x:y;
}
inline Gmpz max BOOST_PREVENT_MACRO_SUBSTITUTION(const Gmpz& x,const Gmpz& y){
  return (x>=y)?x:y;
}


} //namespace CGAL

#endif // CGAL_GMPZ_TYPE_H
