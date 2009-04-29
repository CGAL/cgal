// Copyright (c) 1999,2003,2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Andreas Fabri, Stefan Schirra, Sylvain Pion


#ifndef CGAL_GMPZ_TYPE_H
#define CGAL_GMPZ_TYPE_H

#include <CGAL/basic.h>
#include <gmp.h>
#include <mpfr.h>

#include <boost/operators.hpp>
#include <CGAL/Handle_for.h>

#include <string>
#include <locale>

CGAL_BEGIN_NAMESPACE

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
    > >
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

  Gmpz operator+() const;
  Gmpz operator-() const;

  Gmpz& operator+=(const Gmpz &z);
  Gmpz& operator+=(int i);

  Gmpz& operator-=(const Gmpz &z);
  Gmpz& operator-=(int i);

  Gmpz& operator*=(const Gmpz &z);
  Gmpz& operator*=(int i);

  Gmpz& operator%=(const Gmpz &z);

  Gmpz& operator/=(const Gmpz &z);
  Gmpz& operator/=(int i);

  size_t approximate_decimal_length() const;

  double to_double() const;
  Sign sign() const;

  const mpz_t & mpz() const { return Ptr()->mpZ; }
  mpz_t & mpz() { return ptr()->mpZ; }

  #ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
  int tam() const { return 0; }  // put here a code
                                 // measuring the number of digits
                                 // of the Gmpz
#endif

  // Gives the memory size in bytes. (not documented yet)
  std::size_t size() const
  {
    return mpz_size(mpz()) / (mp_bits_per_limb/8);
  }
};


inline
bool
operator<(const Gmpz &a, const Gmpz &b)
{ return mpz_cmp(a.mpz(), b.mpz()) < 0; }

inline
bool
operator==(const Gmpz &a, const Gmpz &b)
{ return mpz_cmp(a.mpz(), b.mpz()) == 0; }


// mixed operators.
inline
bool
operator<(const Gmpz &a, int b)
{ return mpz_cmp_si(a.mpz(), b) < 0; }

inline
bool
operator==(const Gmpz &a, int b)
{ return mpz_cmp_si(a.mpz(), b) == 0; }

inline
bool
operator>(const Gmpz &a, int b)
{ return mpz_cmp_si(a.mpz(), b) > 0; }


inline
Gmpz
Gmpz::operator+() const {
  return Gmpz( mpz() );
}

inline
Gmpz
Gmpz::operator-() const
{
    Gmpz Res;
    mpz_neg(Res.mpz(), mpz());
    return Res;
}


inline
Gmpz&
Gmpz::operator+=(const Gmpz &z)
{
    Gmpz Res;
    mpz_add(Res.mpz(), mpz(), z.mpz());
    swap(Res);
    return *this;
}

inline
Gmpz&
Gmpz::operator+=(int i)
{
    Gmpz Res;
    if (i >= 0)
        mpz_add_ui(Res.mpz(), mpz(), i);
    else
        mpz_sub_ui(Res.mpz(), mpz(), -i);
    swap(Res);
    return *this;
}

inline
Gmpz&
Gmpz::operator-=(const Gmpz &z)
{
    Gmpz Res;
    mpz_sub(Res.mpz(), mpz(), z.mpz());
    swap(Res);
    return *this;
}

inline
Gmpz&
Gmpz::operator-=(int i)
{
    Gmpz Res;
    if (i >= 0)
        mpz_sub_ui(Res.mpz(), mpz(), i);
    else
        mpz_add_ui(Res.mpz(), mpz(), -i);
    swap(Res);
    return *this;
}

inline
Gmpz&
Gmpz::operator*=(const Gmpz &z)
{
    Gmpz Res;
    mpz_mul(Res.mpz(), mpz(), z.mpz());
    swap(Res);
    return *this;
}

inline
Gmpz&
Gmpz::operator*=(int i)
{
    Gmpz Res;
    mpz_mul_si(Res.mpz(), mpz(), i);
    swap(Res);
    return *this;
}

inline
Gmpz&
Gmpz::operator/=(const Gmpz &z)
{
    CGAL_precondition(z != 0);
    Gmpz Res;
    mpz_tdiv_q(Res.mpz(), mpz(), z.mpz());
    swap(Res);
    return *this;
}

inline
Gmpz&
Gmpz::operator/=(int b)
{
    if (b>0)
    {
        Gmpz Res;
        mpz_tdiv_q_ui(Res.mpz(), mpz(), b);
        swap(Res);
        return *this;
    }
    return *this /= Gmpz(b);
}

inline
Gmpz&
Gmpz::operator%=(const Gmpz &z)
{
    Gmpz Res;
    mpz_tdiv_r(Res.mpz(), mpz(), z.mpz());
    swap(Res);
    return *this;
}


inline
double
Gmpz::to_double() const
{ return mpz_get_d(mpz()); }

inline
Sign
Gmpz::sign() const
{ return static_cast<Sign>(mpz_sgn(mpz())); }

inline
size_t
Gmpz::approximate_decimal_length() const
{ return mpz_sizeinbase(mpz(),10); }

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
void gmpz_eat_white_space(std::istream &is)
{
  std::istream::int_type c;
  do {
    c= is.peek();
    if (c== std::istream::traits_type::eof())
      return;
    else {
      std::istream::char_type cc= c;
      if ( std::isspace(cc, std::locale::classic()) ) {
        is.get();
        // since peek succeeded, this should too
        CGAL_assertion(!is.fail());
      } else {
        return;
      }
    }
  } while (true);
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
  gmpz_eat_white_space(is);

  c=is.peek();
  if (c=='-' || c=='+'){
      is.get();
      CGAL_assertion(!is.fail());
      negative=(c=='-');
      gmpz_eat_white_space(is);
      c=is.peek();
  }

  std::istream::char_type cc= c;

  if (c== std::istream::traits_type::eof() ||
      !std::isdigit(cc, std::locale::classic() ) ){
    is.setstate(std::ios_base::failbit);
  } else {
    CGAL_assertion(cc==c);
    r= cc-zero;
    is.get();
    CGAL_assertion(!is.fail());
    while (true) {
      c=is.peek();
      if (c== std::istream::traits_type::eof()) {
	break;
      }
      cc=c;
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


#include <CGAL/auto_link/GMP.h>
#include <CGAL/auto_link/MPFR.h>


CGAL_END_NAMESPACE

#endif // CGAL_GMPZ_TYPE_H
