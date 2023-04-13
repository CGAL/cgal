// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)        :  Andreas Fabri

#ifndef CGAL_CPP_FLOAT_H
#define CGAL_CPP_FLOAT_H

//#define CGAL_CPPF

#include <CGAL/boost_mp.h>
#include <CGAL/assertions.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>

namespace CGAL {


namespace internal {
  // Only used with an argument known not to be 0.
  inline int low_bit (boost::uint64_t x) {
#if defined(_MSC_VER)
    unsigned long ret;
    _BitScanForward64(&ret, x);
    return (int)ret;
#elif defined(__xlC__)
    return __cnttz8 (x);
#else
    // Assume long long is 64 bits
    return __builtin_ctzll (x);
#endif
  }
  inline int high_bit (boost::uint64_t x) {
#if defined(_MSC_VER)
    unsigned long ret;
    _BitScanReverse64(&ret, x);
    return (int)ret;
#elif defined(__xlC__)
    // Macro supposedly not defined on z/OS.
    return 63 - __cntlz8 (x);
#else
    return 63 - __builtin_clzll (x);
#endif
  }

} // namespace internal

#if 0 // needs C++20
  template <class T>
  void fmt(const T& t)
  {
    std::cout << std::format("{:b}", t) << std::endl;
  }
#endif

class cpp_float {
#ifdef CGAL_CPPF
  boost::multiprecision::cpp_rational rat;
#endif

  typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<512> > Mantissa;
  Mantissa man;
  int exp; /* The number man (an integer) * 2 ^ exp  */

public:
  cpp_float()
    :
#ifdef CGAL_CPPF
    rat(),
#endif
    man(), exp()
  {}

  cpp_float(short i)
    :
#ifdef CGAL_CPPF
    rat(i),
#endif
    man(i),exp(0)
  {}

  cpp_float(int i)
    :
#ifdef CGAL_CPPF
    rat(i),
#endif
    man(i),exp(0)
  {}

  cpp_float(long i)
    :
#ifdef CGAL_CPPF
    rat(i),
#endif
    man(i),exp(0)
  {}
#ifdef CGAL_CPPF
  cpp_float(const Mantissa& man,  int exp, const boost::multiprecision::cpp_rational& rat)
    : rat(rat), man(man),exp(exp)
  {}
#else

  cpp_float(const Mantissa& man, int exp)
      : man(man), exp(exp)
  {}

#ifndef CGAL_CPPF
  template <typename Expression>
  cpp_float(const Expression& man, int exp)
    :man(man), exp(exp)
  {}
#endif

#endif
  cpp_float(double d)

#ifdef CGAL_CPPF
   : rat(d)
#endif
  {
    // std::cout << "\ndouble = " << d << std::endl;
    using boost::uint64_t;
    union {
#ifdef CGAL_LITTLE_ENDIAN
      struct { uint64_t man:52; uint64_t exp:11; uint64_t sig:1; } s;
#else /* CGAL_BIG_ENDIAN */
      //WARNING: untested!
      struct { uint64_t sig:1; uint64_t exp:11; uint64_t man:52; } s;
#endif
      double d;
    } u;
    u.d = d;

    uint64_t m;
    uint64_t dexp = u.s.exp;
    CGAL_assertion_msg(dexp != 2047, "Creating an cpp_float from infinity or NaN.");
    if (dexp == 0) {
      if (d == 0) { exp=0; return; }
      else { // denormal number
        m = u.s.man;
        ++dexp;
      }
    } else {
      m = (1LL << 52) |  u.s.man;
    }


    int idexp = (int)dexp;
    idexp -= 1023;

    // std::cout << "m    = "  << m << std::endl;
    // std::cout << "idexp = " << idexp << std::endl;

    int shifted = internal::low_bit(m);

    m >>= shifted;

    int nbits = internal::high_bit(m);
    // std::cout << "nbits = " << nbits << std::endl;

    exp = idexp - nbits;
    man = m;
    if(u.s.sig){
      man = -man;
    }
#ifdef CGAL_CPPF
    assert(rat.sign() == man.sign());
#endif
    // std::cout << "m = " << m << " * 2^" << exp  << std::endl;
    // fmt(m);
  }

  friend std::ostream& operator<<(std::ostream& os, const cpp_float& m)
  {
    return os << m.man << " * 2 ^ " << m.exp << " ( " << m.to_double() << ") "
#ifdef CGAL_CPPF
              << "  " << m.rat
#endif
      ;
  }


  friend cpp_float operator-(cpp_float const&x)
  {
#ifdef CGAL_CPPF
    return cpp_float(-x.man,x.exp, -x.rat);
#else
    return cpp_float(-x.man,x.exp);
#endif
  }

  cpp_float& operator*=(const cpp_float& other)
  {
#ifdef CGAL_CPPF
    rat *= other.rat;
#endif
    man *= other.man;
    exp += other.exp;
    return *this;
  }


  friend
  cpp_float operator*(const cpp_float& a, const cpp_float&b){
#ifdef CGAL_CPPF
    return cpp_float(a.man*b.man, a.exp+b.exp, a.rat * b.rat);
#else
    return cpp_float(a.man*b.man, a.exp+b.exp);
#endif
  }


  cpp_float operator+=(const cpp_float& other)
  {
#ifdef CGAL_CPPF
    rat += other.rat;
#endif
    int shift = exp - other.exp;
    if(shift > 0){
      man = (man << shift) + other.man;
      exp = other.exp;
    }else if(shift < 0){
      man = man + (other.man << -shift);
    }else{
      man += other.man;
    }
    return *this;
  }


#ifdef CGAL_CPPF
  friend
  cpp_float operator+(const cpp_float& a, const cpp_float&b){
    int shift = a.exp - b.exp;
    if(shift > 0){
      return cpp_float((a.man << shift) + b.man, b.exp, a.rat+b.rat);
    }else if(shift < 0){
      return cpp_float(a.man + (b.man << -shift), a.exp, a.rat+b.rat);
    }
    return cpp_float(a.man + b.man, a.exp, a.rat+b.rat);
  }
#else
  friend
  cpp_float operator+(const cpp_float& a, const cpp_float&b){
    int shift = a.exp - b.exp;
    if(shift > 0){
      return cpp_float((a.man << shift) + b.man, b.exp);
    }else if(shift < 0){
      return cpp_float(a.man + (b.man << -shift), a.exp);
    }
    return cpp_float(a.man + b.man, a.exp);
  }
#endif


  cpp_float operator-=(const cpp_float& other)
  {

#ifdef CGAL_CPPF
    rat -= other.rat;
#endif
    int shift = exp - other.exp;
    if(shift > 0){
      man <<= shift;
      man -= other.man;
      exp = other.exp;
    }else if(shift < 0){
      man -= (other.man << -shift);
    }else{
      man -= other.man;
    }
    return *this;
  }

  #ifdef CGAL_CPPF
  friend
  cpp_float operator-(const cpp_float& a, const cpp_float&b){

    int shift = a.exp - b.exp;
    if(shift > 0){
      return cpp_float((a.man << shift) - b.man, b.exp, a.rat-b.rat);
    }else if(shift < 0){
      return cpp_float(a.man - (b.man << -shift), a.exp, a.rat-b.rat);
    }
    return cpp_float(a.man - b.man, a.exp, a.rat-b.rat);
  }
#else
  friend
  cpp_float operator-(const cpp_float& a, const cpp_float&b){

    int shift = a.exp - b.exp;
    if(shift > 0){
      return cpp_float((a.man << shift) - b.man, b.exp);
    }else if(shift < 0){
      return cpp_float(a.man - (b.man << -shift), a.exp);
    }
    return cpp_float(a.man - b.man, a.exp);
  }
#endif

  bool is_positive() const
  {
    return CGAL::is_positive(man);
  }

  bool is_negative() const
  {
    return CGAL::is_negative(man);
  }

  friend bool operator<(const cpp_float& a, const cpp_float& b)
  {
    if(((! a.is_positive()) && b.is_positive())
       || a.is_negative()&& b.is_zero())return true;
    if(((! b.is_positive()) && a.is_positive())
       ||b.is_negative()&& a.is_zero())return false;

#ifdef CGAL_CPPF
    bool qres = a.rat < b.rat;
#endif
    cpp_float d = b-a;
#ifdef CGAL_CPPF
    assert(qres == d.is_positive());
#endif
    return d.is_positive();
  }

  friend bool operator>(cpp_float const&a, cpp_float const&b){
    return b<a;
  }
  friend bool operator>=(cpp_float const&a, cpp_float const&b){
    return !(a<b);
  }
  friend bool operator<=(cpp_float const&a, cpp_float const&b){
    return !(a>b);
  }


  friend bool operator==(cpp_float const&a, cpp_float const&b){

#ifdef CGAL_CPPF
    bool qres = a.rat == b.rat;
#endif
   int shift = a.exp - b.exp;
    if(shift > 0){
      Mantissa ac(a.man);
      ac <<= shift;
#ifdef CGAL_CPPF
      assert( qres == (ac == b.man));
#endif
      return ac == b.man;
    }else if(shift < 0){
      Mantissa  bc(b.man);
      bc <<= -shift;
#ifdef CGAL_CPPF
      assert(qres == (a.man == bc));
#endif
      return a.man == bc;
    }
#ifdef CGAL_CPPF
    assert(qres == (a.man == b.man));
#endif
    return a.man==b.man;
  }

  Comparison_result compare(const cpp_float& other) const
  {
    if(*this < other) return SMALLER;
    if(*this > other) return LARGER;
    return EQUAL;
  }

  friend bool operator!=(cpp_float const&a, cpp_float const&b){
    return !(a==b);
  }

  double to_double() const
  {
    if(exp == 0){
      return CGAL::to_double(man);
    }
    if(exp > 0){
      Mantissa as(man);
      as <<= exp;
      return CGAL::to_double(as);
    }
    Mantissa pow(1);
    pow <<= -exp;
    boost::multiprecision::cpp_rational r(man, pow);
    return CGAL::to_double(r);
  }

  std::pair<double,double> to_interval() const
  {
      if(exp == 0){
      return CGAL::to_interval(man);
    }
    if(exp > 0){
      Mantissa as = man << exp;
      return CGAL::to_interval(as);
    }
    Mantissa pow(1);
    pow <<= -exp;
    boost::multiprecision::cpp_rational r(man, pow);
    return CGAL::to_interval(r);
  }


  bool is_zero () const {
    return CGAL::is_zero(man);
  }


  bool is_one () const {
    return *this == cpp_float(1);
  }


  CGAL::Sign sign () const
  {
    return CGAL::sign(man);
  }

};


  template <> struct Algebraic_structure_traits< cpp_float >
    : public Algebraic_structure_traits_base< cpp_float, Integral_domain_without_division_tag >  {
      typedef Tag_true            Is_exact;
      typedef Tag_false            Is_numerical_sensitive;

      struct Is_zero
        : public CGAL::cpp98::unary_function< Type, bool > {
          bool operator()( const Type& x ) const {
            return x.is_zero();
          }
        };

      struct Is_one
        : public CGAL::cpp98::unary_function< Type, bool > {
          bool operator()( const Type& x ) const {
            return x.is_one();
          }
        };

      struct Gcd
        : public CGAL::cpp98::binary_function< Type, Type, Type > {
          Type operator()(
              const Type& x,
              const Type& y ) const {
            assert(false);
            return Type(); // cpp_float_gcd(x, y);
          }
        };

      struct Square
        : public CGAL::cpp98::unary_function< Type, Type > {
          Type operator()( const Type& x ) const {
            return x*x ; // cpp_float_square(x);
          }
        };

      struct Integral_division
        : public CGAL::cpp98::binary_function< Type, Type, Type > {
          Type operator()(
              const Type& x,
              const Type& y ) const {
            assert(false);
            return Type(); // x / y;
          }
        };

      struct Sqrt
        : public CGAL::cpp98::unary_function< Type, Type > {
          Type operator()( const Type& x) const {
            assert(false);
            return Type(); // cpp_float_sqrt(x);
          }
        };

      struct Is_square
        : public CGAL::cpp98::binary_function< Type, Type&, bool > {
          bool operator()( const Type& x, Type& y ) const {
            // TODO: avoid doing 2 calls.
            assert(false);
            return true;
          }
          bool operator()( const Type& x) const {
            assert(false);
            return true;
          }
        };

    };
  template <> struct Real_embeddable_traits< cpp_float >
    : public INTERN_RET::Real_embeddable_traits_base< cpp_float , CGAL::Tag_true > {
      struct Sgn
        : public CGAL::cpp98::unary_function< Type, ::CGAL::Sign > {
          ::CGAL::Sign operator()( const Type& x ) const {
            return x.sign();
          }
        };

      struct To_double
        : public CGAL::cpp98::unary_function< Type, double > {
            double operator()( const Type& x ) const {
              return x.to_double();
            }
        };

      struct Compare
        : public CGAL::cpp98::binary_function< Type, Type, Comparison_result > {
            Comparison_result operator()(
                const Type& x,
                const Type& y ) const {
              return x.compare(y);
            }
        };

      struct To_interval
        : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
            std::pair<double, double> operator()( const Type& x ) const {
              return x.to_interval();
            }
        };

    };


CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(cpp_float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,cpp_float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,cpp_float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,cpp_float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float    ,cpp_float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double   ,cpp_float)

} // namespace CGAL


#endif // CGAL_CPP_FLOAT_H
