
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


#include <CGAL/boost_mp_type.h>
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
    return (int)ret; // AF:  was 63 - (int)ret;  The others have to be changed too
#elif defined(__xlC__)
    // Macro supposedly not defined on z/OS.
    return __cntlz8 (x);
#else
    return __builtin_clzll (x);
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

  boost::multiprecision::cpp_int man;
  int exp; /* The number man (an integer) * 2 ^ exp  */

public:
  cpp_float()
    : man(), exp()
  {}

  cpp_float(int i)
    : man(i),exp(0)
  {}

  cpp_float(const boost::multiprecision::cpp_int& man,  int exp)
    : man(man),exp(exp)
  {}


  cpp_float(double d)
  {
      std::cout << "\ndouble = " << d << std::endl;
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

    std::cout << "m    = "  << m << std::endl;
    std::cout << "idexp = " << idexp << std::endl;

    int shifted = internal::low_bit(m);

    m >>= shifted;

    int nbits = internal::high_bit(m);
    std::cout << "nbits = " << nbits << std::endl;

    exp = idexp - nbits;
    if(u.s.sig){
      m = -m;
    }
    std::cout << "m = " << m << " * 2^" << exp  << std::endl;
    // fmt(m);
    man = boost::multiprecision::cpp_int(m);
  }

  friend std::ostream& operator<<(std::ostream& os, const cpp_float& m)
  {
    return os << m.man << " * 2 ^ " << m.exp << " ( " << m.to_double() << ") ";
  }

  friend cpp_float operator-(cpp_float const&x)
  {
    return cpp_float(-x.man,x.exp);
  }

  cpp_float& operator*=(const cpp_float& other)
  {
    man *= other.man;
    exp += other.exp;
    return *this;
  }

  cpp_float operator+=(const cpp_float& other)
  {
    int shift = exp - other.exp;
    if(shift > 0){
      man <<= shift;
      man += other.man;
      exp = other.exp;
    }else if(shift < 0){
      boost::multiprecision::cpp_int cpy(other.man);
      cpy << shift;
      man += cpy;
    }else{
      man += other.man;
    }
    return *this;
  }

  cpp_float operator-=(const cpp_float& other)
  {
    int shift = exp - other.exp;
    if(shift > 0){
      man <<= shift;
      man -= other.man;
      exp = other.exp;
    }else if(shift < 0){
      boost::multiprecision::cpp_int cpy(other.man);
      cpy << shift;
      man -= cpy;
    }else{
      man -= other.man;
    }
    return *this;
  }

  bool positive() const
  {
    return is_positive(man);
  }


  friend bool operator<(const cpp_float& a, const cpp_float& b)
  {
    cpp_float d(b);
    d -= a;
    return d.positive();
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
   int shift = a.exp - b.exp;
    if(shift > 0){
      boost::multiprecision::cpp_int ac(a.man);
      ac <<= shift;
      return ac == b.man;
    }else if(shift < 0){
      boost::multiprecision::cpp_int bc(b.man);
      bc <<= shift;
      return a.man == bc;
    }
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
      boost::multiprecision::cpp_int as(man);
      as <<= exp;
      return CGAL::to_double(as);
    }
    boost::multiprecision::cpp_int pow(1);
    pow <<= -exp;
    boost::multiprecision::cpp_rational rat(man, pow);
    return CGAL::to_double(rat);
  }

  std::pair<double,double> to_interval() const
  {
    assert(false);
    double zero(0);
    return std::make_pair(zero,zero);
  }


  bool is_zero () const {
    assert(false);
    return man==0 && exp == 0;
  }


  bool is_one () const {
    assert(false);
    return true;
    //  return exp==0 && size==1 && data()[0]==1;
  }


  CGAL::Sign sign () const
  {
    return CGAL::sign(man);
  }

};

  inline
  cpp_float operator+(const cpp_float& a, const cpp_float&b){
    cpp_float ret(a);
    return ret += b;
  }

  inline
  cpp_float operator-(const cpp_float& a, const cpp_float&b){
    cpp_float ret(a);
    return ret -= b;
  }

  inline
  cpp_float operator*(const cpp_float& a, const cpp_float&b){
    cpp_float ret(a);
    return ret *= b;
  }


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
            assert(false);
            return false; // x.is_one();
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
