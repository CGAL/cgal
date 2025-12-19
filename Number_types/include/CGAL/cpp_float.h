// Copyright (c) 2023 GeometryFactory (France), INRIA Saclay - Ile de France (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)        :  Andreas Fabri, Marc Glisse

#ifndef CGAL_CPP_FLOAT_H
#define CGAL_CPP_FLOAT_H

//#define CGAL_DEBUG_CPPF

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

// It would have made sense to make this
// boost::multiprecision::number<some_new_backend>, but we keep that
// for later when we contribute to boost::mp

class cpp_float {
#ifdef CGAL_DEBUG_CPPF
  boost::multiprecision::cpp_rational rat;
#endif

  typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<512> > Mantissa;
  Mantissa man;
  int exp; /* The number man (an integer) * 2 ^ exp  */

public:
  cpp_float()
    :
#ifdef CGAL_DEBUG_CPPF
    rat(),
#endif
    man(), exp()
  {}

  cpp_float(short i)
    :
#ifdef CGAL_DEBUG_CPPF
    rat(i),
#endif
    man(i),exp(0)
  {}

  cpp_float(int i)
    :
#ifdef CGAL_DEBUG_CPPF
    rat(i),
#endif
    man(i),exp(0)
  {}

  cpp_float(long i)
    :
#ifdef CGAL_DEBUG_CPPF
    rat(i),
#endif
    man(i),exp(0)
  {}
#ifdef CGAL_DEBUG_CPPF
  cpp_float(const Mantissa& man,  int exp, const boost::multiprecision::cpp_rational& rat)
    : rat(rat), man(man),exp(exp)
  {}
#else

  cpp_float(const Mantissa& man, int exp)
      : man(man), exp(exp)
  {
    CGAL_HISTOGRAM_PROFILER("size (man/exp)", static_cast<int>(man.backend().size()));
  }

#ifndef CGAL_DEBUG_CPPF
  template <typename Expression>
  cpp_float(const Expression& man, int exp)
    :man(man), exp(exp)
  {
    CGAL_HISTOGRAM_PROFILER("size (expression/exp)", static_cast<int>(this->man.backend().size()));
  }
#endif

#endif
  cpp_float(double d)

#ifdef CGAL_DEBUG_CPPF
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
      man.backend().negate();
    }
#ifdef CGAL_DEBUG_CPPF
    CGAL_assertion(rat.sign() == man.sign());
#endif
    // std::cout << "m = " << m << " * 2^" << exp  << std::endl;
    // fmt(m);

    CGAL_HISTOGRAM_PROFILER("size when constructed from double", static_cast<int>(man.backend().size()));
  }


  friend std::ostream& operator<<(std::ostream& os, const cpp_float& m)
  {
    return os << m.to_double();
#if 0 // dehug output
    return os << m.man << " * 2 ^ " << m.exp << " ( " << m.to_double() << ") "
#ifdef CGAL_DEBUG_CPPF
              << "  " << m.rat
#endif
      ;
#endif
  }


  friend std::istream& operator>>(std::istream& is, cpp_float& m)
  {
    double d;
    is >> d;
    m = cpp_float(d);
    return is;
  }


  friend cpp_float operator-(cpp_float const&x)
  {
#ifdef CGAL_DEBUG_CPPF
    return cpp_float(-x.man,x.exp, -x.rat);
#else
    return cpp_float(-x.man,x.exp);
#endif
  }

  cpp_float& operator*=(const cpp_float& other)
  {
#ifdef CGAL_DEBUG_CPPF
    rat *= other.rat;
#endif
    man *= other.man;
    exp += other.exp;
    return *this;
  }

  friend
  cpp_float operator*(const cpp_float& a, const cpp_float&b){
#ifdef CGAL_DEBUG_CPPF
    return cpp_float(a.man*b.man, a.exp+b.exp, a.rat * b.rat);
#else
      return cpp_float(a.man*b.man, a.exp+b.exp);
#endif
  }

  // Marc Glisse commented on github:
  // We can sometimes end up with a mantissa that has quite a few zeros at the end,
  // but the cases where the mantissa is too long by more than 1 limb should be negligible,
  // and normalizing so the mantissa is always odd would have a cost.
  cpp_float operator+=(const cpp_float& other)
  {
#ifdef CGAL_DEBUG_CPPF
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

#ifdef CGAL_DEBUG_CPPF
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
    CGAL_HISTOGRAM_PROFILER("shift+", CGAL::abs(shift));
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

#ifdef CGAL_DEBUG_CPPF
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

  #ifdef CGAL_DEBUG_CPPF
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
    CGAL_HISTOGRAM_PROFILER("shift-", CGAL::abs(shift));
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

  // Would it make sense to compare the sign of the exponent?
  // to distinguish the interval between ]-1,1[  from the rest.
  friend bool operator<(const cpp_float& a, const cpp_float& b)
  {
    if(((! a.is_positive()) && b.is_positive())
       || (a.is_negative() && b.is_zero()))return true;
    if(((! b.is_positive()) && a.is_positive())
       || (b.is_negative() && a.is_zero()))return false;

#ifdef CGAL_DEBUG_CPPF
    bool qres = a.rat < b.rat;
#endif
    cpp_float d = b-a;
#ifdef CGAL_DEBUG_CPPF
    CGAL_assertion(qres == d.is_positive());
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

#ifdef CGAL_DEBUG_CPPF
    bool qres = a.rat == b.rat;
#endif
   int shift = a.exp - b.exp;
   CGAL_HISTOGRAM_PROFILER("shift==", CGAL::abs(shift));
    if(shift > 0){
      Mantissa ac = a.man << shift;
#ifdef CGAL_DEBUG_CPPF
      CGAL_assertion( qres == (ac == b.man));
#endif
      return ac == b.man;
    }else if(shift < 0){
      Mantissa  bc = b.man << -shift;
#ifdef CGAL_DEBUG_CPPF
      CGAL_assertion(qres == (a.man == bc));
#endif
      return a.man == bc;
    }
#ifdef CGAL_DEBUG_CPPF
    CGAL_assertion(qres == (a.man == b.man));
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
    if(! is_positive()) return false;

    int msb = static_cast<int>(boost::multiprecision::msb(man));
    if (msb != -exp) return false;
    int lsb = static_cast<int>(boost::multiprecision::lsb(man));
    return (msb == lsb);
  }


  CGAL::Sign sign () const
  {
    return CGAL::sign(man);
  }

  std::size_t size() const
  {
    return man.backend().size();
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
              const Type& /* x */,
              const Type& /* y */ ) const {
            CGAL_assertion(false);
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
              const Type& /* x */,
              const Type& /* y */ ) const {
            CGAL_assertion(false);
            return Type(); // x / y;
          }
        };

      struct Sqrt
        : public CGAL::cpp98::unary_function< Type, Type > {
          Type operator()( const Type& /* x */) const {
            CGAL_assertion(false);
            return Type(); // cpp_float_sqrt(x);
          }
        };

      struct Is_square
        : public CGAL::cpp98::binary_function< Type, Type&, bool > {
          bool operator()( const Type& /* x */, Type& /* y */ ) const {
            // TODO: avoid doing 2 calls.
            CGAL_assertion(false);
            return true;
          }
          bool operator()( const Type& /* x */) const {
            CGAL_assertion(false);
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
