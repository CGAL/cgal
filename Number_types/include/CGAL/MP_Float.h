// Copyright (c) 2001-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_MP_FLOAT_H
#define CGAL_MP_FLOAT_H

#include <CGAL/number_type_basic.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Quotient.h>
#include <CGAL/Sqrt_extension.h>

#include <CGAL/utils.h>

#include <CGAL/Interval_nt.h>
#include <iostream>
#include <vector>
#include <algorithm>

// MP_Float : multiprecision scaled integers.

// Some invariants on the internal representation :
// - zero is represented by an empty vector, and whatever exp.
// - no leading or trailing zero in the vector => unique

// The main algorithms are :
// - Addition/Subtraction
// - Multiplication
// - Integral division div(), gcd(), operator%().
// - Comparison
// - to_double() / to_interval()
// - Construction from a double.
// - IOs

// TODO :
// - The exponent really overflows sometimes -> make it multiprecision.
// - Write a generic wrapper that adds an exponent to be used by MP integers.
// - Karatsuba (or other) ?  Would be fun to implement at least.
// - Division, sqrt... : different options :
//   - nothing
//   - convert to double, take approximation, compute over double, reconstruct

namespace CGAL {

class MP_Float;

template < typename > class Quotient; // Needed for overloaded To_double

namespace INTERN_MP_FLOAT {

Comparison_result compare(const MP_Float&, const MP_Float&);

MP_Float square(const MP_Float&);

// to_double() returns, not the closest double, but a one bit error is allowed.
// We guarantee : to_double(MP_Float(double d)) == d.

double to_double(const MP_Float&);

double to_double(const Quotient<MP_Float>&);

std::pair<double,double> to_interval(const MP_Float &);

std::pair<double,double> to_interval(const Quotient<MP_Float>&);

MP_Float div(const MP_Float& n1, const MP_Float& n2);

MP_Float gcd(const MP_Float& a, const MP_Float& b);
  
} //namespace INTERN_MP_FLOAT

std::pair<double, int>
to_double_exp(const MP_Float &b);

// Returns (first * 2^second), an interval surrounding b.
std::pair<std::pair<double, double>, int>
to_interval_exp(const MP_Float &b);

std::ostream &
operator<< (std::ostream & os, const MP_Float &b);

// This one is for debug.
std::ostream &
print (std::ostream & os, const MP_Float &b);

std::istream &
operator>> (std::istream & is, MP_Float &b);

MP_Float operator+(const MP_Float &a, const MP_Float &b);

MP_Float operator-(const MP_Float &a, const MP_Float &b);

MP_Float operator*(const MP_Float &a, const MP_Float &b);

MP_Float operator%(const MP_Float &a, const MP_Float &b);


class MP_Float
{
public:
  typedef short          limb;
  typedef unsigned short unsigned_limb;
  typedef int            limb2;
  typedef double         exponent_type;

  typedef std::vector<limb>  V;
  typedef V::const_iterator  const_iterator;
  typedef V::iterator        iterator;

private:

  void remove_leading_zeros()
  {
    while (!v.empty() && v.back() == 0)
      v.pop_back();
  }

  void remove_trailing_zeros()
  {
    if (v.empty() || v.front() != 0)
      return;

    iterator i = v.begin();
    for (++i; *i == 0; ++i)
      ;
    exp += static_cast<exponent_type>(i-v.begin());
    v.erase(v.begin(), i);
  }

  // The constructors from float/double/long_double are factorized in the
  // following template :
  template < typename T >
  void construct_from_builtin_fp_type(T d);

public:
#ifdef CGAL_ROOT_OF_2_ENABLE_HISTOGRAM_OF_NUMBER_OF_DIGIT_ON_THE_COMPLEX_CONSTRUCTOR
    int tam() const { return v.size(); }
#endif

  // Splits a limb2 into 2 limbs (high and low).
  static
  void split(limb2 l, limb & high, limb & low)
  {
    const unsigned int sizeof_limb=8*sizeof(limb);
    const limb2 mask = 0x0000ffff;
   
    //Note: For Integer type, if the destination type is signed, the value is unchanged 
    //if it can be represented in the destination type)
    low = static_cast<limb>(l & mask); //extract low bits from l
    high= static_cast<limb>((l - low) >> sizeof_limb); //extract high bits from l
    
    CGAL_postcondition ( l == low + ( static_cast<limb2>(high) << sizeof_limb ) );
  }

  // Given a limb2, returns the higher limb.
  static
  limb higher_limb(limb2 l)
  {
      limb high, low;
      split(l, high, low);
      return high;
  }

  void canonicalize()
  {
    remove_leading_zeros();
    remove_trailing_zeros();
  }

  MP_Float()
      : exp(0)
  {
    CGAL_assertion(sizeof(limb2)==4); // so that the above 0x0000ffff is correct
    CGAL_assertion(sizeof(limb2) == 2*sizeof(limb));
    CGAL_assertion(v.empty());
    // Creates zero.
  }

#if 0
  // Causes ambiguities
  MP_Float(limb i)
  : v(1,i), exp(0)
  {
    remove_leading_zeros();
  }
#endif

  MP_Float(limb2 i)
  : v(2), exp(0)
  {
    split(i, v[1], v[0]);
    canonicalize();
  }

  MP_Float(float d);

  MP_Float(double d);

  MP_Float(long double d);

  MP_Float operator+() const {
    return *this;
  }

  MP_Float operator-() const
  {
    return MP_Float() - *this;
  }

  MP_Float& operator+=(const MP_Float &a) { return *this = *this + a; }
  MP_Float& operator-=(const MP_Float &a) { return *this = *this - a; }
  MP_Float& operator*=(const MP_Float &a) { return *this = *this * a; }
  MP_Float& operator%=(const MP_Float &a) { return *this = *this % a; }

  exponent_type max_exp() const
  {
    return exponent_type(v.size()) + exp;
  }

  exponent_type min_exp() const
  {
    return exp;
  }

  limb of_exp(exponent_type i) const
  {
    if (i < exp || i >= max_exp())
      return 0;
    return v[static_cast<int>(i-exp)];
  }

  bool is_zero() const
  {
    return v.empty();
  }

  Sign sign() const
  {
    if (v.empty())
      return ZERO;
    if (v.back()>0)
      return POSITIVE;
    CGAL_assertion(v.back()<0);
    return NEGATIVE;
  }

  void clear()
  {
    v.clear();
    exp = 0;
  }

  void swap(MP_Float &m)
  {
    std::swap(v, m.v);
    std::swap(exp, m.exp);
  }

  // Converts to a rational type (e.g. Gmpq).
  template < typename T >
  T to_rational() const
  {
    const unsigned log_limb = 8 * sizeof(MP_Float::limb);

    if (is_zero())
      return 0;

    MP_Float::const_iterator i;
    exponent_type exp2 = min_exp() * log_limb;
    T res = 0;

    for (i = v.begin(); i != v.end(); ++i)
    {
      res += T(std::ldexp(static_cast<double>(*i),static_cast<int>(exp2)));
      exp2 += log_limb;
    }

    return res;
  }

  std::size_t size() const
  {
    return v.size();
  }

  // Returns a scaling factor (in limbs) which would be good to extract to get
  // a value with an exponent close to 0.
  exponent_type find_scale() const
  {
    return exp + exponent_type(v.size());
  }

  // Rescale the value by some factor (in limbs).  (substract the exponent)
  void rescale(exponent_type scale)
  {
    if (v.size() != 0)
      exp -= scale;
  }

  // Accessory function that finds the least significant bit set (its position).
  static unsigned short 
  lsb(limb l)
  {
    unsigned short nb = 0;
    for (; (l&1)==0; ++nb, l=(limb)(l>>1) )
      ;
    return nb;
  }

  // This one is needed for normalizing gcd so that the mantissa is odd
  // and non-negative, and the exponent is 0.
  void gcd_normalize()
  {
    const unsigned log_limb = 8 * sizeof(MP_Float::limb);
    if (is_zero())
      return;
    // First find how many least significant bits are 0 in the last digit.
    unsigned short nb = lsb(v[0]);
    if (nb != 0)
      *this = *this * (1<<(log_limb-nb));
    CGAL_assertion((v[0]&1) != 0);
    exp=0;
    if (sign() == NEGATIVE)
      *this = - *this;
  }

  MP_Float unit_part() const
  {
    if (is_zero())
      return 1;
    MP_Float r = (sign() > 0) ? *this : - *this;
    CGAL_assertion(r.v.begin() != r.v.end());
    unsigned short nb = lsb(r.v[0]);
    r.v.clear();
    r.v.push_back((limb)(1<<nb));
    return (sign() > 0) ? r : -r;
  }

  bool is_integer() const
  {
    return is_zero() || (exp >= 0);
  }

  V v;
  exponent_type exp;
};

namespace internal{
std::pair<MP_Float, MP_Float> // <quotient, remainder>
division(const MP_Float & n, const MP_Float & d);
} // namespace internal

inline
void swap(MP_Float &m, MP_Float &n)
{ m.swap(n); }

inline
bool operator<(const MP_Float &a, const MP_Float &b)
{ return INTERN_MP_FLOAT::compare(a, b) == SMALLER; }

inline
bool operator>(const MP_Float &a, const MP_Float &b)
{ return b < a; }

inline
bool operator>=(const MP_Float &a, const MP_Float &b)
{ return ! (a < b); }

inline
bool operator<=(const MP_Float &a, const MP_Float &b)
{ return ! (a > b); }

inline
bool operator==(const MP_Float &a, const MP_Float &b)
{ return (a.v == b.v) && (a.v.empty() || (a.exp == b.exp)); }

inline
bool operator!=(const MP_Float &a, const MP_Float &b)
{ return ! (a == b); }

MP_Float
approximate_sqrt(const MP_Float &d);

MP_Float
approximate_division(const MP_Float &n, const MP_Float &d);



// Algebraic structure traits specialization
template <> class Algebraic_structure_traits< MP_Float >
  : public Algebraic_structure_traits_base< MP_Float,
                                            Unique_factorization_domain_tag
	    // with some work on mod/div it could be Euclidean_ring_tag
                                          >  {
  public:

    typedef Tag_true            Is_exact;
    typedef Tag_true            Is_numerical_sensitive;

    struct Unit_part
      : public std::unary_function< Type , Type >
    {
      Type operator()(const Type &x) const {
        return x.unit_part();
      }
    };

    struct Integral_division
        : public std::binary_function< Type,
                                 Type,
                                 Type > {
    public:
        Type operator()(
                const Type& x,
                const Type& y ) const {
            std::pair<MP_Float, MP_Float> res = internal::division(x, y);
            CGAL_assertion_msg(res.second == 0,
                "exact_division() called with operands which do not divide");
            return res.first;
        }
    };


    class Square
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return INTERN_MP_FLOAT::square(x);
        }
    };

    class Gcd
      : public std::binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y ) const {
          return INTERN_MP_FLOAT::gcd( x, y );
        }
    };

    class Div
      : public std::binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y ) const {
          return INTERN_MP_FLOAT::div( x, y );
        }
    };

  typedef INTERN_AST::Mod_per_operator< Type > Mod;
// Default implementation of Divides functor for unique factorization domains
  // x divides y if gcd(y,x) equals x up to inverses 
  class Divides 
    : public std::binary_function<Type,Type,bool>{ 
  public:
    bool operator()( const Type& x,  const Type& y) const {  
      return internal::division(y,x).second == 0 ;
    }
    // second operator computing q = x/y 
    bool operator()( const Type& x,  const Type& y, Type& q) const {    
      std::pair<Type,Type> qr = internal::division(y,x);
      q=qr.first;
      return qr.second == 0;
      
    }
    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type , bool)
  };
};

// Real embeddable traits
template <> class Real_embeddable_traits< MP_Float >
  : public INTERN_RET::Real_embeddable_traits_base< MP_Float , CGAL::Tag_true > {
  public:

    class Sgn
      : public std::unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
          return x.sign();
        }
    };

    class Compare
      : public std::binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return INTERN_MP_FLOAT::compare( x, y );
        }
    };

    class To_double
      : public std::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          return INTERN_MP_FLOAT::to_double( x );
        }
    };

    class To_interval
      : public std::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          return INTERN_MP_FLOAT::to_interval( x );
        }
    };
};



namespace INTERN_MP_FLOAT{

//Sqrt_extension internally uses Algebraic_structure_traits
template <class ACDE_TAG_, class FP_TAG>  
double
to_double(const Sqrt_extension<MP_Float,MP_Float,ACDE_TAG_,FP_TAG> &x)
{
  typedef MP_Float RT;
  typedef Quotient<RT> FT;
  typedef CGAL::Rational_traits< FT > Rational;
  Rational r;
  const RT r1 = r.numerator(x.a0());
  const RT d1 = r.denominator(x.a0());

  if(x.is_rational()) {
    std::pair<double, int> n = to_double_exp(r1);
    std::pair<double, int> d = to_double_exp(d1);
    double scale = std::ldexp(1.0, n.second - d.second);
    return (n.first / d.first) * scale;
  }

  const RT r2 = r.numerator(x.a1());
  const RT d2 = r.denominator(x.a1());
  const RT r3 = r.numerator(x.root());
  const RT d3 = r.denominator(x.root());

  std::pair<double, int> n1 = to_double_exp(r1);
  std::pair<double, int> v1 = to_double_exp(d1);
  double scale1 = std::ldexp(1.0, n1.second - v1.second);

  std::pair<double, int> n2 = to_double_exp(r2);
  std::pair<double, int> v2 = to_double_exp(d2);
  double scale2 = std::ldexp(1.0, n2.second - v2.second);

  std::pair<double, int> n3 = to_double_exp(r3);
  std::pair<double, int> v3 = to_double_exp(d3);
  double scale3 = std::ldexp(1.0, n3.second - v3.second);

  return ((n1.first / v1.first) * scale1) + 
         ((n2.first / v2.first) * scale2) *
         std::sqrt((n3.first / v3.first) * scale3);
}

} //namespace INTERN_MP_FLOAT


namespace internal {
// This compares the absolute values of the odd-mantissa.
// (take the mantissas, get rid of all powers of 2, compare
// the absolute values)
inline
Sign
compare_bitlength(const MP_Float &a, const MP_Float &b)
{
  if (a.is_zero())
    return b.is_zero() ? EQUAL : SMALLER;
  if (b.is_zero())
    return LARGER;

  //Real_embeddable_traits<MP_Float>::Abs abs;

  MP_Float aa = CGAL_NTS abs(a);
  MP_Float bb = CGAL_NTS abs(b);

  if (aa.size() > (bb.size() + 2)) return LARGER;
  if (bb.size() > (aa.size() + 2)) return SMALLER;

  // multiply by 2 till last bit is 1.
  while (((aa.v[0]) & 1) == 0) // last bit is zero
    aa = aa + aa;

  while (((bb.v[0]) & 1) == 0) // last bit is zero
    bb = bb + bb;

  // sizes might have changed
  if (aa.size() > bb.size()) return LARGER;
  if (aa.size() < bb.size()) return SMALLER;

  for (std::size_t i = aa.size(); i > 0; --i)
  {
    if (aa.v[i-1] > bb.v[i-1]) return LARGER;
    if (aa.v[i-1] < bb.v[i-1]) return SMALLER;
  }
  return EQUAL;
}

inline // Move it to libCGAL once it's stable.
std::pair<MP_Float, MP_Float> // <quotient, remainder>
division(const MP_Float & n, const MP_Float & d)
{
  typedef MP_Float::exponent_type  exponent_type;

  MP_Float remainder = n, divisor = d;

  CGAL_precondition(divisor != 0);

  // Rescale d to have a to_double() value with reasonnable exponent.
  exponent_type scale_d = divisor.find_scale();
  divisor.rescale(scale_d);
  const double dd = INTERN_MP_FLOAT::to_double(divisor);

  MP_Float res = 0;
  exponent_type scale_remainder = 0;

  bool first_time_smaller_than_divisor = true;

  // School division algorithm.

  while ( remainder != 0 )
  {
    // We have to rescale, since remainder can diminish towards 0.
    exponent_type tmp_scale = remainder.find_scale();
    remainder.rescale(tmp_scale);
    res.rescale(tmp_scale);
    scale_remainder += tmp_scale;

    // Compute a double approximation of the quotient
    // (imagine school division with base ~2^53).
    double approx = INTERN_MP_FLOAT::to_double(remainder) / dd;
    CGAL_assertion(approx != 0);
    res += approx;
    remainder -= approx * divisor;

    if (remainder == 0)
      break;

    // Then we need to fix it up by checking if neighboring double values
    // are closer to the exact result.
    // There should not be too many iterations, because approx is only a few ulps
    // away from the optimal.
    // If we don't do the fixup, then spurious bits can be introduced, which
    // will require an unbounded amount of additional iterations to be eliminated.

    // The direction towards which we need to try to move from "approx".
    double direction = (CGAL_NTS sign(remainder) == CGAL_NTS sign(dd))
                     ?  std::numeric_limits<double>::infinity()
                     : -std::numeric_limits<double>::infinity();

    while (true)
    {
      const double approx2 = nextafter(approx, direction);
      const double delta = approx2 - approx;
      MP_Float new_remainder = remainder - delta * divisor;
      if (CGAL_NTS abs(new_remainder) < CGAL_NTS abs(remainder)) {
        remainder = new_remainder;
        res += delta;
        approx = approx2;
      }
      else {
        break;
      }
    }

    if (remainder == 0)
      break;

    // Test condition for non-exact division (with remainder).
    if (compare_bitlength(remainder, divisor) == SMALLER)
    {
      if (! first_time_smaller_than_divisor)
      {
        // Scale back.
        res.rescale(scale_d - scale_remainder);
        remainder.rescale(- scale_remainder);
        CGAL_postcondition(res * d  + remainder == n);
        return std::make_pair(res, remainder);
      }
      first_time_smaller_than_divisor = false;
    }
  }

  // Scale back the result.
  res.rescale(scale_d - scale_remainder);
  CGAL_postcondition(res * d == n);
  return std::make_pair(res, MP_Float(0));
}

inline // Move it to libCGAL once it's stable.
bool
divides(const MP_Float & d, const MP_Float & n)
{
  return internal::division(n, d).second == 0;
}

} // namespace internal

inline
bool
is_integer(const MP_Float &m)
{
  return m.is_integer();
}



inline
MP_Float
operator%(const MP_Float& n1, const MP_Float& n2)
{
  return internal::division(n1, n2).second;
}


// The namespace INTERN_MP_FLOAT contains global functions like square or sqrt
// which collide with the global functor adapting functions provided by the new
// AST/RET concept.
//
// TODO: IMHO, a better solution would be to put the INTERN_MP_FLOAT-functions
//       into the MP_Float-class... But there is surely a reason why this is not
//       the case..?


namespace INTERN_MP_FLOAT {
  inline
  MP_Float
  div(const MP_Float& n1, const MP_Float& n2)
  {
    return internal::division(n1, n2).first;
  }

  inline
  MP_Float
  gcd( const MP_Float& a, const MP_Float& b)
  {
    if (a == 0) {
      if (b == 0)
        return 0;
      MP_Float tmp=b;
      tmp.gcd_normalize();
      return tmp;
    }
    if (b == 0) {
      MP_Float tmp=a;
      tmp.gcd_normalize();
      return tmp;
    }

    MP_Float x = a, y = b;
    while (true) {
      x = x % y;
      if (x == 0) {
        CGAL_postcondition(internal::divides(y, a) & internal::divides(y, b));
        y.gcd_normalize();
        return y;
      }
      swap(x, y);
    }
  }

} // INTERN_MP_FLOAT


inline
void
simplify_quotient(MP_Float & numerator, MP_Float & denominator)
{
  // Currently only simplifies the two exponents.
#if 0
  // This better version causes problems as the I/O is changed for
  // Quotient<MP_Float>, which then does not appear as rational 123/345,
  // 1.23/3.45, this causes problems in the T2 test-suite (to be investigated).
  numerator.exp -= denominator.exp
                    + (MP_Float::exponent_type) denominator.v.size();
  denominator.exp = - (MP_Float::exponent_type) denominator.v.size();
#elif 1
  numerator.exp -= denominator.exp;
  denominator.exp = 0;
#else
  if (numerator != 0 && denominator != 0) {
    numerator.exp -= denominator.exp;
    denominator.exp = 0;
    const MP_Float g = gcd(numerator, denominator);
    numerator = integral_division(numerator, g);
    denominator = integral_division(denominator, g);
  }
  numerator.exp -= denominator.exp;
  denominator.exp = 0;
#endif
}

inline void simplify_root_of_2(MP_Float &/*a*/, MP_Float &/*b*/, MP_Float&/*c*/) {
#if 0
  if(is_zero(a)) {
  	simplify_quotient(b,c); return;
  } else if(is_zero(b)) {
  	simplify_quotient(a,c); return;
  } else if(is_zero(c)) {
  	simplify_quotient(a,b); return;
  }
  MP_Float::exponent_type va = a.exp +
    (MP_Float::exponent_type) a.v.size();
  MP_Float::exponent_type vb = b.exp +
    (MP_Float::exponent_type) b.v.size();
  MP_Float::exponent_type vc = c.exp +
    (MP_Float::exponent_type) c.v.size();
  MP_Float::exponent_type min = (std::min)((std::min)(va,vb),vc);
  MP_Float::exponent_type max = (std::max)((std::max)(va,vb),vc);
  MP_Float::exponent_type med = (min+max)/2.0;
  a.exp -= med;
  b.exp -= med;
  c.exp -= med;
#endif
}

namespace internal {
  inline void simplify_3_exp(int &a, int &b, int &c) {
    int min = (std::min)((std::min)(a,b),c);
    int max = (std::max)((std::max)(a,b),c);
    int med = (min+max)/2;
    a -= med;
    b -= med;
    c -= med;
  }
}


// specialization of to double functor
template<>
class Real_embeddable_traits< Quotient<MP_Float> >
    : public INTERN_QUOTIENT::Real_embeddable_traits_quotient_base<
Quotient<MP_Float> >{
public:
    struct To_double: public std::unary_function<Quotient<MP_Float>, double>{
         inline
         double operator()(const Quotient<MP_Float>& q) const {
            return INTERN_MP_FLOAT::to_double(q);
        }
    };
    struct To_interval
        : public std::unary_function<Quotient<MP_Float>, std::pair<double,double> > {
        inline
        std::pair<double,double> operator()(const Quotient<MP_Float>& q) const {
            return INTERN_MP_FLOAT::to_interval(q);
        }
    };
};

inline MP_Float min BOOST_PREVENT_MACRO_SUBSTITUTION(const MP_Float& x,const MP_Float& y){
  return (x<=y)?x:y; 
}
inline MP_Float max BOOST_PREVENT_MACRO_SUBSTITUTION(const MP_Float& x,const MP_Float& y){
  return (x>=y)?x:y; 
}


// Coercion_traits
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(MP_Float)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int, MP_Float)


} //namespace CGAL

namespace Eigen {
  template<class> struct NumTraits;
  template<> struct NumTraits<CGAL::MP_Float>
  {
    typedef CGAL::MP_Float Real;
    typedef CGAL::Quotient<CGAL::MP_Float> NonInteger;
    typedef CGAL::MP_Float Nested;
    typedef CGAL::MP_Float Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 1, // Is this lie right?
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 40,
      MulCost = 40
    };
  };
}

#include <CGAL/MP_Float_impl.h>

//specialization for Get_arithmetic_kernel
#include <CGAL/MP_Float_arithmetic_kernel.h>

#endif // CGAL_MP_FLOAT_H
