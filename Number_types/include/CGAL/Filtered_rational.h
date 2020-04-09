// Copyright (c) 2019-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri,
//                 Mael Rouxel-Labbé

#ifndef CGAL_FILTERED_RATIONAL_H
#define CGAL_FILTERED_RATIONAL_H

#include <CGAL/Gmpq.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/number_type_basic.h>

#include <boost/mpl/or.hpp>
#include <boost/utility/enable_if.hpp>

#include <iostream>
#include <sstream>
#include <type_traits>
#include <utility>

// #define CGAL_NT_CHECK_DEBUG(s) std::cerr << s << std::endl
#define CGAL_NT_CHECK_DEBUG(s)

namespace CGAL {

// A number type class, parameterized by 2 number types NT1 and NT2.
// It runs all operations on parallel over NT1 and NT2.
// It is also parameterized by a comparator which compares the values
// of NT1 and NT2 after all arithmetic operations.
template <typename NT1 = Interval_nt<false>,
          typename NT2 = Gmpq>
class Filtered_rational
{
  typedef Filtered_rational<NT1, NT2>                   Self;

  // NT1 must be CGAL::Interval_nt, but the protection (the template) is up to the user
  CGAL_static_assertion((std::is_same<NT1, CGAL::Interval_nt<true> >::value) ||
                        (std::is_same<NT1, CGAL::Interval_nt<false> >::value));

public:
  Filtered_rational() {}

  Filtered_rational(const int i)
    : _n1(i)
#ifndef CGAL_LAZY_FILTERED_RATIONAL_KERNEL // note that this is IF_NOT_DEF
    , _n2(i)
#endif
  {
    CGAL_assertion(is_valid());
  }

  Filtered_rational(const double d)
    : _n1(d)
#ifndef CGAL_LAZY_FILTERED_RATIONAL_KERNEL // note that this is IF_NOT_DEF
    , _n2(d)
#endif
  {
    CGAL_assertion(is_valid());
  }

  Filtered_rational(const NT1& n1, const NT2& n2)
    : _n1(n1), _n2(n2)
#ifdef CGAL_LAZY_FILTERED_RATIONAL_KERNEL
    , eii(true)
#endif
  {
    CGAL_assertion(is_valid());
  }

  Filtered_rational(const std::pair<NT1, NT2>& np)
    : _n1(np.first), _n2(np.second)
#ifdef CGAL_LAZY_FILTERED_RATIONAL_KERNEL
    , eii(true)
#endif
  {
    CGAL_assertion(is_valid());
  }

  Filtered_rational(const NT1& n1)
    : _n1(n1)
  {
#ifndef CGAL_LAZY_FILTERED_RATIONAL_KERNEL // note that this is IF_NOT_DEF
    CGAL_assertion(_n1.is_point());
    _n2 = NT2(_n1.inf());
#endif
  }

  // disabled if NT2 == NT1, int, or double
  template <typename T> // the template is only for SFINAE, this must be NT2
  Filtered_rational(const T& n2,
                    typename boost::disable_if<
                               boost::mpl::or_<
                                 boost::is_same<T, NT1>,
                                 boost::is_same<T, int>,
                                 boost::is_same<T, double> > >::type* = 0)
    : _n1(to_interval(_n2)), _n2(n2)
#ifdef CGAL_LAZY_FILTERED_RATIONAL_KERNEL
    , eii(true)
#endif
  {
    CGAL_static_assertion((std::is_same<T, NT2>::value));
    CGAL_assertion(is_valid());
  }

  Self& operator=(const std::pair<NT1, NT2>& np)
  {
    _n1 = np.first;
    _n2 = np.second;
#ifdef CGAL_LAZY_FILTERED_RATIONAL_KERNEL
    eii = true;
#endif
  }

  Self& operator+=(const Self& a) { n2() += a.n2(); n1() = to_interval(n2()); return *this; }
  Self& operator-=(const Self& a) { n2() -= a.n2(); n1() = to_interval(n2()); return *this; }
  Self& operator*=(const Self& a) { n2() *= a.n2(); n1() = to_interval(n2()); return *this; }
  Self& operator/=(const Self& a) { n2() /= a.n2(); n1() = to_interval(n2()); return *this; }

  // Accessors and setters.
  const NT1& approx() const { return _n1; }
  const NT1& n1() const { return _n1; }
  NT1& n1() { return _n1; }

  const NT2& exact() const { return n2(); }
  const NT2& n2() const
  {
#ifdef CGAL_LAZY_FILTERED_RATIONAL_KERNEL
    if(!eii)
    {
      CGAL_assertion(_n1.is_point());
      _n2 = NT2(_n1.inf());
      eii = true;
    }
#endif
    return _n2;
  }

  NT2& n2()
  {
#ifdef CGAL_LAZY_FILTERED_RATIONAL_KERNEL
    if(!eii)
    {
      CGAL_assertion(_n1.is_point());
      _n2 = NT2(_n1.inf());
      eii = true;
    }
#endif
    return _n2;
  }

  bool is_valid() const
  {
    if(!CGAL::is_valid(_n1))
      return false;

#ifdef CGAL_LAZY_FILTERED_RATIONAL_KERNEL
    if(!eii)
      return true; // _n2 is not yet built, so nothing more can be checked
#endif

    if(!CGAL::is_valid(_n2))
      return false;

    // the exact number must be in the interval
    const NT1 n2i = to_interval(n2());
    return (n2i.inf() >= n1().inf()) && (n2i.sup() <= n1().sup());
  }

private:
  NT1 _n1;
  mutable NT2 _n2;
#ifdef CGAL_LAZY_FILTERED_RATIONAL_KERNEL
  mutable bool eii = false; // eii means "exact is initialized"
#endif
};

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator+(const Filtered_rational<NT1, NT2>& a)
{
  CGAL_NT_CHECK_DEBUG("operator+");
  return Filtered_rational<NT1, NT2>(+ a.n1(), + a.n2() );
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator+(const Filtered_rational<NT1, NT2>& a,
          const Filtered_rational<NT1, NT2>& b)
{
  CGAL_NT_CHECK_DEBUG("operator+");
  return Filtered_rational<NT1, NT2>(a.n2() + b.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator+(int a, const Filtered_rational<NT1, NT2>& b)
{
  CGAL_NT_CHECK_DEBUG("operator+");
  return Filtered_rational<NT1, NT2>(a + b.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator+(const Filtered_rational<NT1, NT2>& a, int b)
{
  CGAL_NT_CHECK_DEBUG("operator+");
  return Filtered_rational<NT1, NT2>(a.n2() + b);
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator-(const Filtered_rational<NT1, NT2>& a,
          const Filtered_rational<NT1, NT2>& b)
{
  CGAL_NT_CHECK_DEBUG("operator-");
  return Filtered_rational<NT1, NT2>(a.n2() - b.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator-(int a, const Filtered_rational<NT1, NT2>& b)
{
  CGAL_NT_CHECK_DEBUG("operator-");
  return Filtered_rational<NT1, NT2>(a - b.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator-(const Filtered_rational<NT1, NT2>& a, int b)
{
  CGAL_NT_CHECK_DEBUG("operator-");
  return Filtered_rational<NT1, NT2>(a.n2() - b);
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator-(const Filtered_rational<NT1, NT2>& a)
{
  CGAL_NT_CHECK_DEBUG("unary operator-");
  return Filtered_rational<NT1, NT2>(-a.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator*(const Filtered_rational<NT1, NT2>& a,
          const Filtered_rational<NT1, NT2>& b)
{
  CGAL_NT_CHECK_DEBUG("operator*");
  return Filtered_rational<NT1, NT2>(a.n2() * b.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator*(int a, const Filtered_rational<NT1, NT2>& b)
{
  CGAL_NT_CHECK_DEBUG("operator*");
  return Filtered_rational<NT1, NT2>(a * b.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator*(const Filtered_rational<NT1, NT2>& a, int b)
{
  CGAL_NT_CHECK_DEBUG("operator*");
  return Filtered_rational<NT1, NT2>(a.n2() * b);
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator/(const Filtered_rational<NT1, NT2>& a,
          const Filtered_rational<NT1, NT2>& b)
{
  CGAL_NT_CHECK_DEBUG("operator/");
  return Filtered_rational<NT1, NT2>(a.n2() / b.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator/(int a, const Filtered_rational<NT1, NT2>& b)
{
  CGAL_NT_CHECK_DEBUG("operator/");
  return Filtered_rational<NT1, NT2>(a / b.n2());
}

template <typename NT1, typename NT2>
Filtered_rational<NT1, NT2>
operator/(const Filtered_rational<NT1, NT2>& a, int b)
{
  CGAL_NT_CHECK_DEBUG("operator/");
  return Filtered_rational<NT1, NT2>(a.n2() / b);
}

// arithmetic operators end
// compare operators begin

template <typename NT1, typename NT2>
bool
operator==(const Filtered_rational<NT1, NT2>& a,
           const Filtered_rational<NT1, NT2>& b)
{
  if(a.n1().is_point() && (b.n1().is_point()))
    return a.n1().inf() == b.n1().inf();

  if(a.n1().do_overlap(b.n1()))
    return a.n2() == b.n2();

  return false;
}

template <typename NT1, typename NT2>
bool
operator!=(const Filtered_rational<NT1, NT2>& a,
           const Filtered_rational<NT1, NT2>& b)
{
  return ! (a == b);
}

template <typename NT1, typename NT2>
bool
operator<(const Filtered_rational<NT1, NT2>& a,
          const Filtered_rational<NT1, NT2>& b)
{
  if(a.n1().sup() < b.n1().inf())
    return true;

  if(a.n1().inf() > b.n1().sup())
    return false;

  return a.n2() < b.n2();
}

template <typename NT1, typename NT2>
bool
operator>(const Filtered_rational<NT1, NT2>& a,
          const Filtered_rational<NT1, NT2>& b)
{
  return b < a;
}

template <typename NT1, typename NT2>
bool
operator<=(const Filtered_rational<NT1, NT2>& a,
           const Filtered_rational<NT1, NT2>& b)
{
  if(a.n1().sup() <= b.n1().inf())
    return true;

  if(a.n1().inf() >= b.n1().sup())
    return false;

  return a.n2() <= b.n2();
}

template <typename NT1, typename NT2>
bool
operator>=(const Filtered_rational<NT1, NT2>& a,
           const Filtered_rational<NT1, NT2>& b)
{
  return b <= a;
}

template <typename NT1, typename NT2>
bool
operator==(const Filtered_rational<NT1, NT2>& a, int i)
{
  bool b1 = a.n1() == i;
  CGAL_assertion(b1 == ( a.n2() == i ) );
  return b1;
}

template <typename NT1, typename NT2>
bool
operator!=(const Filtered_rational<NT1, NT2>& a, int i)
{
  Uncertain<bool> ub = a.n1() != i;
  if(is_certain(ub))
    return get_certain(ub);

  return a.n2() != i;
}

template <typename NT1, typename NT2>
bool
operator<(const Filtered_rational<NT1, NT2>& a, int i)
{
  Uncertain<bool> ub = a.n1() < i;
  if(is_certain(ub))
    return get_certain(ub);

  return a.n2() < i;
}

template <typename NT1, typename NT2>
bool
operator>(const Filtered_rational<NT1, NT2>& a, int i)
{
  Uncertain<bool> ub = a.n1() > i;
  if(is_certain(ub))
    return get_certain(ub);

  return a.n2() > i;
}

template <typename NT1, typename NT2>
bool
operator<=(const Filtered_rational<NT1, NT2>& a, int i)
{
  Uncertain<bool> ub = a.n1() <= i;
  if(is_certain(ub))
    return get_certain(ub);

  return a.n2() <= i;
}

template <typename NT1, typename NT2>
bool
operator>=(const Filtered_rational<NT1, NT2>& a, int i)
{
  Uncertain<bool> ub = a.n1() >= i;
  if(is_certain(ub))
    return get_certain(ub);

  return a.n2() >= i;
}

template <typename NT1, typename NT2>
bool
operator==(int i, const Filtered_rational<NT1, NT2>& b)
{
  return b == i;
}

template <typename NT1, typename NT2>
bool
operator!=(int i, const Filtered_rational<NT1, NT2>& b)
{
  return b != i;
}

template <typename NT1, typename NT2>
bool
operator<(int i, const Filtered_rational<NT1, NT2>& b)
{
  return b > i;
}

template <typename NT1, typename NT2>
bool
operator>(int i, const Filtered_rational<NT1, NT2>& b)
{
  return b < i;
}

template <typename NT1, typename NT2>
bool
operator<=(int i, const Filtered_rational<NT1, NT2>& b)
{
  return b >= i;
}

template <typename NT1, typename NT2>
bool
operator>=(int i, const Filtered_rational<NT1, NT2>& b)
{
  return b <= i;
}

// compare operators end
// Functors Is_valid

template <typename NT1, typename NT2>
class Is_valid< Filtered_rational<NT1, NT2> >
    : public CGAL::cpp98::unary_function< Filtered_rational<NT1, NT2> , bool > {
public :
  bool operator()(const  Filtered_rational<NT1, NT2>& a ) const {
    bool b1 = is_valid(a.n1());
    CGAL_assertion(b1 == is_valid(a.n2()) );
    return b1;
  }
};

namespace FR_INTERN{
// -----------------------------

// fwd
template < typename Filtered_rational, typename Algebraic_category>
class FR_AST_base
  : public Algebraic_structure_traits_base< Filtered_rational , Null_tag>
{ };

template <typename NT1, typename NT2>
class FR_AST_base<Filtered_rational<NT1, NT2>, Integral_domain_without_division_tag>
  : public Algebraic_structure_traits_base<Filtered_rational<NT1, NT2>,
                                           Integral_domain_without_division_tag>
{
private:
  typedef Algebraic_structure_traits<NT1>   AST1;
  typedef Algebraic_structure_traits<NT2>   AST2;
  typedef Filtered_rational<NT1, NT2>       Type;

public:
  //CGAL::Algebraic_structure_traits<>::Simplify
  class Simplify
    : public CGAL::cpp98::unary_function< Type& , void >
  {
  public:
    void operator()( Type& a) const
    {
      typename AST1::Simplify()(a.n1());
      typename AST2::Simplify()(a.n2());
      CGAL_assertion(a.is_valid());
    }
  };

  //CGAL::Algebraic_structure_traits< >::Is_zero
  class Is_zero
    : public CGAL::cpp98::unary_function< Type, bool >
  {
  public:
    bool operator()(const  Type& a ) const
    {
      Uncertain<bool> ub = typename AST1::Is_zero()(a.n1());
      if(is_certain(ub)){
        return get_certain(ub);
      }
      return typename AST2::Is_zero()(a.n2());
    }
  };

  // CGAL::Algebraic_structure_traits< >::Is_one
  class Is_one
    : public CGAL::cpp98::unary_function< Type, bool >
  {
  public:
    bool operator()(const Type& a) const
    {
      Uncertain<bool> ub = typename AST1::Is_one()(a.n1());
      if(is_certain(ub))
        return get_certain(ub);

      return typename AST2::Is_one()(a.n2());
    }
  };

  // CGAL::Algebraic_structure_traits<  >::Square
  class Square
    : public CGAL::cpp98::unary_function< Type , Type >
  {
  public:
    Type operator()(const Type& a) const
    {
      return Type(typename AST2::Square()(a.n2()));
    }
  };


  // CGAL::Algebraic_structure_traits<  >::Unit_part
  class Unit_part
    : public CGAL::cpp98::unary_function< Type , Type >
  {
  public:
    Type operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Unit_part");
      return Type(typename AST2::Unit_part()(a.n2()));
    }
  };
};

template <typename NT1, typename NT2>
class FR_AST_base<Filtered_rational<NT1, NT2> , Integral_domain_tag>
  : public FR_AST_base<Filtered_rational<NT1, NT2>,
                       Integral_domain_without_division_tag>
{
private:
  typedef Algebraic_structure_traits<NT1> AST1;
  typedef Algebraic_structure_traits<NT2> AST2;
  typedef Filtered_rational<NT1, NT2> Type;

public:
  // CGAL::Algebraic_structure_traits< >::Integral_division
  class Integral_division
    : public CGAL::cpp98::binary_function< Type, Type, Type >
  {
  public:
    Type operator()( const Type& a, const Type& b) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Integral_division");
      return Type(typename AST2::Integral_division()(a.n2(),b.n2()));
    }
  };

  class Divides
    : public CGAL::cpp98::binary_function< Type, Type, bool >
  {
  public:
    bool operator()( const Type& a, const Type& b) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Divides");
      Uncertain<bool> ub =  typename AST1::Divides()(a.n1(),b.n1());
      if(is_certain(ub))
        return get_certain(ub);

      return typename AST2::Divides()(a.n2(),b.n2());
    }

    bool operator()( const Type& a, const Type& b, Type& q) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Divides");
      NT1 q1;
      bool result1 =  typename AST1::Divides()(a.n1(), b.n1(), q1);
      NT2 q2;
      CGAL_assertion_code( bool result2 = ) // needed for CGAL_assert only
      typename AST2::Divides()(a.n2(), b.n2(), q2);
      q = Type(q1,q2);
      CGAL_assertion(result1 == result2);
      return result1;
    }

    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(Type, bool)
  };
};


template <typename NT1, typename NT2>
class FR_AST_base<Filtered_rational<NT1, NT2>,
                  Unique_factorization_domain_tag>
  : public FR_AST_base<Filtered_rational<NT1, NT2>, Integral_domain_tag >
{
private:
  typedef Algebraic_structure_traits<NT1> AST1;
  typedef Algebraic_structure_traits<NT2> AST2;
  typedef Filtered_rational<NT1, NT2> Type;
public:
  // CGAL::Algebraic_structure_traits< >::Gcd
  class Gcd
    : public CGAL::cpp98::binary_function< Type, Type, Type >
  {
  public:
    Type operator()(const Type& a,
                    const Type& b) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Gcd");
      return Type(typename AST1::Gcd()(a.n1(),b.n1()),
                  typename AST2::Gcd()(a.n2(),b.n2()));
    }
  };
};

template <typename NT1, typename NT2>
class FR_AST_base<Filtered_rational<NT1, NT2>, Euclidean_ring_tag>
  : public FR_AST_base<Filtered_rational<NT1, NT2> , Unique_factorization_domain_tag >
{
private:
  typedef Algebraic_structure_traits<NT1> AST1;
  typedef Algebraic_structure_traits<NT2> AST2;
  typedef Filtered_rational<NT1, NT2> Type;

public:
  // CGAL::Algebraic_structure_traits< >::Div
  class Div
    : public CGAL::cpp98::binary_function< Type, Type, Type >
  {
  public:
    Type operator()(const Type& a,
                    const Type& b) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Div");
      return Type(typename AST1::Div()(a.n1(),b.n1()),
                  typename AST2::Div()(a.n2(),b.n2()));
    }
  };

  // CGAL::Algebraic_structure_traits< >::Mod
  class Mod
    : public CGAL::cpp98::binary_function< Type, Type, Type >
  {
  public:
    Type operator()(const Type& a,
                    const Type& b) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Mod");
      return Type(typename AST1::Mod()(a.n1(),b.n1()),
                  typename AST2::Mod()(a.n2(),b.n2()));
    }
  };

  // CGAL::Algebraic_structure_traits< >::Div_mod
  class Div_mod
  {
  public:
    typedef Type    first_argument_type;
    typedef Type    second_argument_type;
    typedef Type&   third_argument_type;
    typedef Type&   fourth_argument_type;
    typedef void  result_type;

    void operator()(const Type& a,
                    const Type& b,
                    Type& q,
                    Type& r) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Div_mod");
      NT1 q1,r1;
      NT2 q2,r2;
      typename AST1::Div_mod()(a.n1(),b.n1(),q1,r1);
      typename AST2::Div_mod()(a.n2(),b.n2(),q2,r2);
      q = Type(q1,q2);
      r = Type(r1,r2);
    }
  };
};

template <typename NT1, typename NT2>
class FR_AST_base<Filtered_rational<NT1, NT2> , Field_tag >
  : public FR_AST_base< Filtered_rational<NT1, NT2>, Integral_domain_tag>
{
private:
  typedef Algebraic_structure_traits<NT1> AST1;
  typedef Algebraic_structure_traits<NT2> AST2;
  typedef Filtered_rational<NT1, NT2> Type;

public:
  class Inverse
    : public CGAL::cpp98::unary_function< Type, Type >
  {
  public:
    Type operator()( const Type& a ) const {
      NT1 r1 = typename AST1::Inverse()(a.n1());
      NT2 r2 = typename AST2::Inverse()(a.n2());
      return Type(r1,r2);
    }
  };
};

template <typename NT1, typename NT2>
class FR_AST_base< Filtered_rational<NT1, NT2> , Field_with_sqrt_tag >
  : public FR_AST_base< Filtered_rational<NT1, NT2> , Field_tag >
{
private:
  typedef Algebraic_structure_traits<NT1> AST1;
  typedef Algebraic_structure_traits<NT2> AST2;
  typedef Filtered_rational<NT1, NT2> Type;

public:
  // CGAL::Algebraic_structure_traits<  >::Sqrt
  class Sqrt
    : public CGAL::cpp98::unary_function< Type , Type >
  {
  public:
    Type operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("AST::Sqrt");
      return Type(typename AST1::Sqrt()(a.n1()),
                  typename AST2::Sqrt()(a.n2()));
    }
  };
};

} // namespace FR_INTERN

template <typename NT1, typename NT2>
class Algebraic_structure_traits <Filtered_rational<NT1, NT2> >
  : public FR_INTERN::FR_AST_base<Filtered_rational<NT1, NT2> ,
                                  typename Algebraic_structure_traits<NT2>::Algebraic_category >
{
  typedef Algebraic_structure_traits<NT2> AST1;

public:
  typedef Filtered_rational<NT1, NT2> Type;
  typedef typename AST1::Algebraic_category Algebraic_category;
  typedef typename AST1::Is_exact Is_exact;
};

namespace FR_INTERN {

template < typename Filtered_rational, typename Is_real_embeddable >
class FR_RET_base;

template < typename NT >
class FR_RET_base<NT,Tag_false> : public Real_embeddable_traits<NT>
{ };

template <typename NT1, typename NT2>
class FR_RET_base<Filtered_rational<NT1, NT2> , Tag_true>
  : public INTERN_RET::Real_embeddable_traits_base< Filtered_rational<NT1, NT2> , CGAL::Tag_true >
{
private:
  typedef Real_embeddable_traits<NT1> RET1;
  typedef Real_embeddable_traits<NT2> RET2;
  typedef Filtered_rational<NT1, NT2> Type;

public:
  // CGAL::Real_embeddable_traits<  >::Abs
  class Abs
    : public CGAL::cpp98::unary_function< Type , Type >
  {
  public:
    Type operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("RET::Abs");
      return Type(typename RET1::Abs()(a.n1()),
                  typename RET2::Abs()(a.n2()));
    }
  };

  // CGAL::Real_embeddable_traits<  >::Sign
  class Sgn
    : public CGAL::cpp98::unary_function< Type , ::CGAL::Sign >
  {
  public:
    ::CGAL::Sign operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("RET::Sign");
      Uncertain<Sign> ub =  typename RET1::Sgn()(a.n1());
      if(is_certain(ub))
        return get_certain(ub);

      return typename RET2::Sgn()(a.n2() );
    }
  };

  // CGAL::Real_embeddable_traits<  >::Is_finite
  class Is_finite
    : public CGAL::cpp98::unary_function< Type , bool >
  {
  public:
    bool operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("RET::Is_finite");
      bool r1 =  typename RET1::Is_finite()(a.n1());
      // CGAL_assertion( r1 == typename RET2::Is_finite()(a.n2()) );
      // AF: what to do??
      // See what Gmpq does
      return r1;
    }
  };

  // CGAL::Real_embeddable_traits<  >::Is_positive
  class Is_positive
    : public CGAL::cpp98::unary_function< Type , bool >
  {
  public:
    bool operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("RET::Is_positive");
      Uncertain<bool> ub =  typename RET1::Is_positive()(a.n1());
      if(is_certain(ub))
        return get_certain(ub);

      return typename RET2::Is_positive()(a.n2());
    }
  };

  // CGAL::Real_embeddable_traits<  >::Is_negative
  class Is_negative
    : public CGAL::cpp98::unary_function< Type , bool >
  {
  public:
    bool operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("RET::Is_negative");
      bool r1 =  typename RET1::Is_negative()(a.n1());
      CGAL_assertion( r1 == typename RET2::Is_negative()(a.n2()) );
      return r1;
    }
  };

  // CGAL::Real_embeddable_traits<  >::Is_zero
  class Is_zero
    : public CGAL::cpp98::unary_function< Type , bool >
  {
  public:
    bool operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("RET::Is_zero");
      Uncertain<bool> ub =  typename RET1::Is_zero()(a.n1());
      if(is_certain(ub))
        return get_certain(ub);

      return typename RET2::Is_zero()(a.n2() );
    }
  };

  // CGAL::Real_embeddable_traits<  >::Compare
  class Compare
    : public CGAL::cpp98::binary_function< Type , Type, Comparison_result >
  {
  public:
    Comparison_result operator()(const Type& a, const Type& b) const
    {
      CGAL_NT_CHECK_DEBUG("RET::Compare");
      Uncertain<Comparison_result> ucr =  typename RET1::Compare()(a.n1(),b.n1());
      if(is_certain(ucr))
        return get_certain(ucr);

      return typename RET2::Compare()(a.n2(),b.n2());
    }
  };

  // CGAL::Real_embeddable_traits<  >::To_double
  class To_double
    : public CGAL::cpp98::unary_function< Type , double >
  {
  public:
    double operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("RET::To_double");
      double r1 =  typename RET1::To_double()(a.n1());
      // AF: Do we have to check something ??
      // YES See Lazy_exact_nt
      // r1 == typename RET2::To_double()(a.n2()) );
      return r1;
    }
  };

  // CGAL::Real_embeddable_traits<  >::To_interval
  class To_interval
    : public CGAL::cpp98::unary_function< Type , std::pair<double, double> >
  {
  public:
    std::pair<double, double> operator()(const Type& a) const
    {
      CGAL_NT_CHECK_DEBUG("RET::To_interval");
      std::pair<double, double> r1 =  typename RET1::To_interval()(a.n1());
      return r1;
    }
  };
};

} // namespace FR_INTERN

template <typename NT1, typename NT2>
class Real_embeddable_traits< Filtered_rational<NT1, NT2> >
  : public FR_INTERN::FR_RET_base<Filtered_rational<NT1, NT2>,
                                  typename Real_embeddable_traits<NT1>::Is_real_embeddable >
{
  typedef Real_embeddable_traits<NT1> RET1;

public:
  typedef Filtered_rational<NT1, NT2>       Type;
  typedef typename RET1::Is_real_embeddable Is_real_embeddable;
};

template <typename NT1, typename NT2>
struct Coercion_traits< Filtered_rational<NT1, NT2>, Filtered_rational<NT1, NT2> >
{
  typedef Tag_true  Are_explicit_interoperable;
  typedef Tag_true  Are_implicit_interoperable;
  typedef Filtered_rational<NT1, NT2> Type;

  struct Cast
  {
    typedef Type result_type;
    Type operator()(const Type& x) const { return x;}
  };
};

template <typename NT1, typename NT2>
struct Coercion_traits< Filtered_rational<NT1, NT2>, int >
{
  typedef Tag_true  Are_explicit_interoperable;
  typedef Tag_true  Are_implicit_interoperable;
  typedef Filtered_rational<NT1, NT2> Type;

  struct Cast
  {
    typedef Type result_type;
    Type operator()(const Type& x) const { return x;}
    Type operator()(int x) const { return Type(x);}
  };
};

template <typename NT1, typename NT2>
struct Coercion_traits< int , Filtered_rational<NT1, NT2> >
  : public Coercion_traits< Filtered_rational<NT1, NT2> , int >
{ };

namespace FR_INTERN {

template < typename NT_checker, typename Tag = Tag_false >
struct Coercion_traits_double
{
  typedef Tag_false  Are_explicit_interoperable;
  typedef Tag_false  Are_implicit_interoperable;
  typedef Null_tag Type;
};

template <typename NT1, typename NT2>
struct Coercion_traits_double< Filtered_rational<NT1, NT2> , Tag_true >
{
  typedef Tag_true  Are_explicit_interoperable;
  typedef Tag_true  Are_implicit_interoperable;
  typedef Filtered_rational<NT1, NT2> Type;

  struct Cast
  {
    typedef Type result_type;
    Type operator()(const Type& x) const { return x;}
    Type operator()(const double& x) const {
      return Type(x);
    }
  };
};

} // namespace FR_INTERN

template <typename NT1, typename NT2>
struct Coercion_traits< Filtered_rational<NT1, NT2>, double >
  : public FR_INTERN::Coercion_traits_double<Filtered_rational<NT1, NT2>,
                                             typename Coercion_traits<NT1, double>::Are_implicit_interoperable>
{ };

template <typename NT1, typename NT2>
struct Coercion_traits< double , Filtered_rational<NT1, NT2> >
  : public Coercion_traits< Filtered_rational<NT1, NT2>, double >
{ };

// What to do with the IO ?
// - either pick one as the main, and convert to the second.
// - output/input pairs, and read both (=> problems with FP values).

template <typename NT1, typename NT2>
std::ostream&
operator<< (std::ostream& os, const Filtered_rational<NT1, NT2>& b)
{
  return os << to_double(b.n1());
}

template <typename NT1, typename NT2>
std::istream&
operator>> (std::istream& is, Filtered_rational<NT1, NT2>& b)
{
  is >> b.n2();
  b.n1() = to_interval(b.n2());
  return is;
}

template <typename NT1, typename NT2>
bool fit_in_double(const Filtered_rational<NT1, NT2>& a, double& r)
{
  return fit_in_double(a.n1(),r);
}

} // namespace CGAL

#endif // CGAL_FILTERED_RATIONAL_H
