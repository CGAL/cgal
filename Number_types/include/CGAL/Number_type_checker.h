// Copyright (c) 2005-2007
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
// Author(s)     : Sylvain Pion, Michael Hemmer

#ifndef CGAL_NUMBER_TYPE_CHECKER_H
#define CGAL_NUMBER_TYPE_CHECKER_H

#include <CGAL/number_type_basic.h>
#include <sstream>

// A number type class, parameterized by 2 number types NT1 and NT2.
// It runs all operations on parallel over NT1 and NT2.
// It is also parameterized by a comparator which compares the values
// of NT1 and NT2 after all arithmetic operations.

// #define CGAL_NT_CHECK_DEBUG(s) std::cerr << s << std::endl
#define CGAL_NT_CHECK_DEBUG(s)

namespace CGAL {

namespace internal {

// std::equal_to<T> sticks things to one type, so I define my own.
struct Equal_to
{
  template < typename NT1, typename NT2 >
  bool operator()(const NT1& a, const NT2& b) const
  { return a == b; }
};

} // namespace internal


template < typename NT1, typename NT2, typename Cmp = internal::Equal_to >
class Number_type_checker
{
  NT1 _n1;
  NT2 _n2;
  Cmp cmp;
public:

  typedef Number_type_checker<NT1, NT2, Cmp> Self;

  Number_type_checker() {}
  Number_type_checker(int i)
    : _n1(i), _n2(i) { CGAL_assertion(is_valid()); }
  Number_type_checker(double d)
    : _n1(d), _n2(d) { CGAL_assertion(is_valid()); }
  Number_type_checker(const NT1 &n1, const NT2 &n2)
    : _n1(n1), _n2(n2) { CGAL_assertion(is_valid()); }

  // The following need to be dependant on NT1 != {NT2,int,double} ...
  //Number_type_checker(const NT1 &n1) : _n1(n1), _n2(n1) {}
  //Number_type_checker(const NT2 &n2) : _n1(n2), _n2(n2) {}

  Self& operator+=(const Self &a)
  { n1() += a.n1(); n2() += a.n2(); CGAL_assertion(is_valid()); return *this; }

  Self& operator-=(const Self &a)
  { n1() -= a.n1(); n2() -= a.n2(); CGAL_assertion(is_valid()); return *this; }

  Self& operator*=(const Self &a)
  { n1() *= a.n1(); n2() *= a.n2(); CGAL_assertion(is_valid()); return *this; }

  Self& operator/=(const Self &a)
  { n1() /= a.n1(); n2() /= a.n2(); CGAL_assertion(is_valid()); return *this; }

  // Accessors and setters.
  const NT1& n1() const { return _n1; }
  const NT2& n2() const { return _n2; }

  NT1& n1() { return _n1; }
  NT2& n2() { return _n2; }

  // Validity checking.
  bool is_valid() const
  {
    CGAL_NT_CHECK_DEBUG("Checking...");
    bool b = cmp(_n1, _n2);
    if (!b)
    {
      CGAL_NT_CHECK_DEBUG("Different values :");
      CGAL_NT_CHECK_DEBUG("n1 = " << _n1);
      CGAL_NT_CHECK_DEBUG("n2 = " << _n2);
    }
    return b;
  }
};



template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator+(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   CGAL_NT_CHECK_DEBUG("operator+");
   return Number_type_checker<NT1, NT2, Cmp>(+ a.n1(), + a.n2() );
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator+(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator+");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() + b.n1(), a.n2() + b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator+(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator+");
   return Number_type_checker<NT1, NT2, Cmp>(a + b.n1(), a + b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator+(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   CGAL_NT_CHECK_DEBUG("operator+");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() + b, a.n2() + b);
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator-(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator-");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() - b.n1(), a.n2() - b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator-(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator-");
   return Number_type_checker<NT1, NT2, Cmp>(a - b.n1(), a - b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator-(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   CGAL_NT_CHECK_DEBUG("operator-");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() - b, a.n2() - b);
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator-(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   CGAL_NT_CHECK_DEBUG("unary operator-");
   return Number_type_checker<NT1, NT2, Cmp>(-a.n1(), -a.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator*(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator*");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() * b.n1(), a.n2() * b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator*(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator*");
   return Number_type_checker<NT1, NT2, Cmp>(a * b.n1(), a * b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator*(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   CGAL_NT_CHECK_DEBUG("operator*");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() * b, a.n2() * b);
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator/(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator/");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() / b.n1(), a.n2() / b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator/(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator/");
   return Number_type_checker<NT1, NT2, Cmp>(a / b.n1(), a / b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator/(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   CGAL_NT_CHECK_DEBUG("operator/");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() / b, a.n2() / b);
}

// arithmetic operators end
// compare operators begin

template < typename NT1, typename NT2, typename Cmp >
bool
operator==(const Number_type_checker<NT1, NT2, Cmp> &a,
           const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() == b.n1();
   CGAL_assertion(b1 == ( a.n2() == b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator!=(const Number_type_checker<NT1, NT2, Cmp> &a,
           const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() != b.n1();
   CGAL_assertion(b1 == ( a.n2() != b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() < b.n1();
   CGAL_assertion(b1 == ( a.n2() < b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() > b.n1();
   CGAL_assertion(b1 == ( a.n2() > b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<=(const Number_type_checker<NT1, NT2, Cmp> &a,
           const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() <= b.n1();
   CGAL_assertion(b1 == ( a.n2() <= b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>=(const Number_type_checker<NT1, NT2, Cmp> &a,
           const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() >= b.n1();
   CGAL_assertion(b1 == ( a.n2() >= b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator==(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() == i;
   CGAL_assertion(b1 == ( a.n2() == i ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator!=(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() != i;
   CGAL_assertion(b1 == ( a.n2() != i ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() < i;
   CGAL_assertion(b1 == ( a.n2() < i ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() > i;
   CGAL_assertion(b1 == ( a.n2() > i ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<=(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() <= i;
   CGAL_assertion(b1 == ( a.n2() <= i ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>=(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() >= i;
   CGAL_assertion(b1 == ( a.n2() >= i ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator==(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i == b.n1();
   CGAL_assertion(b1 == ( i == b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator!=(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i != b.n1();
   CGAL_assertion(b1 == ( i != b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i < b.n1();
   CGAL_assertion(b1 == ( i < b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i > b.n1();
   CGAL_assertion(b1 == ( i > b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<=(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i <= b.n1();
   CGAL_assertion(b1 == ( i <= b.n2() ) );
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>=(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i >= b.n1();
   CGAL_assertion(b1 == ( i >= b.n2() ) );
   return b1;
}

// compare operators end
// Functors Is_valid

template < typename NT1, typename NT2, typename Cmp >
class Is_valid< Number_type_checker<NT1, NT2, Cmp> >
    : public CGAL::cpp98::unary_function< Number_type_checker<NT1, NT2, Cmp> , bool > {
public :
    bool operator()(const  Number_type_checker<NT1, NT2, Cmp>& a ) const {
        bool b1 = is_valid(a.n1());
        CGAL_assertion(b1 == is_valid(a.n2()) );
        // Should we also call a.is_valid() ?
        return b1;
    }
};

namespace NTC_INTERN{
// -----------------------------

// fwd
template < typename Number_type_checker, typename Algebraic_category>
class NTC_AST_base
    :public Algebraic_structure_traits_base< Number_type_checker , Null_tag>{
};

template < typename NT1, typename NT2, typename Cmp >
class NTC_AST_base
< Number_type_checker<NT1, NT2, Cmp> , Integral_domain_without_division_tag>
:public Algebraic_structure_traits_base<Number_type_checker<NT1, NT2, Cmp>,
Integral_domain_without_division_tag>
{
private:
    typedef Algebraic_structure_traits<NT1> AST1;
    typedef Algebraic_structure_traits<NT2> AST2;
    typedef Number_type_checker<NT1, NT2, Cmp> Type;

public:
    //CGAL::Algebraic_structure_traits<>::Simplify
    class Simplify
        : public CGAL::cpp98::unary_function< Type& , void > {
    public:
        void operator()( Type& a) const {
            typename AST1::Simplify()(a.n1());
            typename AST2::Simplify()(a.n2());
            CGAL_assertion(a.is_valid());
        }
    };

    //CGAL::Algebraic_structure_traits< >::Is_zero
    class Is_zero
        : public CGAL::cpp98::unary_function< Type, bool > {
    public:
        bool operator()(const  Type& a ) const {
            bool b1 = typename AST1::Is_zero()(a.n1());
            CGAL_assertion(b1 == typename AST2::Is_zero()(a.n2()) );
            CGAL_assertion(a.is_valid());
            return b1;
        }
    };

    // CGAL::Algebraic_structure_traits< >::Is_one
    class Is_one
        : public CGAL::cpp98::unary_function< Type, bool > {
    public:
        bool operator()(const Type& a) const {
            bool b1 = typename AST1::Is_one()(a.n1());
            CGAL_assertion(b1 == typename AST2::Is_one()(a.n2()) );
            CGAL_assertion(a.is_valid());
            return b1;
        }
    };
    // CGAL::Algebraic_structure_traits<  >::Square
    class Square
        : public CGAL::cpp98::unary_function< Type , Type > {
    public:
        Type operator()(const Type& a) const {
            return Type(
                    typename AST1::Square()(a.n1()),
                    typename AST2::Square()(a.n2()));
        }
    };


    // CGAL::Algebraic_structure_traits<  >::Unit_part
    class Unit_part
        : public CGAL::cpp98::unary_function< Type , Type > {
    public:
        Type operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("AST::Unit_part");
            return Type(
                    typename AST1::Unit_part()(a.n1()),
                    typename AST2::Unit_part()(a.n2()));
        }
    };
};

template < typename NT1, typename NT2, typename Cmp >
class NTC_AST_base
      < Number_type_checker< NT1, NT2, Cmp> , Integral_domain_tag >
      :public  NTC_AST_base
      < Number_type_checker< NT1, NT2, Cmp> , Integral_domain_without_division_tag >
{
private:
    typedef Algebraic_structure_traits<NT1> AST1;
    typedef Algebraic_structure_traits<NT2> AST2;
    typedef Number_type_checker<NT1, NT2, Cmp> Type;

public:

// CGAL::Algebraic_structure_traits< >::Integral_division
    class Integral_division
        : public CGAL::cpp98::binary_function< Type, Type, Type > {
    public:
        Type operator()( const Type& a, const Type& b) const {
            CGAL_NT_CHECK_DEBUG("AST::Integral_division");
            return Type(
                    typename AST1::Integral_division()(a.n1(),b.n1()),
                    typename AST2::Integral_division()(a.n2(),b.n2()));
        }
    };

  class Divides
    : public CGAL::cpp98::binary_function< Type, Type, bool > {
  public:
    bool operator()( const Type& a, const Type& b) const {
      CGAL_NT_CHECK_DEBUG("AST::Divides");
      bool result =  typename AST1::Divides()(a.n1(),b.n1());
      CGAL_assertion(result == typename AST2::Divides()(a.n2(),b.n2()));
      return result;
    }
    bool operator()( const Type& a, const Type& b, Type& q) const {
      CGAL_NT_CHECK_DEBUG("AST::Divides");
      NT1 q1;
      bool result1 =  typename AST1::Divides()(a.n1(),b.n1(),q1);
      NT2 q2;
      CGAL_assertion_code( bool result2 = ) // needed for CGAL_assert only
        typename AST2::Divides()(a.n2(),b.n2(),q2);
      q = Type(q1,q2);
      CGAL_assertion(result1 == result2);
      return result1;
    }
    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(Type,bool)
  };
};


template < typename NT1, typename NT2, typename Cmp >
class NTC_AST_base
< Number_type_checker< NT1, NT2, Cmp> , Unique_factorization_domain_tag >
    :public  NTC_AST_base
< Number_type_checker< NT1, NT2, Cmp> , Integral_domain_tag >
{
private:
    typedef Algebraic_structure_traits<NT1> AST1;
    typedef Algebraic_structure_traits<NT2> AST2;
    typedef Number_type_checker<NT1, NT2, Cmp> Type;
public:
    // CGAL::Algebraic_structure_traits< >::Gcd
    class Gcd
        : public CGAL::cpp98::binary_function< Type,
                                  Type,
                                  Type > {
    public:
        Type operator()(
                const Type& a,
                const Type& b) const {
            CGAL_NT_CHECK_DEBUG("AST::Gcd");
            return Type(
                    typename AST1::Gcd()(a.n1(),b.n1()),
                    typename AST2::Gcd()(a.n2(),b.n2()));
        }
    };
};

template < typename NT1, typename NT2, typename Cmp >
class NTC_AST_base
< Number_type_checker< NT1, NT2, Cmp> , Euclidean_ring_tag >
    :public  NTC_AST_base
< Number_type_checker< NT1, NT2, Cmp> , Unique_factorization_domain_tag >
{
private:
    typedef Algebraic_structure_traits<NT1> AST1;
    typedef Algebraic_structure_traits<NT2> AST2;
    typedef Number_type_checker<NT1, NT2, Cmp> Type;
public:
    // CGAL::Algebraic_structure_traits< >::Div
    class Div
        : public CGAL::cpp98::binary_function< Type,
                                  Type,
                                  Type > {
    public:
        Type operator()(
                const Type& a,
                const Type& b) const {
            CGAL_NT_CHECK_DEBUG("AST::Div");
            return Type(
                    typename AST1::Div()(a.n1(),b.n1()),
                    typename AST2::Div()(a.n2(),b.n2()));
        }
    };
    // CGAL::Algebraic_structure_traits< >::Mod
    class Mod
        : public CGAL::cpp98::binary_function< Type,
                                  Type,
                                  Type > {
    public:
        Type operator()(
                const Type& a,
                const Type& b) const {
            CGAL_NT_CHECK_DEBUG("AST::Mod");
            return Type(
                    typename AST1::Mod()(a.n1(),b.n1()),
                    typename AST2::Mod()(a.n2(),b.n2()));
        }
    };

    // CGAL::Algebraic_structure_traits< >::Div_mod
    class Div_mod {
    public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef Type&   fourth_argument_type;
        typedef void  result_type;

        void operator()(
                const Type& a,
                const Type& b,
                Type& q,
                Type& r) const {
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


template < typename NT1, typename NT2, typename Cmp >
class NTC_AST_base
      < Number_type_checker< NT1, NT2, Cmp> , Field_tag >
      :public  NTC_AST_base
      < Number_type_checker< NT1, NT2, Cmp> , Integral_domain_tag >
{
private:
  typedef Algebraic_structure_traits<NT1> AST1;
  typedef Algebraic_structure_traits<NT2> AST2;
  typedef Number_type_checker<NT1, NT2, Cmp> Type;
public:
  class Inverse
    : public CGAL::cpp98::unary_function< Type, Type > {
  public:
    Type operator()( const Type& a ) const {
      NT1 r1 = typename AST1::Inverse()(a.n1());
      NT2 r2 = typename AST2::Inverse()(a.n2());
      return Type(r1,r2);
    }
  };


};

template < typename NT1, typename NT2, typename Cmp >
class NTC_AST_base
      < Number_type_checker< NT1, NT2, Cmp> , Field_with_sqrt_tag >
      :public  NTC_AST_base
      < Number_type_checker< NT1, NT2, Cmp> , Field_tag >
{
private:
    typedef Algebraic_structure_traits<NT1> AST1;
    typedef Algebraic_structure_traits<NT2> AST2;
    typedef Number_type_checker<NT1, NT2, Cmp> Type;

public:

    // CGAL::Algebraic_structure_traits<  >::Sqrt
    class Sqrt
        : public CGAL::cpp98::unary_function< Type , Type > {
    public:
        Type operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("AST::Sqrt");
            return Type(
                    typename AST1::Sqrt()(a.n1()),
                    typename AST2::Sqrt()(a.n2()));
        }
    };
};

} // namespace NTC_INTERN


template < typename NT1, typename NT2, typename Cmp >
class Algebraic_structure_traits <Number_type_checker<NT1, NT2, Cmp> >
  :public NTC_INTERN::NTC_AST_base< Number_type_checker< NT1, NT2, Cmp> ,
                                 typename Algebraic_structure_traits<NT1>::Algebraic_category >
{
    typedef Algebraic_structure_traits<NT1> AST1;
public:
    typedef Number_type_checker< NT1, NT2, Cmp> Type;
    typedef typename AST1::Algebraic_category Algebraic_category;
    typedef typename AST1::Is_exact Is_exact;

};


namespace NTC_INTERN{
template < typename Number_type_checker, typename Is_real_embeddable >
class NTC_RET_base;

template < typename NT >
class NTC_RET_base<NT,Tag_false> : public Real_embeddable_traits<NT>
{};

template < typename NT1, typename NT2, typename Cmp >
class NTC_RET_base
< Number_type_checker<NT1, NT2, Cmp> , Tag_true>
  :public INTERN_RET::Real_embeddable_traits_base< Number_type_checker< NT1, NT2, Cmp > , CGAL::Tag_true >
{
private:
    typedef Real_embeddable_traits<NT1> RET1;
    typedef Real_embeddable_traits<NT2> RET2;
    typedef Number_type_checker<NT1, NT2, Cmp> Type;
public:

    // CGAL::Real_embeddable_traits<  >::Abs
    class Abs
        : public CGAL::cpp98::unary_function< Type , Type > {
    public:
        Type operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("RET::Abs");
            return Type(
                    typename RET1::Abs()(a.n1()),
                    typename RET2::Abs()(a.n2()));
        }
    };

    // CGAL::Real_embeddable_traits<  >::Sign
    class Sgn
        : public CGAL::cpp98::unary_function< Type , ::CGAL::Sign > {
    public:
        ::CGAL::Sign operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("RET::Sign");
            ::CGAL::Sign r1 =  typename RET1::Sgn()(a.n1());
            CGAL_assertion( r1 == typename RET2::Sgn()(a.n2()) );
            return r1;
        }
    };

    // CGAL::Real_embeddable_traits<  >::Is_finite
    class Is_finite
        : public CGAL::cpp98::unary_function< Type , bool > {
    public:
        bool operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("RET::Is_finite");
            bool r1 =  typename RET1::Is_finite()(a.n1());
            CGAL_assertion( r1 == typename RET2::Is_finite()(a.n2()) );
            return r1;
        }
    };

    // CGAL::Real_embeddable_traits<  >::Is_positive
    class Is_positive
        : public CGAL::cpp98::unary_function< Type , bool > {
    public:
        bool operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("RET::Is_positive");
            bool r1 =  typename RET1::Is_positive()(a.n1());
            CGAL_assertion( r1 == typename RET2::Is_positive()(a.n2()) );
            return r1;
        }
    };

    // CGAL::Real_embeddable_traits<  >::Is_negative
    class Is_negative
        : public CGAL::cpp98::unary_function< Type , bool > {
    public:
        bool operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("RET::Is_negative");
            bool r1 =  typename RET1::Is_negative()(a.n1());
            CGAL_assertion( r1 == typename RET2::Is_negative()(a.n2()) );
            return r1;
        }
    };

    // CGAL::Real_embeddable_traits<  >::Is_zero
    class Is_zero
        : public CGAL::cpp98::unary_function< Type , bool > {
    public:
        bool operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("RET::Is_zero");
            bool r1 =  typename RET1::Is_zero()(a.n1());
            CGAL_assertion( r1 == typename RET2::Is_zero()(a.n2()) );
            return r1;
        }
    };

    // CGAL::Real_embeddable_traits<  >::Compare
    class Compare
        : public CGAL::cpp98::binary_function< Type , Type, Comparison_result > {
    public:
        Comparison_result operator()(const Type& a, const Type& b) const {
            CGAL_NT_CHECK_DEBUG("RET::Compare");
            Comparison_result r1 =  typename RET1::Compare()(a.n1(),b.n1());
            CGAL_assertion( r1 == typename RET2::Compare()(a.n2(),b.n2()) );
            return r1;
        }
    };

    // CGAL::Real_embeddable_traits<  >::To_double
    class To_double
        : public CGAL::cpp98::unary_function< Type , double > {
    public:
        double operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("RET::To_double");
            double r1 =  typename RET1::To_double()(a.n1());
            CGAL_assertion( r1 == typename RET2::To_double()(a.n2()) );
            return r1;
        }
    };

    // CGAL::Real_embeddable_traits<  >::To_interval
    class To_interval
        : public CGAL::cpp98::unary_function< Type , std::pair<double, double> > {
    public:
        std::pair<double, double> operator()(const Type& a) const {
            CGAL_NT_CHECK_DEBUG("RET::To_interval");
            std::pair<double, double> r1 =  typename RET1::To_interval()(a.n1());
            CGAL_assertion( r1 == typename RET2::To_interval()(a.n2()) );
            return r1;
        }
    };
};

} // namespace NTC_INTERN

template < typename NT1, typename NT2, typename Cmp >
class Real_embeddable_traits
< Number_type_checker<NT1, NT2, Cmp> >
    :public NTC_INTERN::NTC_RET_base< Number_type_checker< NT1, NT2, Cmp > ,
      typename Real_embeddable_traits<NT1>::Is_real_embeddable >
{
    typedef Real_embeddable_traits<NT1> RET1;
public:
    typedef Number_type_checker< NT1, NT2, Cmp >  Type;
    typedef typename RET1::Is_real_embeddable     Is_real_embeddable;
};

template < typename NT1, typename NT2, typename Cmp >
struct Coercion_traits< Number_type_checker<NT1,NT2,Cmp>,
                        Number_type_checker<NT1,NT2,Cmp> >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef Number_type_checker<NT1,NT2,Cmp> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const Type& x) const { return x;}
    };
};

template < typename NT1, typename NT2, typename Cmp >
struct Coercion_traits< Number_type_checker<NT1,NT2,Cmp>, int >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef Number_type_checker<NT1,NT2,Cmp> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const Type& x) const { return x;}
        Type operator()(int x) const { return Type(x);}
    };
};

template < typename NT1, typename NT2, typename Cmp >
struct Coercion_traits< int , Number_type_checker<NT1,NT2,Cmp> >
    :public Coercion_traits< Number_type_checker<NT1,NT2,Cmp> , int >{};

namespace NTC_INTERN {


template < typename NT_checker, typename Tag = Tag_false >
struct Coercion_traits_double{
    typedef Tag_false  Are_explicit_interoperable;
    typedef Tag_false  Are_implicit_interoperable;
    typedef Null_tag Type;
};

template < typename NT1, typename NT2, typename Cmp>
struct Coercion_traits_double< Number_type_checker<NT1,NT2,Cmp> ,
                               Tag_true >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef Number_type_checker<NT1,NT2,Cmp> Type;
     struct Cast{
        typedef Type result_type;
         Type operator()(const Type& x) const { return x;}
         Type operator()(const double& x) const {
             return Type(x);
         }
    };
};
} // namespace NTC_INTERN

template < typename NT1, typename NT2, typename Cmp >
struct Coercion_traits< Number_type_checker<NT1,NT2,Cmp>, double >
    :public NTC_INTERN::Coercion_traits_double< Number_type_checker<NT1,NT2,Cmp>,
            typename Coercion_traits<NT1,double>::Are_implicit_interoperable >
{};

template < typename NT1, typename NT2, typename Cmp >
struct Coercion_traits< double , Number_type_checker<NT1,NT2,Cmp> >
    :public Coercion_traits< Number_type_checker<NT1,NT2,Cmp>, double > {};


// What to do with the IO ?
// - either pick one as the main, and convert to the second.
// - output/input pairs, and read both (=> problems with FP values).

template < typename NT1, typename NT2, typename Cmp >
std::ostream &
operator<< (std::ostream & os, const Number_type_checker<NT1, NT2, Cmp> &b)
{
  return os << b.n1();
}

template < typename NT1, typename NT2, typename Cmp >
std::istream &
operator>> (std::istream & is, Number_type_checker<NT1, NT2, Cmp> &b)
{
  is >> b.n1();
  std::stringstream ss;
  ss << b.n1();
  ss >> b.n2();
  return is;
}

} //namespace CGAL

#endif // CGAL_NUMBER_TYPE_CHECKER_H
