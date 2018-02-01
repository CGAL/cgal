// Copyright (c) 2001,2007  
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
// Author(s)     : Geert-Jan Giezeman,
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>

#ifndef CGAL_COUNTED_NUMBER_H
#define CGAL_COUNTED_NUMBER_H

#include <CGAL/number_type_basic.h>
#include <CGAL/atomic.h>
#include <CGAL/boost/iterator/transform_iterator.hpp> // for Root_of_selector
#include <iostream>

namespace CGAL {

template <class NT>
class Counted_number {
#ifdef CGAL_NO_ATOMIC
    static unsigned long
#else
    static CGAL::cpp11::atomic<unsigned long>
#endif
                         s_neg_count, s_add_count, s_sub_count,
                         s_mul_count, s_div_count,
                  			 s_eq_count, s_comp_count,
                         s_simplify_count,
                         s_unit_part_count,
                         s_is_zero_count,
                         s_is_one_count,
                         s_square_count,
                         s_integral_division_count,
                         s_is_square_count,
                         s_sqrt_count,
                         s_kth_root_count,
                         s_root_of_count,
                         s_gcd_count,
                         s_div_mod_count,
                         s_mod_count;
    NT m_rep;
  public:
    typedef NT Rep_type;
    static void reset()
            { s_neg_count=0; s_add_count=0; s_sub_count=0;
              s_mul_count=0; s_div_count=0;
      	      s_eq_count=0; s_comp_count = 0;
              s_simplify_count = 0; s_unit_part_count = 0; s_is_zero_count = 0;
              s_is_one_count = 0; s_square_count = 0;
              s_integral_division_count = 0; s_is_square_count = 0;
              s_sqrt_count = 0; s_kth_root_count = 0; s_root_of_count = 0;
              s_gcd_count = 0; s_div_mod_count = 0; s_mod_count = 0;
            }
    static void inc_neg_count() {++s_neg_count;}
    static void inc_add_count() {++s_add_count;}
    static void inc_sub_count() {++s_sub_count;}
    static void inc_mul_count() {++s_mul_count;}
    static void inc_div_count() {++s_div_count;}
    static void inc_eq_count() {++s_eq_count;}
    static void inc_comp_count() {++s_comp_count;}
    static void inc_simplify_count() {++s_simplify_count;}
    static void inc_unit_part_count() {++s_unit_part_count;}
    static void inc_is_zero_count() {++s_is_zero_count;}
    static void inc_is_one_count() {++s_is_one_count;}
    static void inc_square_count() {++s_square_count;}
    static void inc_integral_division_count() {++s_integral_division_count;}
    static void inc_is_square_count() {++s_is_square_count;}
    static void inc_sqrt_count() {++s_sqrt_count;}
    static void inc_kth_root_count() {++s_kth_root_count;}
    static void inc_root_of_count() {++s_root_of_count;}
    static void inc_gcd_count() {++s_gcd_count;}
    static void inc_div_mod_count() {++s_div_mod_count;}
    static void inc_mod_count() {++s_mod_count;}

    static unsigned long neg_count() {return s_neg_count;}
    static unsigned long add_count() {return s_add_count;}
    static unsigned long sub_count() {return s_sub_count;}
    static unsigned long mul_count() {return s_mul_count;}
    static unsigned long div_count() {return s_div_count;}
    static unsigned long eq_count() {return s_eq_count;}
    static unsigned long comp_count() {return s_comp_count;}
    static unsigned long simplify_count() {return s_simplify_count;}
    static unsigned long unit_part_count() {return s_unit_part_count;}
    static unsigned long is_zero_count() {return s_is_zero_count;}
    static unsigned long is_one_count() {return s_is_one_count;}
    static unsigned long square_count() {return s_square_count;}
    static unsigned long integral_division_count() {
      return s_integral_division_count;
    }
    static unsigned long is_square_count() {return s_is_square_count;}
    static unsigned long sqrt_count() {return s_sqrt_count;}
    static unsigned long kth_root_count() {return s_kth_root_count;}
    static unsigned long root_of_count() {return s_root_of_count;}
    static unsigned long gcd_count() {return s_gcd_count;}
    static unsigned long div_mod_count() {return s_div_mod_count;}
    static unsigned long mod_count() {return s_mod_count;}

    static unsigned long count()
            { return s_neg_count + s_add_count + s_sub_count +
                     s_mul_count + s_div_count +
      	             s_eq_count + s_comp_count +
                     s_simplify_count + s_unit_part_count + s_is_zero_count +
                     s_is_one_count + s_square_count +
                     s_integral_division_count + s_is_square_count +
                     s_sqrt_count + s_kth_root_count + s_root_of_count +
                     s_gcd_count + s_div_mod_count + s_mod_count;
            }

    static void report(std::ostream &os);
    NT rep() const {return m_rep;}
    Counted_number() {}
    //explicit Counted_number(int n) :m_rep(n){}
    explicit Counted_number(NT n) :m_rep(n){}
    Counted_number operator-() const
            {inc_neg_count();return Counted_number(-m_rep);}
    Counted_number const & operator+=(Counted_number const &n)
            {
		inc_add_count();
		m_rep += n.m_rep;
		return *this;}
    Counted_number const & operator-=(Counted_number const &n)
            {inc_sub_count(); m_rep -= n.m_rep; return *this;}
    Counted_number const & operator*=(Counted_number const &n)
            {inc_mul_count(); m_rep *= n.m_rep; return *this;}
    Counted_number const & operator/=(Counted_number const &n)
            {inc_div_count(); m_rep /= n.m_rep; return *this;}

    // Counted operations
    void simplify() {
      inc_simplify_count();
      CGAL_NTS simplify( m_rep );
    }

    Counted_number unit_part() const {
      inc_unit_part_count();
      return Counted_number( CGAL_NTS unit_part( rep() ) );
    }

    bool is_zero() const {
      inc_is_zero_count();
      return CGAL_NTS is_zero( rep() );
    }

    bool is_one() const {
      inc_is_one_count();
      return CGAL_NTS is_one( rep() );
    }

    Counted_number square() const {
      inc_square_count();
      return Counted_number( CGAL_NTS square( rep() ) );
    }

    Counted_number integral_division( const Counted_number& n ) const {
      inc_integral_division_count();
      return Counted_number( CGAL_NTS integral_division( rep(), n.rep() ) );
    }

    bool is_square( Counted_number& result ) const {
      inc_is_square_count();
      NT result_as_nt;
      bool is_integral = CGAL_NTS is_square( rep(), result_as_nt );
      result = Counted_number( result_as_nt );
      return is_integral;
    }

    Counted_number sqrt() const {
      inc_sqrt_count();
      return Counted_number( CGAL_NTS sqrt( rep() ) );
    }

    Counted_number kth_root( int k ) const {
      inc_kth_root_count();
      return Counted_number( CGAL_NTS kth_root( k, rep() ) );
    }

    Counted_number gcd( const Counted_number& n ) const {
      inc_gcd_count();
      return Counted_number( CGAL_NTS gcd( rep(), n.rep() ) );
    }

    Counted_number div( const Counted_number& n ) const {
      inc_div_count();
      return Counted_number( CGAL_NTS div( rep(), n.rep() ) );
    }

    Counted_number mod( const Counted_number& n ) const {
      inc_mod_count();
      return Counted_number( CGAL_NTS mod( rep(), n.rep() ) );
    }

    void div_mod( const Counted_number& n, Counted_number& q,
                  Counted_number& r ) const {
      inc_div_mod_count();
      NT q_as_nt, r_as_nt;
      CGAL_NTS div_mod( rep(), n.rep(), q_as_nt, r_as_nt );
      q = Counted_number( q_as_nt );
      r = Counted_number( r_as_nt );
    }

    // Other operations
    inline double to_double() const {
      return CGAL_NTS to_double( rep() );
    }

    inline std::pair<double, double> to_interval() const {
      return CGAL_NTS to_interval( rep() );
    }
};

#ifdef CGAL_NO_ATOMIC
template <class NT>
unsigned long Counted_number<NT>::s_neg_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_add_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_sub_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_mul_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_div_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_eq_count=0;

template <class NT>
unsigned long Counted_number<NT>::s_comp_count=0;

template< class NT >
unsigned long Counted_number<NT>::s_simplify_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_unit_part_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_is_zero_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_is_one_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_square_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_integral_division_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_is_square_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_sqrt_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_kth_root_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_root_of_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_gcd_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_div_mod_count = 0;

template< class NT >
unsigned long Counted_number<NT>::s_mod_count = 0;
#else 
template <class NT>
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_neg_count;

template <class NT>
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_add_count;

template <class NT>
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_sub_count;

template <class NT>
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_mul_count;

template <class NT>
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_div_count;

template <class NT>
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_eq_count;

template <class NT>
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_comp_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_simplify_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_unit_part_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_is_zero_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_is_one_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_square_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_integral_division_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_is_square_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_sqrt_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_kth_root_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_root_of_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_gcd_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_div_mod_count;

template< class NT >
CGAL::cpp11::atomic<unsigned long> Counted_number<NT>::s_mod_count;
#endif


//unary +
template <class NT> Counted_number<NT>
operator + (const Counted_number<NT>& n1){
    return n1;
}

template <class NT>
Counted_number<NT>
operator+(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_add_count();
    return Counted_number<NT>(n1.rep() + n2.rep());
}

template <class NT>
Counted_number<NT>
operator-(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_sub_count();
    return Counted_number<NT>(n1.rep() - n2.rep());
}

template <class NT>
Counted_number<NT>
operator*(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_mul_count();
    return Counted_number<NT>(n1.rep() * n2.rep());
}

template <class NT>
Counted_number<NT>
operator/(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_div_count();
    return Counted_number<NT>(n1.rep() / n2.rep());
}

template< class NT >
Counted_number<NT>
operator%( const Counted_number<NT>& x, const Counted_number<NT>& y ) {
  Counted_number<NT>::inc_mod_count();
  return Counted_number<NT>( x.rep() % y.rep() );
}

template <class NT>
bool
operator==(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_eq_count();
    return (n1.rep() == n2.rep());
}

template <class NT>
bool
operator!=(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_eq_count();
    return (n1.rep() != n2.rep());
}

template <class NT>
bool
operator<(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_comp_count();
    return (n1.rep() < n2.rep());
}

template <class NT>
bool
operator>(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_comp_count();
    return (n1.rep() > n2.rep());
}

template <class NT>
bool
operator<=(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_comp_count();
    return (n1.rep() <= n2.rep());
}

template <class NT>
bool
operator>=(Counted_number<NT> const &n1, Counted_number<NT> const &n2)
{
    Counted_number<NT>::inc_comp_count();
    return (n1.rep() >= n2.rep());
}

template <class NT>
class Is_valid< Counted_number<NT> >
  : public CGAL::unary_function< Counted_number<NT>, bool > {
  public:
    bool operator()( const Counted_number<NT>& x ) {
      return is_valid( x.rep() );
    }
};

template <class NT>
void Counted_number<NT>::report(std::ostream &os)
{
    os << count() << " operations\n";
    if (neg_count() > 0)
        os << "  " << neg_count() << " negations\n";
    if (add_count() > 0)
        os << "  " << add_count() << " additions\n";
    if (sub_count() > 0)
        os << "  " << sub_count() << " subtractions\n";
    if (mul_count() > 0)
        os << "  " << mul_count() << " multiplications\n";
    if (div_count() > 0)
        os << "  " << div_count() << " divisions\n";
    if (eq_count() > 0)
        os << "  " << eq_count() << " equality tests\n";
    if (comp_count() > 0)
        os << "  " << comp_count() << " comparisons\n";
    if (simplify_count() > 0)
        os << "  " << simplify_count() << " simplify-calls\n";
    if (unit_part_count() > 0)
        os << "  " << unit_part_count() << " unit_part-calls\n";
    if (is_zero_count() > 0)
        os << "  " << is_zero_count() << " is_zero-calls\n";
    if (is_one_count() > 0)
        os << "  " << is_one_count() << " is_one-calls\n";
    if (square_count() > 0)
        os << "  " << square_count() << " square-calls\n";
    if (integral_division_count() > 0)
        os << "  " << integral_division_count() << " integral_division-calls\n";
    if (is_square_count() > 0)
        os << "  " << is_square_count() << " is_square-calls\n";
    if (kth_root_count() > 0)
        os << "  " << kth_root_count() << " kth_root-calls\n";
    if (root_of_count() > 0)
        os << "  " << root_of_count() << " root_of-calls\n";
    if (gcd_count() > 0)
        os << "  " << gcd_count() << " gcd-calls\n";
    if (div_mod_count() > 0)
        os << "  " << div_mod_count() << " div_mod-calls\n";
    if (mod_count() > 0)
        os << "  " << mod_count() << " mod-calls\n";
}

template <class NT>
std::ostream& operator<<(std::ostream &os, Counted_number<NT> const &n)
{
    return os << ::CGAL::oformat( n.rep() )<< std::endl;
}

template <class NT>
std::istream& operator>>(std::istream &is, Counted_number<NT> &n)
{
    NT num;
    is >> ::CGAL::iformat(num);
    if (is) n = Counted_number<NT>(num);
    return is;
}

namespace INTERN_COUNTED_NUMBER{

template< class NT, class Functor >
struct Simplify_selector {
  struct Simplify : public CGAL::unary_function<NT&, void> {
    void operator()( NT& x ) const {
      x.simplify();
    }
  };
};

template< class NT >
struct Simplify_selector< NT, Null_functor > {
  typedef Null_functor Simplify;
};

template< class NT, class Functor >
struct Unit_part_selector {
  struct Unit_part : public CGAL::unary_function<NT, NT > {
    NT operator()( const NT& x ) const {
      return x.unit_part();
    }
  };
};

template< class NT >
struct Unit_part_selector< NT, Null_functor > {
  typedef Null_functor Unit_part;
};

template< class NT, class Functor >
struct Is_zero_selector {
  struct Is_zero : public CGAL::unary_function<NT, bool > {
    bool operator()( const NT& x ) const {
      return x.is_zero();
    }
  };
};

template< class NT >
struct Is_zero_selector< NT, Null_functor > {
  typedef Null_functor Is_zero;
};

template< class NT, class Functor >
struct Is_one_selector {
  struct Is_one : public CGAL::unary_function<NT, bool > {
    bool operator()( const NT& x ) const {
      return x.is_one();
    }
  };
};

template< class NT >
struct Is_one_selector< NT, Null_functor > {
  typedef Null_functor Is_one;
};

template< class NT, class Functor >
struct Square_selector {
  struct Square : public CGAL::unary_function<NT, NT > {
    NT operator()( const NT& x ) const {
      return x.square();
    }
  };
};

template< class NT >
struct Square_selector< NT, Null_functor > {
  typedef Null_functor Square;
};

template< class NT, class Functor >
struct Integral_division_selector {
  struct Integral_division : public CGAL::binary_function<NT, NT, NT > {
    NT operator()( const NT& x, const NT& y ) const {
      return x.integral_division( y );
    }
  };
};

template< class NT >
struct Integral_division_selector< NT, Null_functor > {
  typedef Null_functor Integral_division;
};

template< class NT, class Functor >
struct Is_square_selector {
  struct Is_square : public CGAL::binary_function<NT, NT&, bool > {
      bool operator()( const NT& x, NT& y ) const {
          return x.is_square( y );
      }
      bool operator()( const NT& x) const {
          NT y;
          return x.is_square( y );
      }
  };
};

template< class NT >
struct Is_square_selector< NT, Null_functor > {
  typedef Null_functor Is_square;
};


template <class NT, class AlgebraicStructureTag>
struct Sqrt_selector{
    struct Sqrt : public CGAL::unary_function<NT,NT> {
        NT operator ()(const NT& x) const {
            return x.sqrt();
        }
    };
};
template <class NT>
struct Sqrt_selector<NT,Null_functor> {
    typedef Null_functor Sqrt;
};

template< class NT, class Functor >
struct Kth_root_selector {
  struct Kth_root : public CGAL::binary_function<int, NT, NT > {
    NT operator()( int k, const NT& x ) const {
      return x.kth_root( k );
    }
  };
};

template< class NT >
struct Kth_root_selector< NT, Null_functor > {
  typedef Null_functor Kth_root;
};

template< class NT, class Functor >
struct Root_of_selector {
  private:
    typedef typename NT::Rep_type Rep_type;
    struct Cast{
      typedef Rep_type result_type;
      result_type operator()(const NT& counted_number) const {
        return counted_number.rep();
      }
    };

  public:
    struct Root_of {
//      typedef typename Functor::Boundary Boundary;
      typedef NT result_type;
      template< class Input_iterator >
      NT operator()( int k, Input_iterator begin, Input_iterator end ) const {
        NT::inc_root_of_count();
        Cast cast;
        return NT( Functor()( k,
                              ::boost::make_transform_iterator( begin, cast ),
                              ::boost::make_transform_iterator( end, cast ) ) );
      }

      // TODO: Why are the arguments not const-ref?
/*      template< class Input_iterator >
      NT operator()( Boundary lower, Boundary upper,
                     Input_iterator begin, Input_iterator end ) const {
        NT::inc_root_of_count();
        Cast cast;
        return NT( Functor()( lower, upper,
                             ::boost::make_transform_iterator( begin, cast ),
                             ::boost::make_transform_iterator( end, cast ) ) );
      }*/
    };
};

template< class NT >
struct Root_of_selector< NT, Null_functor > {
  typedef Null_functor Root_of;
};

template< class NT, class Functor >
struct Gcd_selector {
  struct Gcd : public CGAL::binary_function<NT, NT, NT > {
    NT operator()( const NT& x, const NT& y ) const {
      return x.gcd( y );
    }
  };
};

template< class NT >
struct Gcd_selector< NT, Null_functor > {
  typedef Null_functor Gcd;
};

template< class NT, class Functor >
struct Div_selector {
  struct Div : public CGAL::binary_function<NT, NT, NT > {
    NT operator()( const NT& x, const NT& y ) const {
      return x.div( y );
    }
  };
};

template< class NT >
struct Div_selector< NT, Null_functor > {
  typedef Null_functor Div;
};

template< class NT, class Functor >
struct Mod_selector {
  struct Mod : public CGAL::binary_function<NT, NT, NT > {
    NT operator()( const NT& x, const NT& y ) const {
      return x.mod( y );
    }
  };
};

template< class NT >
struct Mod_selector< NT, Null_functor > {
  typedef Null_functor Mod;
};

template< class NT, class Functor >
struct Div_mod_selector {
  struct Div_mod {
    typedef void result_type;
    typedef NT   first_argument_type;
    typedef NT   second_argument_type;
    typedef NT&  third_argument_type;
    typedef NT&  fourth_argument_type;

    void operator()( const NT& x, const NT& y, NT& q, NT& r ) const {
      x.div_mod( y, q, r );
    }
  };
};

template< class NT >
struct Div_mod_selector< NT, Null_functor >{
  typedef Null_functor Div_mod;
};

} // end namespace INTERN_COUNTED_NUMBER

template <class NT>
class Algebraic_structure_traits<Counted_number<NT> >
    :public Algebraic_structure_traits_base
      <Counted_number<NT>,
       typename Algebraic_structure_traits<NT>::Algebraic_category >
{
private:
    typedef Algebraic_structure_traits<NT> AST_NT;
    typedef typename AST_NT::Algebraic_category NT_as_tag;

public:
    typedef typename Algebraic_structure_traits<NT>::Is_exact Is_exact;
    typedef typename AST_NT::Is_numerical_sensitive Is_numerical_sensitive;

    typedef typename INTERN_COUNTED_NUMBER::Simplify_selector
    <Counted_number<NT>, typename AST_NT::Simplify > ::Simplify Simplify;

    typedef typename INTERN_COUNTED_NUMBER::Unit_part_selector
    <Counted_number<NT>, typename AST_NT::Unit_part > ::Unit_part Unit_part;

    typedef typename INTERN_COUNTED_NUMBER::Is_zero_selector
    <Counted_number<NT>, typename AST_NT::Is_zero > ::Is_zero Is_zero;

    typedef typename INTERN_COUNTED_NUMBER::Is_one_selector
    <Counted_number<NT>, typename AST_NT::Is_one > ::Is_one Is_one;

    typedef typename INTERN_COUNTED_NUMBER::Square_selector
    <Counted_number<NT>, typename AST_NT::Square > ::Square Square;

    typedef typename INTERN_COUNTED_NUMBER::Integral_division_selector
    <Counted_number<NT>, typename AST_NT::Integral_division> ::Integral_division Integral_division;

    typedef typename INTERN_COUNTED_NUMBER::Is_square_selector
    <Counted_number<NT>, typename AST_NT::Is_square > ::Is_square Is_square;

    typedef typename INTERN_COUNTED_NUMBER::Sqrt_selector
    <Counted_number<NT>, typename AST_NT::Sqrt> ::Sqrt Sqrt;

    typedef typename INTERN_COUNTED_NUMBER::Kth_root_selector
    <Counted_number<NT>, typename AST_NT::Kth_root > ::Kth_root Kth_root;

    typedef typename INTERN_COUNTED_NUMBER::Root_of_selector
    <Counted_number<NT>, typename AST_NT::Root_of > ::Root_of Root_of;

    typedef typename INTERN_COUNTED_NUMBER::Gcd_selector
    <Counted_number<NT>, typename AST_NT::Gcd > ::Gcd Gcd;

    typedef typename INTERN_COUNTED_NUMBER::Div_selector
    <Counted_number<NT>, typename AST_NT::Div > ::Div Div;

    typedef typename INTERN_COUNTED_NUMBER::Mod_selector
    <Counted_number<NT>, typename AST_NT::Mod > ::Mod Mod;

    typedef typename INTERN_COUNTED_NUMBER::Div_mod_selector
    <Counted_number<NT>, typename AST_NT::Div_mod > ::Div_mod Div_mod;
};

template <class NT>
class Real_embeddable_traits<Counted_number<NT> >
  : public INTERN_RET::Real_embeddable_traits_base <Counted_number<NT> , 
   typename Real_embeddable_traits<NT>::Is_real_embeddable > 
{
    typedef Real_embeddable_traits<NT> RET_NT;

public:
    typedef typename INTERN_COUNTED_NUMBER::Is_zero_selector
    <Counted_number<NT>, typename RET_NT::Is_zero > ::Is_zero Is_zero;

    class Is_finite
      : public CGAL::unary_function< Counted_number<NT>, bool > {
      public:
        bool operator()( const Counted_number<NT>& x ) const {
          return CGAL_NTS is_finite( x.rep() );
        }
    };

    struct To_double : public CGAL::unary_function< Counted_number<NT>, double > {
        double operator()(const Counted_number<NT>& x) const {
            return x.to_double();
        }
    };

    struct To_interval: public CGAL::unary_function< Counted_number<NT>, std::pair<double,double> > {
        std::pair<double,double>
        operator()(const Counted_number<NT>& x) const {
            return x.to_interval();
        }
    };
};

template<typename NT> inline 
Counted_number<NT> min BOOST_PREVENT_MACRO_SUBSTITUTION(
const Counted_number<NT> & x,
const Counted_number<NT> & y){
  return CGAL::Min<Counted_number<NT> > ()(x,y);
}
template<typename NT> inline 
Counted_number<NT> max BOOST_PREVENT_MACRO_SUBSTITUTION(
const Counted_number<NT> & x,
const Counted_number<NT> & y){
  return CGAL::Max<Counted_number<NT> > ()(x,y);
}

} //namespace CGAL

#endif
