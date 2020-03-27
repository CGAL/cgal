// Copyright (c) 1999,2007
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
// Author(s)     : Stefan Schirra, Michael Hemmer


#ifndef CGAL_INT_H
#define CGAL_INT_H

#include <CGAL/Interval_nt.h>
#include <CGAL/Modular_traits.h>
#include <CGAL/Modular_arithmetic/Residue_type.h>

namespace CGAL {



namespace INTERN_INT {
    template< class Type >
    class Is_square_per_double_conversion
      : public CGAL::cpp98::binary_function< Type, Type&,
                                bool > {
      public:
        bool operator()( const Type& x,
                         Type& y ) const {
          y = (Type) std::sqrt( (double)x );
          return x == y * y;
        }
        bool operator()( const Type& x ) const {
            Type y =
                (Type) std::sqrt( (double)x );
            return x == y * y;
        }

    };
} // INTERN_INT

// int
template<> class Algebraic_structure_traits< int >
  : public Algebraic_structure_traits_base< int, Euclidean_ring_tag > {

  public:
    typedef Tag_false            Is_exact;
    typedef Tag_true             Is_numerical_sensitive;

    typedef INTERN_AST::Div_per_operator< Type >  Div;
    typedef INTERN_AST::Mod_per_operator< Type >  Mod;

    typedef INTERN_INT::
       Is_square_per_double_conversion< Type > Is_square;
};

template <> class Real_embeddable_traits< int >
  : public INTERN_RET::Real_embeddable_traits_base< int , CGAL::Tag_true > {};

/*! \ingroup CGAL_Modular_traits_spec
  \brief Specialization of CGAL::Modular_traits for \c int.

  A model of concept ModularTraits, supports \c int.
*/
  template <typename T>
  class Modular_traits;

template<>
class Modular_traits<int>{
public:
    typedef int NT;
    typedef ::CGAL::Tag_true Is_modularizable;
    typedef Residue Residue_type;

    struct Modular_image{
        Residue_type operator()(int i){
            return Residue_type(i);
        }
    };
    struct Modular_image_representative{
        NT operator()(const Residue_type& x){
            return x.get_value();
        }
    };
};

// long

template<> class Algebraic_structure_traits< long int >
  : public Algebraic_structure_traits_base< long int,
                                            Euclidean_ring_tag > {

  public:
    typedef Tag_false            Is_exact;
    typedef Tag_true           Is_numerical_sensitive;

    typedef INTERN_AST::Div_per_operator< Type >  Div;
    typedef INTERN_AST::Mod_per_operator< Type >  Mod;

    typedef INTERN_INT::
       Is_square_per_double_conversion< Type > Is_square;
};

template <> class Real_embeddable_traits< long int >
  : public INTERN_RET::Real_embeddable_traits_base< long int , CGAL::Tag_true > {
public:

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          return Interval_nt<true>(x).pair();
        }
    };
};


/*! \ingroup CGAL_Modular_traits_spec
  \brief Specialization of CGAL::Modular_traits for \c long.

  A model of concept ModularTraits, supports \c long.
*/
template<>
class Modular_traits<long>{
public:
    typedef long NT;
    typedef ::CGAL::Tag_true Is_modularizable;
    typedef Residue Residue_type;

    struct Modular_image{
        Residue_type operator()(long i){
            return Residue_type(i);
        }
    };
    struct Modular_image_representative{
        NT operator()(const Residue_type& x){
            return NT(x.get_value());
        }
    };
};

// short

template<> class Algebraic_structure_traits< short int >
  : public Algebraic_structure_traits_base< short int,
                                            Euclidean_ring_tag > {

  public:
    typedef Tag_false            Is_exact;
    typedef Tag_true             Is_numerical_sensitive;

    // Explicitly defined functors which have no support for implicit
    //  interoperability. This is nescessary because of the implicit conversion
    //  to int for binary operations between short ints.
    class Integral_division
      : public CGAL::cpp98::binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y) const {
          Algebraic_structure_traits<Type>::Div actual_div;
          CGAL_precondition_msg( actual_div( x, y) * y == x,
                  "'x' must be divisible by 'y' in "
                  "Algebraic_structure_traits<...>::Integral_div()(x,y)" );
          return actual_div( x, y);
        }
    };

    class Gcd
      : public CGAL::cpp98::binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y) const {
          Algebraic_structure_traits<Type>::Mod mod;
          Algebraic_structure_traits<Type>::Unit_part unit_part;
          Algebraic_structure_traits<Type>::Integral_division integral_div;
          // First: the extreme cases and negative sign corrections.
          if (x == Type(0)) {
              if (y == Type(0))
                  return Type(0);
              return integral_div( y, unit_part(y) );
          }
          if (y == Type(0))
              return integral_div(x, unit_part(x) );
          Type u = integral_div( x, unit_part(x) );
          Type v = integral_div( y, unit_part(y) );
          // Second: assuming mod is the most expensive op here, we don't compute it
          // unnecessarily if u < v
          if (u < v) {
              v = mod(v,u);
              // maintain invariant of v > 0 for the loop below
              if ( v == Type(0) )
                  return u;
          }

          Type w;
          do {
              w = mod(u,v);
              if ( w == Type(0))
                  return v;
              u = mod(v,w);
              if ( u == Type(0))
                  return w;
              v = mod(w,u);
          } while (v != Type(0));
          return u;
        }
    };

    class Div_mod {
      public:
        typedef Type    first_argument_type;
        typedef Type    second_argument_type;
        typedef Type&   third_argument_type;
        typedef Type&   fourth_argument_type;
        typedef void  result_type;
        void operator()( const Type& x,
                         const Type& y,
                         Type& q, Type& r) const {
          q = Type(x / y);
          r = Type(x % y);
          CGAL_assertion(x == q * y + r);
          return;
        }
    };

    // based on \c Div_mod.
    class Div
      : public CGAL::cpp98::binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                         const Type& y) const {
          return Type(x / y);
        };
    };

    // based on \c Div_mod.
    class Mod
      : public CGAL::cpp98::binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                         const Type& y) const {
          return Type(x % y);
        };
    };

    typedef INTERN_INT::
       Is_square_per_double_conversion< Type > Is_square;
};

template <> class Real_embeddable_traits< short int >
  : public INTERN_RET::Real_embeddable_traits_base< short int , CGAL::Tag_true > {};

// unsigned int

template <> class Real_embeddable_traits< unsigned int >
   : public INTERN_RET::Real_embeddable_traits_base< unsigned int , CGAL::Tag_true > {};

// unsigned long

template <> class Real_embeddable_traits< unsigned long >
   : public INTERN_RET::Real_embeddable_traits_base< unsigned long , CGAL::Tag_true > {
public:

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          return Interval_nt<true>(x).pair();
        }
    };
};

// Note : "long long" support is in <CGAL/long_long.h>

} //namespace CGAL

#endif // CGAL_INT_H
