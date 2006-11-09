// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Stefan Schirra, Michael Hemmer
 

#ifndef CGAL_INT_H
#define CGAL_INT_H

#include <CGAL/number_type_basic.h>

CGAL_BEGIN_NAMESPACE

namespace INTERN_INT {
    template< class Algebraic_structure >
    class Is_square_per_double_conversion 
      : public Binary_function< Algebraic_structure, Algebraic_structure&,
                                bool > {
      public:
        bool operator()( const Algebraic_structure& x,
                         Algebraic_structure& y ) const {
          y = (Algebraic_structure) CGAL_CLIB_STD::sqrt( (double)x );
          return x == y * y;
        }
        bool operator()( const Algebraic_structure& x ) const {
            Algebraic_structure y = 
                (Algebraic_structure) CGAL_CLIB_STD::sqrt( (double)x );
            return x == y * y;
        }
        
    };
} // INTERN_INT


// int

template <> struct Number_type_traits<int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;
};

template<> class Algebraic_structure_traits< int >
  : public Algebraic_structure_traits_base< int, Euclidean_ring_tag > {

  public:
    typedef Tag_true            Is_exact;
    
    typedef INTERN_AST::Div_per_operator< Algebraic_structure >  Div;
    typedef INTERN_AST::Mod_per_operator< Algebraic_structure >  Mod;
    
    typedef INTERN_INT::
       Is_square_per_double_conversion< Algebraic_structure > Is_square;
};

template <> class Real_embeddable_traits< int > 
  : public Real_embeddable_traits_base< int > {
  public:
          
    typedef INTERN_RET::To_double_by_conversion< Real_embeddable >
                                                                      To_double;
    typedef INTERN_RET::To_interval_by_conversion< Real_embeddable >
                                                                    To_interval;
};

inline
io_Read_write
io_tag(int)
{ return io_Read_write(); }

// long

template <> struct Number_type_traits<long int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;
};

template<> class Algebraic_structure_traits< long int >
  : public Algebraic_structure_traits_base< long int, 
                                            Euclidean_ring_tag > {

  public:
    typedef Tag_true            Is_exact;
    
    typedef INTERN_AST::Div_per_operator< Algebraic_structure >  Div;
    typedef INTERN_AST::Mod_per_operator< Algebraic_structure >  Mod;       

    typedef INTERN_INT::
       Is_square_per_double_conversion< Algebraic_structure > Is_square;
};

template <> class Real_embeddable_traits< long int > 
  : public Real_embeddable_traits_base< long int > {
  public:
          
    typedef INTERN_RET::To_double_by_conversion< Real_embeddable >
                                                                      To_double;
    typedef INTERN_RET::To_interval_by_conversion< Real_embeddable >
                                                                    To_interval;
};

inline
io_Operator
io_tag(long int)
{ return io_Operator(); }

// short

template <> struct Number_type_traits<short int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;
};

template<> class Algebraic_structure_traits< short int >
  : public Algebraic_structure_traits_base< short int, 
                                            Euclidean_ring_tag > {

  public:
    typedef Tag_true            Is_exact;

    // Explicitly defined functors which have no support for implicit
    //  interoperability. This is nescessary because of the implicit conversion
    //  to int for binary operations between short ints.
    class Integral_division 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const { 
          Algebraic_structure_traits<Algebraic_structure>::Div actual_div;
          CGAL_precondition_msg( !is_exact(x) || actual_div( x, y) * y == x,
            "'x' must be divisible by 'y' in "
            "Algebraic_structure_traits<...>::Integral_div()(x,y)" );
          return actual_div( x, y);          
        }      
    };

    class Gcd 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const {
          Algebraic_structure_traits<Algebraic_structure>::Mod mod;
          Algebraic_structure_traits<Algebraic_structure>::Unit_part unit_part;
          Algebraic_structure_traits<Algebraic_structure>::Integral_division integral_div;
          // First: the extreme cases and negative sign corrections.
          if (x == Algebraic_structure(0)) {
              if (y == Algebraic_structure(0))  
                  return Algebraic_structure(0);
              return integral_div( y, unit_part(y) );
          }
          if (y == Algebraic_structure(0))
              return integral_div(x, unit_part(x) );
          Algebraic_structure u = integral_div( x, unit_part(x) );
          Algebraic_structure v = integral_div( y, unit_part(y) );
          // Second: assuming mod is the most expensive op here, we don't compute it
          // unnecessarily if u < v
          if (u < v) {
              v = mod(v,u);
              // maintain invariant of v > 0 for the loop below
              if ( v == Algebraic_structure(0) )
                  return u;
          }

          Algebraic_structure w;
          do {
              w = mod(u,v);
              if ( w == Algebraic_structure(0))
                  return v;
              u = mod(v,w);
              if ( u == Algebraic_structure(0))
                  return w;
              v = mod(w,u);
          } while (v != Algebraic_structure(0));
          return u;
        }        
    };

    class Div_mod { 
      public:
        typedef Algebraic_structure    first_argument_type;
        typedef Algebraic_structure    second_argument_type;
        typedef Algebraic_structure&   third_argument_type;
        typedef Algebraic_structure&   fourth_argument_type;
        typedef Arity_tag< 4 >         Arity;
        typedef void  result_type;
        void operator()( const Algebraic_structure& x, 
                         const Algebraic_structure& y, 
                         Algebraic_structure& q, Algebraic_structure& r) const {
          q = x / y;
          r = x % y;          
          return;
        }        
    };

    // based on \c Div_mod.
    class Div 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > {
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const {
          return x / y;
        };        
    };

    // based on \c Div_mod.
    class Mod 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const {
          return x % y;
        };        
    };
    
    typedef INTERN_INT::
       Is_square_per_double_conversion< Algebraic_structure > Is_square;
};

template <> class Real_embeddable_traits< short int > 
  : public Real_embeddable_traits_base< short int > {
  public:
          
    typedef INTERN_RET::To_double_by_conversion< Real_embeddable >
                                                                      To_double;
    typedef INTERN_RET::To_interval_by_conversion< Real_embeddable >
                                                                    To_interval;
};

inline
io_Operator
io_tag(short int)
{ return io_Operator(); }


// Note : "long long" support is in <CGAL/long_long.h>

CGAL_END_NAMESPACE

#endif // CGAL_INT_H
