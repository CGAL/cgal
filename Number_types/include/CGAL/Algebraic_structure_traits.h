// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
// 
//
// Author(s)     : Michael Hemmer  <mhemmer@uni-mainz.de>

#ifndef CGAL_ALGEBRAIC_STRUCTURE_TRAITS_H
#define CGAL_ALGEBRAIC_STRUCTURE_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/functional_base.h> // Unary_function, Binary_function
#include <boost/type_traits/is_same.hpp>


CGAL_BEGIN_NAMESPACE

// REMARK: Some of the following comments and references are just copy & pasted
//         from EXACUS and have to be adapted/removed in the future.

// The tags for Algebra_type corresponding to the number type concepts
// ===================================================================

//! corresponds to the \c IntegralDomainWithoutDiv concept.
struct Integral_domain_without_division_tag {};

//! corresponds to the \c IntegralDomain concept.
struct Integral_domain_tag : public Integral_domain_without_division_tag {};

//! corresponds to the \c UFDomain concept.
struct Unique_factorization_domain_tag : public Integral_domain_tag {};

//! corresponds to the \c EuclideanRing concept.
struct Euclidean_ring_tag : public Unique_factorization_domain_tag {};

//! corresponds to the \c Field concept.
struct Field_tag : public Integral_domain_tag {};

//! corresponds to the \c FieldWithSqrt concept.
struct Field_with_sqrt_tag : public Field_tag {};

//! corresponds to the \c FieldWithKthRoot concept
struct Field_with_kth_root_tag : public Field_with_sqrt_tag {};

//! corresponds to the \c FieldWithRootOF concept.
struct Field_with_root_of_tag : public Field_with_kth_root_tag {};


// The algebraic structure traits template
// =========================================================================
template< class Algebraic_structure_ > 
class Algebraic_structure_traits  {
  public:
    typedef Algebraic_structure_  Algebraic_structure;
    typedef Null_tag       Algebraic_structure_tag;
    typedef Null_tag       Is_exact;

    typedef Null_functor Simplify;
    typedef Null_functor Unit_part;
    typedef Null_functor Integral_division;
    typedef Null_functor Is_square;    
    typedef Null_functor Gcd;
    typedef Null_functor Div_mod;
    typedef Null_functor Div;
    typedef Null_functor Mod;
    typedef Null_functor Square;
    typedef Null_functor Is_zero;
    typedef Null_functor Is_one;    
    typedef Null_functor Sqrt;
    typedef Null_functor Kth_root;
    typedef Null_functor Root_of; 
  
};

// Do we need this function? 
// Currently it causes an error in package Min_sphere_of_spheres
template< class Algebraic_structure >
const bool is_exact( const Algebraic_structure& as ) {
    return check_tag( 
            typename Algebraic_structure_traits< Algebraic_structure >::Is_exact() );
}

// The algebraic structure traits base class
// =========================================================================
template< class Algebraic_structure, class Algebra_type >
class Algebraic_structure_traits_base;

//! The template specialization that can be used for types that are not any
//! of the number type concepts. All functors are set to \c Null_functor
//! or suitable defaults. The \c Simplify functor does nothing by default.
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, Null_tag > {
  public:
    typedef Algebraic_structure_  Algebraic_structure;
    typedef Null_tag       Algebraic_structure_tag;
    typedef Tag_false       Is_exact;

    // does nothing by default
    class Simplify 
      : public Unary_function< Algebraic_structure&, void > {
      public:
        void operator()( Algebraic_structure& ) const {}
    };

    typedef Null_functor Unit_part;
    typedef Null_functor Integral_division;
    typedef Null_functor Is_square;    
    typedef Null_functor Gcd;
    typedef Null_functor Div_mod;
    typedef Null_functor Div;
    typedef Null_functor Mod;
    typedef Null_functor Square;
    typedef Null_functor Is_zero;
    typedef Null_functor Is_one;    
    typedef Null_functor Sqrt;
    typedef Null_functor Kth_root;
    typedef Null_functor Root_of; 
};

//! The template specialization that is used if the number type is
//! a model of the \c IntegralDomainWithoutDiv concept. The \c Simplify
//! does nothing by default and the \c Unit_part is equal to
//! \c Algebraic_structure(-1) for negative numbers and 
//! \c Algebraic_structure(1) otherwise
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, 
                                       Integral_domain_without_division_tag > 
    : public Algebraic_structure_traits_base< Algebraic_structure_, 
                                              Null_tag > {
  public:
    typedef Algebraic_structure_                   Algebraic_structure;
    typedef Integral_domain_without_division_tag  Algebraic_structure_tag;

    // returns Algebraic_structure(1) by default
    class Unit_part 
      : public Unary_function< Algebraic_structure, Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x ) const {
          return( x < Algebraic_structure(0)) ? 
                  Algebraic_structure(-1) : Algebraic_structure(1); 
        }
    };
    
    class Square 
      : public Unary_function< Algebraic_structure, Algebraic_structure > {
      public:        
        Algebraic_structure operator()( const Algebraic_structure& x ) const {
          return x*x;
        }
    };
    
    class Is_zero 
      : public Unary_function< Algebraic_structure, bool > {
      public:        
        bool operator()( const Algebraic_structure& x ) const {
          return x == Algebraic_structure(0);
        }
    };

    class Is_one 
      : public Unary_function< Algebraic_structure, bool > {
      public:        
        bool operator()( const Algebraic_structure& x ) const {
          return x == Algebraic_structure(1);
        }
    };

};


//! The template specialization that is used if the number type is
//! a model of the \c IntegralDomain concept. It is equivalent to the 
//! specialization
//! for the \c IntegralDomainWithoutDiv concept. The additionally required 
//! \c Integral_division functor needs to be implemented in the 
//! \c Algebraic_structure_traits itself.
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, 
                                       Integral_domain_tag >
    : public Algebraic_structure_traits_base< Algebraic_structure_, 
                                       Integral_domain_without_division_tag > {
  public:
    typedef Algebraic_structure_       Algebraic_structure;
    typedef Integral_domain_tag  Algebraic_structure_tag;
};


//! The template specialization that is used if the number type is
//! a model of the \c UFDomain concept. It is equivalent to the specialization
//! for the \c IntegralDomain concept. The additionally required 
//! \c Integral_div functor
//! and \c Gcd functor need to be implemented in the 
//! \c Algebraic_structure_traits itself.
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, 
                                       Unique_factorization_domain_tag >
    : public Algebraic_structure_traits_base< Algebraic_structure_, 
                                              Integral_domain_tag > {
  public:
    typedef Algebraic_structure_  Algebraic_structure;
    typedef Unique_factorization_domain_tag    Algebraic_structure_tag;
};


//! The template specialization that is used if the number type is
//! a model of the \c EuclideanRing concept.
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, 
                                       Euclidean_ring_tag >
    : public Algebraic_structure_traits_base< Algebraic_structure_, 
                                              Unique_factorization_domain_tag > {
  public:
    typedef Algebraic_structure_        Algebraic_structure;
    typedef Euclidean_ring_tag    Algebraic_structure_tag;

    // maps to \c Div by default.
    class Integral_division 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const { 
          typename Algebraic_structure_traits<Algebraic_structure>::Div 
                                                                    actual_div;
          CGAL_precondition_msg( !is_exact(x) || actual_div( x, y) * y == x,
            "'x' must be divisible by 'y' in "
            "Algebraic_structure_traits<...>::Integral_div()(x,y)" );
          return actual_div( x, y);          
        }
      
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )  
    };

    // Algorithm from NiX/euclids_algorithm.h
    class Gcd 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const {
          typename Algebraic_structure_traits<Algebraic_structure>::Mod mod;
          typename Algebraic_structure_traits<Algebraic_structure>::Unit_part unit_part;
          typename Algebraic_structure_traits<Algebraic_structure>::Integral_division integral_div;
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
          // Third: generic case of two positive integer values and u >= v.
          // The standard loop would be:
          //      while ( v != 0) {
          //          int tmp = mod(u,v);
          //          u = v;
          //          v = tmp;
          //      }
          //      return u;
          //
          // But we want to save us all the variable assignments and unroll
          // the loop. Before that, we transform it into a do {...} while()
          // loop to reduce branching statements.
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
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };

    // based on \c Div and \c Mod.
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
          typedef Algebraic_structure_traits<Algebraic_structure> Traits;
          typename Traits::Div  actual_div;
          typename Traits::Mod  actual_mod;
          q = actual_div( x, y );
          r = actual_mod( x, y );          
          return;
        }
        
        template < class NT1, class NT2 >
        void operator()( const NT1& x, const NT2& y,
                         Algebraic_structure& q, Algebraic_structure& r ) const {
          BOOST_STATIC_ASSERT((::boost::is_same<
                  typename Coercion_traits< NT1, NT2 >::Coercion_type, Algebraic_structure
                                               >::value));
          
          typename Coercion_traits< NT1, NT2 >::Cast cast;
          operator()( cast(x), cast(y), q, r );          
        }
    };

    // based on \c Div_mod.
    class Div 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > {
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const {
          typename Algebraic_structure_traits<Algebraic_structure>
                                                    ::Div_mod actual_div_mod;
          Algebraic_structure q;     
          Algebraic_structure r;
          actual_div_mod( x, y, q, r );
          return q;
        };
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };

    // based on \c Div_mod.
    class Mod 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const {
          typename Algebraic_structure_traits<Algebraic_structure>
                                                    ::Div_mod actual_div_mod;
          Algebraic_structure q;     
          Algebraic_structure r;
          actual_div_mod( x, y, q, r );
          return r;
        };
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };
};


//! The template specialization that is used if the number type is
//! a model of the \c Field concept. \c Unit_part ()(x)
//! returns \c NT(1) if the value \c x is equal to \c NT(0) and
//! otherwise the value \c x itself. The \c Integral_div
//! maps to the \c operator/.
//! See also \link NiX_NT_traits_functors concept NT_traits \endlink .
//! \ingroup NiX_NT_traits_bases
//
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, Field_tag >
    : public Algebraic_structure_traits_base< Algebraic_structure_, 
                                              Integral_domain_tag > {
  public:
    typedef Algebraic_structure_        Algebraic_structure;
    typedef Field_tag             Algebraic_structure_tag;

    // returns the argument \a a by default
    class Unit_part 
      : public Unary_function< Algebraic_structure, Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x ) const {
            return( x == Algebraic_structure(0)) ? Algebraic_structure(1) : x;
        }
    };
    // maps to \c operator/ by default.
    class Integral_division 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > { 
      public:
        Algebraic_structure operator()( const Algebraic_structure& x, 
                                        const Algebraic_structure& y) const { 

          CGAL_precondition_msg( !is_exact(x) || (x / y) * y  == x,
            "'x' must be divisible by 'y' in "
            "Algebraic_structure_traits<...>::Integral_div()(x,y)" );
          return x / y;
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
    };
};


//! The template specialization that is used if the number type is a model
//! of the \c FieldWithSqrt concept. It is equivalent to the 
//! specialization for the \c Field concept. The additionally required 
//! \c NiX::NT_traits::Sqrt functor need to be
//! implemented in the \c NT_traits itself.
//! \ingroup NiX_NT_traits_bases
//
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, 
                                       Field_with_sqrt_tag>
    : public Algebraic_structure_traits_base< Algebraic_structure_, 
                                              Field_tag> {
  public:
    typedef Algebraic_structure_        Algebraic_structure;
    typedef Field_with_sqrt_tag   Algebraic_structure_tag;

    struct Is_square
        :public Binary_function<Algebraic_structure,Algebraic_structure&,bool>
    {
        bool operator()(const Algebraic_structure& x) const {return true;}
        bool operator()(
                const Algebraic_structure& x,
                Algebraic_structure      & result) const {
            typename Algebraic_structure_traits<Algebraic_structure>::Sqrt sqrt;
            result = sqrt(x);
            return true;
        }
    };
};

//! The template specialization that is used if the number type is a model
//! of the \c FieldWithKthRoot concept. It is equivalent to the 
//! specialization for the \c Field concept. The additionally required 
//! \c NiX::NT_traits::Kth_root functor need to be
//! implemented in the \c Algebraic_structure_traits itself.
//! \ingroup NiX_NT_traits_bases
//
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, 
                                       Field_with_kth_root_tag>
    : public Algebraic_structure_traits_base< Algebraic_structure_, 
                                              Field_with_sqrt_tag> {
    
    
    
  public:
    typedef Algebraic_structure_        Algebraic_structure;
    typedef Field_with_kth_root_tag   Algebraic_structure_tag;
};




//! The template specialization that is used if the number type is a model
//! of the \c FieldWithRootOf concept. It is equivalent to the 
//! specialization for the \c FieldWithKthRoot concept. The additionally 
//! required \c NiX::NT_traits::Root_of functor need to be
//! implemented in the \c NT_traits itself.
//! \ingroup NiX_NT_traits_bases
//
template< class Algebraic_structure_ >
class Algebraic_structure_traits_base< Algebraic_structure_, 
                                       Field_with_root_of_tag >
    : public Algebraic_structure_traits_base< Algebraic_structure_, 
                                              Field_with_kth_root_tag > {
  public:
    typedef Algebraic_structure_           Algebraic_structure;
    typedef Field_with_root_of_tag   Algebraic_structure_tag;
};

// Some common functors to be used by AST specializations
namespace INTERN_AST {
  template< class Algebraic_structure >
  class Div_per_operator 
    : public Binary_function< Algebraic_structure, Algebraic_structure, 
                              Algebraic_structure > {
    public:      
      Algebraic_structure operator()( const Algebraic_structure& x, 
                                      const Algebraic_structure& y ) const {
        return x / y;
      }
      
      CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
  };
  
  template< class Algebraic_structure >
  class Mod_per_operator 
    : public Binary_function< Algebraic_structure, Algebraic_structure,
                              Algebraic_structure > {
    public:
      Algebraic_structure operator()( const Algebraic_structure& x, 
                                      const Algebraic_structure& y ) const {
        return x % y;
      }
      
      CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
  };
  
  template< class Algebraic_structure >
  class Is_square_per_sqrt
    : public Binary_function< Algebraic_structure, Algebraic_structure&,
                              bool > {
    public:      
      bool operator()( const Algebraic_structure& x, 
                       Algebraic_structure& y ) const {
          typename Algebraic_structure_traits< Algebraic_structure >::Sqrt
              actual_sqrt;
          y = actual_sqrt( x );
          return y * y == x;
      }
      bool operator()( const Algebraic_structure& x) const {
          Algebraic_structure dummy;
          return operator()(x,dummy);
      }
  };
} // INTERN_AST






CGAL_END_NAMESPACE

#endif  // CGAL_ALGEBRAIC_STRUCTURE_TRAITS_H
