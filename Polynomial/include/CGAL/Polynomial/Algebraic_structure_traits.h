// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


#ifndef CGAL_POLYNOMIAL_ALGEBRAIC_STRUCTURE_TRAITS_H
#define CGAL_POLYNOMIAL_ALGEBRAIC_STRUCTURE_TRAITS_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

// Extend to a UFDomain as coefficient range
// Forward declaration for <NiX/polynomial_gcd.h> for NT_traits<Poly...>::Gcd
namespace CGALi {
template <class NT> inline
Polynomial<NT> gcd(const Polynomial<NT>&, const Polynomial<NT>&);
} // namespace CGALi


// Now we wrap up all of this in the actual NT_traits
// specialization for Polynomial<NT>
/*! \ingroup NiX_Polynomial
    \brief \c NiX::NT_traits < \c NiX::Polynomial<NT> >
 *
 *  If \c NT is a model of a number type concept, then so is
 *  \c Polynomial<NT>. A specialization of \c NiX::NT_traits
 *  is provided automatically by NiX/Polynomial.h.
 *
 *  The number type concepts for the coefficient domain NT are
 *  mapped to those for the polynomials as follows:
 *  <PRE>
    IntegralDomainWithoutDiv --> IntegralDomainWithoutDiv
    IntegralDomain           --> IntegralDomain
    UFDomain       --> UFDomain
    EuclideanRing  --> UFDomain
    Field          --> EuclideanRing
    FieldWithSqrt  --> EuclideanRing
    </PRE>
 *
 *  \c Polynomial<NT> is \c RealComparable iff \c NT is.
 *  The ordering is determined by the \c sign() of differences, see ibid.
 *
 *  \c Polynomial<NT> offers a non-<TT>Null_tag</TT> \c To_double
 *  iff \c NT does. If non-<TT>Null_tag</TT>, it returns a coefficient-wise
 *  \c double approximation of the polynomial.
 */



// Algebraic structure traits

template< class POLY, class Algebraic_type >
class Polynomial_algebraic_structure_traits_base;

// The most basic suite of algebraic operations, suitable for the
// most basic kind of coefficient range, viz. a IntegralDomainWithoutDiv.
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, 
                                             Integral_domain_without_division_tag >
  : public Algebraic_structure_traits_base< POLY, 
                                            Integral_domain_without_division_tag > {
  public:
    typedef Integral_domain_without_division_tag Algebraic_category;
    
    class Simplify 
      : public Unary_function< POLY&, void > {
      public:
        void operator()( POLY& p ) const {
          p.simplify_coefficients();
        }        
    };
    
    class Unit_part 
      : public Unary_function< POLY, POLY > {
      public:
        POLY operator()( const POLY& x ) const {
          return POLY( x.unit_part() );
        }
    };
};

// Extend to the case that the coefficient range is a IntegralDomain (with div)
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Integral_domain_tag >
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                            Integral_domain_without_division_tag > {
  public:
    typedef Integral_domain_tag Algebraic_category;
    
    class Integral_division 
      : public Binary_function< POLY, POLY, POLY > {
      public:
        POLY operator()( const POLY& x, const POLY& y ) const {
          return x / y;
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
        
    };
};

template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Unique_factorization_domain_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                            Integral_domain_tag > {
  public:
    typedef Unique_factorization_domain_tag Algebraic_category;
    
    class Gcd 
      : public Binary_function< POLY, POLY, POLY > {
      public:
        POLY operator()( const POLY& x, const POLY& y ) const {
            typedef Algebraic_structure_traits<POLY> AST;
            typename AST::Integral_division idiv;
            typename AST::Unit_part upart; 

          // First: the extreme cases and negative sign corrections.          
          if (x == POLY(0)) {
              if (y == POLY(0))  
                  return POLY(0);
              return idiv(y,upart(y));
          }
          if (y == POLY(0))
              return idiv(x,upart(x));
          
          return CGALi::gcd(x,y);
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
        
    };
};

// Clone this for a EuclideanRing
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Euclidean_ring_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                            Unique_factorization_domain_tag > {
  // nothing new
};

// Extend to a field as coefficient range
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Field_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                            Unique_factorization_domain_tag > {
  public:
    typedef Euclidean_ring_tag Algebraic_category;
    
    class Div_mod {
      public:
        typedef POLY first_argument_type;
        typedef POLY second_argument_type;
        typedef POLY& third_argument_type;
        typedef POLY& fourth_argument_type;
        typedef void result_type;
        typedef Arity_tag< 4 >  Arity;
        
        void operator()( const POLY& x, const POLY& y, 
                         POLY& q, POLY& r ) const {
          POLY::euclidean_division( x, y, q, r );
        }
        
        template < class NT1, class NT2 >
        void operator()( const NT1& x, const NT2& y,
                         POLY& q, POLY& r ) const {
          BOOST_STATIC_ASSERT((::boost::is_same<
                  typename Coercion_traits< NT1, NT2 >::Type, POLY
                                               >::value));
          
          typename Coercion_traits< NT1, NT2 >::Cast cast;
          operator()( cast(x), cast(y), q, r );          
        }
        
    };
    
    class Div 
      : public Binary_function< POLY, POLY, POLY > {
      public:
        POLY operator()(const POLY& a, const POLY& b) const {
          POLY q, r;
          Div_mod()(a, b, q, r);
          return q;
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
    };
    
    class Mod 
      : public Binary_function< POLY, POLY, POLY > {
      public:
        POLY operator () (const POLY& a, const POLY& b) const {
          POLY q, r;
          Div_mod()(a, b, q, r);
          return r;
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
    };
    
};

// Clone this for a FieldWithSqrt
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Field_with_sqrt_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, Field_tag > {
  // nothing new
};

// Clone this for a FieldWithKthRoot
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, 
                                                  Field_with_kth_root_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                                       Field_with_sqrt_tag > {
  // nothing new
};

// Clone this for a FieldWithRootOf
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, 
                                                  Field_with_root_of_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                                    Field_with_kth_root_tag > {
  // nothing new
};

// The actual algebraic structure traits class
template< class NT > class Algebraic_structure_traits< Polynomial< NT > >
  : public Polynomial_algebraic_structure_traits_base< Polynomial< NT >,
         typename Algebraic_structure_traits< NT >::Algebraic_category > {
  public:
    typedef Polynomial<NT> Algebraic_structure;
    typedef typename Algebraic_structure_traits< NT >::Is_exact Is_exact;
};

CGAL_END_NAMESPACE
#endif // CGAL_POLYNOMIAL_ALGEBRAIC_STRUCTURE_TRAITS_H
