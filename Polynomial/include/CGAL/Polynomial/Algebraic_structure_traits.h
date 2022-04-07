// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
#include <CGAL/Coercion_traits.h>
#include <CGAL/Polynomial/modular_filter.h>

namespace CGAL {

// Extend to a UFDomain as coefficient range
// Forward declaration for <NiX/polynomial_gcd.h> for NT_traits<Poly...>::Gcd
namespace internal {
template <class NT> inline
Polynomial<NT> gcd_(const Polynomial<NT>&, const Polynomial<NT>&);
} // namespace internal


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
      : public CGAL::cpp98::unary_function< POLY&, void > {
      public:
        void operator()( POLY& p ) const {
          p.simplify_coefficients();
        }
    };

    class Unit_part
      : public CGAL::cpp98::unary_function< POLY, POLY > {
      public:
        POLY operator()( const POLY& x ) const {
          return POLY( x.unit_part() );
        }
    };

  class Is_zero
    : public CGAL::cpp98::unary_function< POLY, bool > {
  public:
    bool operator()( const POLY& x ) const {
      return x.is_zero();
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
      : public CGAL::cpp98::binary_function< POLY, POLY, POLY > {
      public:
        POLY operator()( const POLY& x, const POLY& y ) const {
          return x / y;
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
    };
private:
  typedef typename POLY::NT COEFF;
  typedef Algebraic_structure_traits<COEFF> AST_COEFF;
  typedef typename AST_COEFF::Divides Divides_coeff;
  typedef typename Divides_coeff::result_type BOOL;
public:
  class Divides
    : public CGAL::cpp98::binary_function<POLY,POLY,BOOL>{
  public:
     BOOL operator()( const POLY& p1, const POLY& p2) const {
       POLY q;
       return (*this)(p1,p2,q);
     }
     BOOL operator()( const POLY& p1, const POLY& p2, POLY& q) const {
       Divides_coeff divides_coeff;

       q=POLY(0);

       COEFF q1;
       BOOL result;
       if (p2.is_zero()) {
         q=POLY(0);
         return true;
       }

       int d1 = p1.degree();
       int d2 = p2.degree();
       if ( d2 < d1 ) {
         q = POLY(0);
         return false;
       }

       typedef std::vector<COEFF> Vector;
       Vector V_R, V_Q;
       V_Q.reserve(d2);
       if(d1==0){
         for(int i=d2;i>=0;--i){
           result=divides_coeff(p1[0],p2[i],q1);
           if(!result) return false;
           V_Q.push_back(q1);
         }
         V_R.push_back(COEFF(0));
       }
       else{
         V_R.reserve(d2);
         V_R=Vector(p2.begin(),p2.end());
         Vector tmp1;
         tmp1.reserve(d1);
         for(int k=0; k<=d2-d1; ++k){
           result=divides_coeff(p1[d1],V_R[d2-k],q1);
           if(!result) return false;
           V_Q.push_back(q1);
           for(int j=0;j<d1;++j){
             tmp1.push_back(p1[j]*V_Q[k]);
           }
           V_R[d2-k]=COEFF(0);
           for(int i=d2-d1-k;i<=d2-k-1;++i){
             V_R[i]=V_R[i]-tmp1[i-(d2-d1-k)];
           }
           tmp1.clear();
         }


       }
       q = POLY(V_Q.rbegin(),V_Q.rend());
       POLY r = POLY(V_R.begin(),V_R.end());

       return (r == POLY(0));
     }
     CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(POLY,BOOL)
   };
};

template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Unique_factorization_domain_tag >
  : public Polynomial_algebraic_structure_traits_base< POLY,
                                            Integral_domain_tag > {
  public:
    typedef Unique_factorization_domain_tag Algebraic_category;

  class Gcd
    : public CGAL::cpp98::binary_function< POLY, POLY, POLY > {
    typedef typename Polynomial_traits_d<POLY>::Multivariate_content Mcontent;
    typedef typename Mcontent::result_type ICoeff;

    ICoeff gcd_help(const ICoeff& , const ICoeff& , Field_tag) const {
      return ICoeff(1);
    }
    ICoeff gcd_help(const ICoeff& x, const ICoeff& y,
        Unique_factorization_domain_tag) const {
      return CGAL::gcd(x,y);
    }
  public:
    POLY operator()( const POLY& x, const POLY& y ) const {
      if(x==y) return x;

      typedef Algebraic_structure_traits<POLY> AST;
      typename AST::Integral_division idiv;
      typename AST::Unit_part upart;

      // First: the extreme cases and negative sign corrections.
      if (CGAL::is_zero(x)) {
        if (CGAL::is_zero(y))
          return POLY(0);
        return idiv(y,upart(y));
      }
      if (CGAL::is_zero(y))
        return idiv(x,upart(x));

      if (internal::may_have_common_factor(x,y))
        return CGAL::internal::gcd_(x,y);
      else{
        typename Algebraic_structure_traits<ICoeff>::Algebraic_category category;
        return POLY(gcd_help(Mcontent()(x),Mcontent()(y), category));
      }
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

        void operator()( const POLY& x, const POLY& y,
                         POLY& q, POLY& r ) const {
          POLY::euclidean_division( x, y, q, r );
        }

        template < class NT1, class NT2 >
        void operator()( const NT1& x, const NT2& y,
                         POLY& q, POLY& r ) const {
          CGAL_static_assertion((::boost::is_same<
                  typename Coercion_traits< NT1, NT2 >::Type, POLY
                                               >::value));

          typename Coercion_traits< NT1, NT2 >::Cast cast;
          operator()( cast(x), cast(y), q, r );
        }

    };

    class Div
      : public CGAL::cpp98::binary_function< POLY, POLY, POLY > {
      public:
        POLY operator()(const POLY& a, const POLY& b) const {
          POLY q, r;
          Div_mod()(a, b, q, r);
          return q;
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
    };

    class Mod
      : public CGAL::cpp98::binary_function< POLY, POLY, POLY > {
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

} //namespace CGAL
#endif // CGAL_POLYNOMIAL_ALGEBRAIC_STRUCTURE_TRAITS_H
