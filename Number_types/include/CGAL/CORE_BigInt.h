// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Michael Hemmer  <mhemmer@uni-mainz.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Hemmer  <mhemmer@uni-mainz.de>
#ifndef CGAL_CORE_BIGINT_H
#define CGAL_CORE_BIGINT_H

#include <CGAL/number_type_basic.h>
#include <CGAL/CORE_coercion_traits.h>
#include <CGAL/CORE_Expr.h> // used for To_interval-functor

#define CORE_LEVEL 4
#include <CORE/CORE.h>

CGAL_BEGIN_NAMESPACE

//
// Algebraic structure traits
//
template <> class Algebraic_structure_traits< CORE::BigInt >
  : public Algebraic_structure_traits_base< CORE::BigInt, 
                                            Euclidean_ring_tag >  {
  public:
    typedef Tag_true            Is_exact;
                
    typedef INTERN_AST::Is_square_per_sqrt< Algebraic_structure >
                                                                 Is_square;            

    typedef INTERN_AST::Div_per_operator< Algebraic_structure > Div;
    typedef INTERN_AST::Mod_per_operator< Algebraic_structure > Mod;
    
    class Sqrt 
      : public Unary_function< Algebraic_structure, Algebraic_structure > {
      public:
        //! computes the largest NT not larger than the square root of \a a.
        Algebraic_structure operator()( const Algebraic_structure& x) const {
          Algebraic_structure result;
          mpz_sqrt(result.get_mp(), x.get_mp());
          return result;
        }
    };


    class Gcd 
      : public Binary_function< Algebraic_structure, Algebraic_structure,
                                Algebraic_structure > {
      public:        
        Algebraic_structure operator()( const Algebraic_structure& x,
                                        const Algebraic_structure& y) const {
          if ( x == Algebraic_structure(0) && y == Algebraic_structure(0) ) 
              return Algebraic_structure(0);
          Algebraic_structure result;
          mpz_gcd(result.get_mp(), x.get_mp(), y.get_mp());
          return result;
        }
    };
};

//
// Real embeddable traits
//
template <> class Real_embeddable_traits< CORE::BigInt > 
  : public Real_embeddable_traits_base< CORE::BigInt > {

  public:

    class Abs 
      : public Unary_function< Real_embeddable, Real_embeddable > {
      public:
        Real_embeddable operator()( const Real_embeddable& x ) const {
          return CORE::abs( x );
        }
    };
    
    class Sign 
      : public Unary_function< Real_embeddable, ::CGAL::Sign > {
      public:        
        ::CGAL::Sign operator()( const Real_embeddable& x ) const {
          return (::CGAL::Sign) CORE::sign( x );
        }        
    };
    
    class Compare 
      : public Binary_function< Real_embeddable, Real_embeddable,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Real_embeddable& x, 
                                            const Real_embeddable& y ) const {
          typedef Real_embeddable_traits<int> Int_traits;
          return Int_traits::Sign()( ::CORE::cmp(x,y));
        }
    };

    class To_double 
      : public Unary_function< Real_embeddable, double > {
      public:        
        double operator()( const Real_embeddable& x ) const {
          // this call is required to get reasonable values for the double
          // approximation 
          return x.doubleValue();
        }
    };
    
    class To_interval 
      : public Unary_function< Real_embeddable, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Real_embeddable& x ) const {
          Real_embeddable_traits<CORE::Expr>::To_interval to_interval;
          CORE::Expr temp(x);

          return to_interval(temp);
        }
    };
};

template<>
struct Needs_parens_as_product<CORE::BigInt>{
    bool operator()(const CORE::BigInt& x){
        return CGAL_NTS is_negative(x);
    }
};


CGAL_END_NAMESPACE

#endif // CGAL_CORE_BIGINT_H
