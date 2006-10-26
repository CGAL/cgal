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
//                 Michael Hemmer  <mhemmer@uni-mainz.de>

#ifndef CGAL_CORE_BIGFLOAT_H
#define CGAL_CORE_BIGFLOAT_H

#include <CGAL/basic.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/utils.h>
#include <CGAL/CORE_coercion_traits.h>
#include <CGAL/functional_base.h> // Unary_function, Binary_function


#define CORE_LEVEL 4
#include <CORE/CORE.h>

CGAL_BEGIN_NAMESPACE

//todo: rm
template <> struct Number_type_traits<CORE::BigFloat> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_false Has_exact_division;
  typedef Tag_false Has_exact_sqrt;
};

//
// Algebraic structure traits
//
template <> class Algebraic_structure_traits< CORE::BigFloat >
  : public Algebraic_structure_traits_base< CORE::BigFloat, 
                                            CGAL::Field_with_kth_root_tag >  {
  public:
    typedef CGAL::Tag_false            Is_exact;
                          
    class Sqrt 
      : public Unary_function< Algebraic_structure, Algebraic_structure > {
      public:        
        Algebraic_structure operator()( const Algebraic_structure& x ) const {
          return CORE::sqrt( x );
        }
    };
    
    class Kth_root 
      : public Binary_function<int, Algebraic_structure, Algebraic_structure> {
      public:
        Algebraic_structure operator()( int k, 
                                        const Algebraic_structure& x) const {
            CGAL_precondition_msg( k > 0, "'k' must be positive for k-th roots");
            // CORE::radical isn't implemented for negative values of x, so we
            //  have to handle this case separately
            if( x < 0 && k%2 != 0) {
              return Algebraic_structure(-CORE::radical( -x, k ) );
            }
  
  
            return Algebraic_structure( CORE::radical( x, k ) );
        }
    };    
};

//
// Real embeddable traits
//
template <> class Real_embeddable_traits< CORE::BigFloat > 
  : public Real_embeddable_traits_base< CORE::BigFloat > {
  public:
      
    class Abs 
      : public Unary_function< Real_embeddable, Real_embeddable > {
      public:
        Real_embeddable operator()( const Real_embeddable& x ) const {
            return CORE::abs( x );
        }
    };
    
    class Sign 
      : public Unary_function< Real_embeddable, CGAL::Sign > {
      public:        
        CGAL::Sign operator()( const Real_embeddable& x ) const {
          return (CGAL::Sign) CORE::sign( x );
        }        
    };
    
    class Compare 
      : public Binary_function< Real_embeddable, Real_embeddable,
                                CGAL::Comparison_result > {
      public:        
        CGAL::Comparison_result operator()( const Real_embeddable& x, 
                                            const Real_embeddable& y ) const {
          return (CGAL::Comparison_result) CORE::cmp( x, y );
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


CGAL_END_NAMESPACE

#endif // CGAL_CORE_BIGFLOAT_H
