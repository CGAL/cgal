// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_KERNEL_1_H
#define CGAL_ALGEBRAIC_KERNEL_1_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_real_pure.h>
#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/Polynomial/ipower.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {
    template< class AlgebraicReal1, class Isolator_ >
    class Algebraic_kernel_d_1_base {
        public:
            
            typedef AlgebraicReal1              Algebraic_real_1;
            typedef Isolator_                   Isolator;
            
            typedef typename Algebraic_real_1::Coefficient      Coefficient;
            typedef Polynomial< Coefficient > Polynomial_1;
            
            class Algebraic_real_traits {
                typedef Algebraic_real_1                      Type;
                typedef typename Algebraic_real_1::Rational   Boundary;
                
                struct Boundary_between 
                    : public Binary_function< Type, Type, Boundary > {
                    Boundary operator()( const Type& t1, 
                                         const Type& t2 ) const {
                        return t1.rational_between( t2 );
                    }
                };
                                
                struct Lower_boundary
                    : public Unary_function< Type, Boundary > {
                    Boundary operator()( const Type& t ) const {
                        return t.low();
                    }
                };
                
                struct Upper_boundary
                    : public Unary_function< Type, Boundary > {
                    Boundary operator()( const Type& t ) const {
                        return t.high();
                    }
                };
                
                struct Refine
                    : public Unary_function< Type, void > {
                    void operator()( const Type& t ) const {
                        t.refine();
                    }
                    
                    void operator()( Type& t, int rel_prec ) const {
                        // If t is zero, we can refine the interval to
                        //  infinite precission
                        if( CGAL::is_zero( t ) ) {
                            t = Type(0);                            
                        } else {
                            // Refine until both boundaries have the same sign
                            while( CGAL::sign( t.high() ) != 
                                   CGAL::sign( t.low() ) )
                                t.refine();
                            
                            CGAL_assertion( CGAL::sign( t.high() ) != CGAL::ZERO &&
                                            CGAL::sign( t.low() ) != CGAL::ZERO );
                            
                            // Calculate the needed precision
                            Boundary prec = Boundary(1) / 
                                            INTERN_POLYNOMIAL::ipower( Boundary(2), rel_prec );
                            
                            // Refine until precision is reached
                            while( CGAL::abs( t.high() - t.low() ) /
                                   CGAL::max( CGAL::abs( t.high() ),
                                              CGAL::abs( t.low() ) ) > prec ) {
                                t.refine();

                                CGAL_assertion( CGAL::sign( t.high() ) != CGAL::ZERO &&
                                                CGAL::sign( t.low() ) != CGAL::ZERO );
                                
                            } 
                        }                        
                    }
                };                
            }; // class Algebraic_real_traits
    };
} // namespace CGALi

#define CGAL_AK_BOUNDARY    Rational
#define CGAL_AK_REPCLASS    CGALi::Algebraic_real_rep
#define CGAL_AK_ISOLATOR    CGALi::Descartes


template< class Coefficient >
class Algebraic_kernel_d_1    
    : public CGALi::Algebraic_kernel_d_1_base< 

    // Template argument #1 (AlgebraicReal1)        
        CGALi::Algebraic_real_pure< 
            Coefficient, 
            typename CGALi::Get_arithmetic_kernel< Coefficient >::CGAL_AK_BOUNDARY,
            ::CGAL::Handle_policy_no_union,     
            CGAL_AK_REPCLASS< Coefficient, typename CGALi::Get_arithmetic_kernel< Coefficient >::CGAL_AK_BOUNDARY > >,
        
    // Template argument #2 (Isolator_)
        CGAL_AK_ISOLATOR< typename CGAL::Polynomial< Coefficient >,
                          typename CGALi::Get_arithmetic_kernel< Coefficient >::CGAL_AK_BOUNDARY > >
            
     {};


CGAL_END_NAMESPACE


#endif // CGAL_ALGEBRAIC_KERNEL_1_H
