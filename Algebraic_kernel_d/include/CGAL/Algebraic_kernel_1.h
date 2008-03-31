// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//                 Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_KERNEL_1_H
#define CGAL_ALGEBRAIC_KERNEL_1_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_pure.h>
#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/Algebraic_kernel_d/Real_roots.h>
#include <CGAL/Algebraic_kernel_d/refine_zero_against.h>
#include <CGAL/Polynomial/ipower.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {
    template< class AlgebraicReal1, class Isolator_ >
    class Algebraic_kernel_1_base {
        public:
            
            typedef AlgebraicReal1              Algebraic_real_1;
            typedef Isolator_                   Isolator;
            
            typedef typename Algebraic_real_1::Coefficient      Coefficient;
            typedef typename Algebraic_real_1::Rational         Boundary;
           typedef Polynomial< Coefficient >                   Polynomial_1;
            
            class Algebraic_real_traits {
              public:
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
                                            CGAL::ipower( Boundary(2), rel_prec );
                            
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
            
            // Functors of Algebraic_kernel_1
            struct Solve_1 {
                template< class OutputIterator >
                OutputIterator operator()( const Polynomial_1& p,
                                           OutputIterator res,
                                           bool known_to_be_square_free = false ) const {
                    CGALi::Real_roots< Algebraic_real_1, Isolator > real_roots;

                    if( known_to_be_square_free ) {
                        real_roots( p, res );
                    } else {
                        std::vector< int > dummy;
                        real_roots( p, res, std::back_inserter( dummy ) );
                    }
                    
                    return res;
                }
                
                template< class OutputIteratorRoots, class OutputIteratorMults >
                std::pair< OutputIteratorRoots, OutputIteratorMults >
                operator()( const Polynomial_1& p, OutputIteratorRoots roots,
                            OutputIteratorMults mults ) const {
                    CGALi::Real_roots< Algebraic_real_1, Isolator > real_roots;
                    
                    real_roots( p, roots, mults );
                    
                    return std::pair< OutputIteratorRoots, OutputIteratorMults >(
                        roots, mults );                
                }                
            };
            
            struct Sign_at_1 
                : public Binary_function< Polynomial_1, Algebraic_real_1, CGAL::Sign > {
                    CGAL::Sign operator()( const Polynomial_1& p, const Algebraic_real_1& ar ) const {
                        typedef typename Algebraic_real_traits::Boundary Boundary;
                        Boundary low = ar.low();
                        Boundary high = ar.high();
                        if( low == high ) {
                            return p.sign_at( low );
                        }
                        
                        if( strong_refine_zero_against( low, high, ar.polynomial(), p ) ) {
                            return CGAL::ZERO;
                        } else {                                
                            CGAL_postcondition( p.sign_at(low) == p.sign_at(high) );                                                

                            return p.sign_at( low );
                        }                        
                    }
            };                
            
            struct Is_square_free_1 
                : public Unary_function< Polynomial_1, bool > {
                    bool operator()( const Polynomial_1& p ) const {
                        typename CGAL::Polynomial_traits_d< Polynomial_1 >::Is_square_free isf;
                        return isf(p);
                    }
            };
            
            struct Is_coprime_1
                : public Binary_function< Polynomial_1, Polynomial_1, bool > {
                    bool operator()( const Polynomial_1& p1, const Polynomial_1& p2 ) const {
                        typename CGAL::Polynomial_traits_d< Polynomial_1 >::Total_degree total_degree;
                        
                        // TODO: Is GCD already filtered? 
                        return( total_degree( gcd_utcf( p1, p2 ) ) == 0 );                        
                    } 
            }; 
            
            struct Make_square_free_1
                : public Unary_function< Polynomial_1, Polynomial_1 > {
                    Polynomial_1 operator()( const Polynomial_1& p ) const {
                        return typename CGAL::Polynomial_traits_d< Polynomial_1 >::Make_square_free()( p );
                    }
            };
            
            struct Make_coprime_1 {
                typedef bool         result_type;
                typedef Polynomial_1 first_argument_type;
                typedef Polynomial_1 second_argument_type;
                typedef Polynomial_1 third_argument_type;
                typedef Polynomial_1 fourth_argument_type;
                typedef Polynomial_1 fifth_argument_type;
                
                bool operator()( const Polynomial_1& p1,
                                 const Polynomial_1& p2,
                                 Polynomial_1& g, // ggT utcf 
                                 Polynomial_1& q1, // Rest utcf
                                 Polynomial_1& q2 ) const {
                    g = typename CGAL::Polynomial_traits_d< Polynomial_1 >::Gcd_up_to_constant_factor()( p1, p2 );
                    q1 = p1 / g;
                    q2 = p2 / g;
                    return Is_coprime_1()( q1, q2 );
                }                 
            };
            
            struct Square_free_factorization_1 {
                typedef int result_type;
                
                template< class OutputIterator1, class OutputIterator2 >
                int operator()( const Polynomial_1& p, OutputIterator1 it1,
                                OutputIterator2 it2 ) const {
                    return typename CGAL::Polynomial_traits_d< Polynomial_1 >::Square_free_factorization_up_to_constant_factor()( p, it1, it2 );
                } 
            };
            
    };
} // namespace CGALi


template< class Coefficient,
          class Boundary = typename CGAL::Get_arithmetic_kernel< Coefficient >::Arithmetic_kernel::Rational,
          class RepClass = CGALi::Algebraic_real_rep< Coefficient, Boundary >,
          class Isolator = CGALi::Descartes< typename CGAL::Polynomial< Coefficient >, Boundary > >
class Algebraic_kernel_1    
    : public CGALi::Algebraic_kernel_1_base< 

    // Template argument #1 (AlgebraicReal1)        
        CGALi::Algebraic_real_pure< 
            Coefficient, 
            Boundary,
            ::CGAL::Handle_policy_no_union,     
            RepClass >,
        
    // Template argument #2 (Isolator_)
        Isolator >
            
     {};


CGAL_END_NAMESPACE


#endif // CGAL_ALGEBRAIC_KERNEL_1_H
