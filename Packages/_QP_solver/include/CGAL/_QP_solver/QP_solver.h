// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/_QP_solver/QP_solver.h
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.6
// revision_date : 2000/09/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Solver for Quadratic Programs
// ============================================================================
                                                                               
#ifndef CGAL_QP_SOLVER_H
#define CGAL_QP_SOLVER_H

// includes
#include <CGAL/Optimisation/basic.h>

#include <CGAL/_QP_solver/Basis_inverse.h>
//#include <CGAL/_QP_solver/Pricing_strategy_base.h>
//#include <CGAL/_QP_solver/Full_exact_pricing.h>

#include <CGAL/_QP_solver/Join_random_access_iterator.h>
#include <CGAL/_QP_solver/Sparse_vector_iterator.h>

#include <CGAL/_QP_solver/Compute_quotient.h>
#include <CGAL/_QP_solver/Access_by_index.h>

#include <CGAL/IO/Verbose_ostream.h>

#include <vector>
#include <functional>
#include <algorithm>
#include <iterator>


CGAL_BEGIN_NAMESPACE
                    

// Class declaration
// =================
template < class Rep_ >
class QP_solver;
                
template < class Rep_ >
class Pricing_strategy_base;
                
template < class Rep_ >
class Full_exact_pricing;
                

// Class interface
// ===============
template < class Rep_ >
class QP_solver {
  public:
    // self
    typedef  Rep_                       Rep;
    typedef  QP_solver<Rep>             Self;

    // types from the representation class
    typedef  typename Rep::NT           NT;
    typedef  typename Rep::ET           ET;

    typedef  typename Rep::A_iterator   A_iterator;
    typedef  typename Rep::B_iterator   B_iterator;
    typedef  typename Rep::C_iterator   C_iterator;
    typedef  typename Rep::D_iterator   D_iterator;

    typedef  typename Rep::Is_lp        Is_lp;

  private:
    // private types
    typedef  CGAL::Tag_true             Tag_true;
    typedef  CGAL::Tag_false            Tag_false;
    
    
    typedef  std::vector<int>           Indices;
    typedef  Indices::const_iterator    Index_iterator;
    
    typedef  std::vector<ET>            Values;
    typedef  typename Values::const_iterator
                                        Value_iterator;
    
    
    typedef  CGAL::Compute_quotient<ET> Compute_quotient;
    typedef  std::binder2nd< Compute_quotient >
                                        Make_quotient;
    
    
    typedef  CGAL::Access_by_index<Value_iterator,true,false>
                                        Access_value_checked;
    
    
    typedef  CGAL::Sparse_vector_iterator<NT>
                                        Artificial_column;
    
    typedef  std::vector<Artificial_column>
                                        A_artificial;
    
    
    typedef  std::vector<NT>            C_auxiliary;
    
    
    typedef  CGAL::Access_by_index< C_iterator,
                 false, false >         Access_c_B;
    typedef  CGAL::Join_random_access_iterator_1<
                 Index_iterator, Access_c_B >
                                        c_B_iterator;
    
    
    typedef  CGAL::Access_by_index<
                 typename std::iterator_traits<D_iterator>::value_type,
                 false, true >          Access_D_Bj;
    typedef  CGAL::Join_random_access_iterator_1<
                 Index_iterator, Access_D_Bj >
                                        D_Bj_iterator;
    
    

  public:
    
    // public types
    enum Status { UPDATE, INFEASIBLE, UNBOUNDED, OPTIMAL };
    
    
    typedef  Index_iterator             Basic_variable_index_iterator;
    typedef  Value_iterator             Basic_variable_numerator_iterator;
    typedef  CGAL::Join_random_access_iterator_1<
                 Basic_variable_numerator_iterator, Make_quotient >
                                        Basic_variable_value_iterator;
    
    
    typedef  CGAL::Join_random_access_iterator_1<
                 Index_iterator, Access_value_checked >
                                        Variable_numerator_iterator;
    typedef  CGAL::Join_random_access_iterator_1<
                 Variable_numerator_iterator, Make_quotient >
                                        Variable_value_iterator;
    
    
    typedef  Value_iterator             Lambda_numerator_iterator;
    typedef  CGAL::Join_random_access_iterator_1<
                 Lambda_numerator_iterator, Make_quotient >
                                        Lambda_value_iterator;
    
    
    typedef  CGAL::Pricing_strategy_base<Rep>
                                        Pricing_strategy;
    typedef  CGAL::Full_exact_pricing<Rep>
                                        Pricing_strategy_default;
    
    
    typedef  typename A_artificial::const_iterator
                                        A_artificial_iterator;
    
    
    typedef  typename C_auxiliary::const_iterator
                                        C_auxiliary_iterator;
    
    

    
    // creation
    QP_solver( int verbose = 0, std::ostream& stream = std::cout);
      /*
        : nt_0( 0), nt_1( 1), nt_minus_1( -nt_1),
          et_0( 0), et_1( 1), et_2( 2),
          vout1( verbose == 1, stream),
          vout2( verbose >= 2, stream),
          vout3( verbose == 3, stream),
          vout ( verbose >  0, stream),
          qp_n( 0), qp_m( 0), inv_M_B( vout3), m_phase( 0),
          d( inv_M_B.denominator()),
          
          strategyP( &strategy_default)
                                       
        {
            
            CGAL_optimisation_debug {
                vout2 << "======================================" << std::endl
                      << "The CGAL Solver for Quadratic Programs" << std::endl
                      << "======================================" << std::endl;
            }
             
        }
      */
    
    // set-up of QP
    void  set( int n, int m, int max_b,
               A_iterator A, B_iterator b, C_iterator c, D_iterator D);
    
    
    // initialization (of phase I)
    void  init( );
    
    
    // access to QP
    int  number_of_variables  ( ) const { return qp_n; }
    int  number_of_constraints( ) const { return qp_m; }
    
    A_iterator  a_begin( ) const { return qp_A;      }
    A_iterator  a_end  ( ) const { return qp_A+qp_n; }
    
    B_iterator  b_begin( ) const { return qp_b;      }
    B_iterator  b_end  ( ) const { return qp_b+qp_m; }
    
    C_iterator  c_begin( ) const { return qp_c;      }
    C_iterator  c_end  ( ) const { return qp_c+qp_n; }
    
    D_iterator  d_begin( ) const { return qp_D;      }
    D_iterator  d_end  ( ) const { return qp_D+qp_n; }
    
    
    // access to current status
    int     phase     ( ) const { return m_phase; }
    Status  status    ( ) const { return m_status; }
    int     iterations( ) const { return m_iterations; }
    bool    is_optimal( ) const { return status() == OPTIMAL; }
    
    
    // access to common denominator
    const ET&  variables_common_denominator( ) const { return d; }
    
    // access to basic variables
    int  number_of_basic_variables( ) const { return B.size(); }
    
    Basic_variable_index_iterator
    basic_variables_index_begin( ) const { return B.begin(); }
    
    Basic_variable_index_iterator
    basic_variables_index_end  ( ) const { return B.end(); }
    
    Basic_variable_numerator_iterator
    basic_variables_numerator_begin( ) const { return x_B.begin(); }
    
    Basic_variable_numerator_iterator
    basic_variables_numerator_end  ( ) const { return x_B.end(); }
    
    Basic_variable_value_iterator
    basic_variables_value_begin( ) const
        { return Basic_variable_value_iterator(
                     basic_variables_numerator_begin(),
                     Make_quotient( Compute_quotient(), d)); }
    
    Basic_variable_value_iterator
    basic_variables_value_end  ( ) const
        { return Basic_variable_value_iterator(
                     basic_variables_numerator_end  (),
                     Make_quotient( Compute_quotient(), d)); }
    
    
    bool  is_basic( int i) const
        { CGAL_optimisation_precondition( i >= 0);
          CGAL_optimisation_precondition( i < ( phase() == 1
              ? number_of_variables()+number_of_constraints()
              : number_of_variables()));
          return ( in_B[ i] >= 0); }
    
    
    // access to variables
    Variable_numerator_iterator
    variables_numerator_begin( ) const
        { return Variable_numerator_iterator(
                     in_B.begin(), Access_value_checked( x_B.begin(), et_0)); }
    
    Variable_numerator_iterator
    variables_numerator_end  ( ) const
        { return Variable_numerator_iterator(
                     in_B.end  (), Access_value_checked( x_B.begin(), et_0)); }
    
    Variable_value_iterator
    variables_value_begin( ) const
        { return Variable_value_iterator(
                     variables_numerator_begin(),
                     Make_quotient( Compute_quotient(), d)); }
    
    Variable_value_iterator
    variables_value_end  ( ) const
        { return Variable_value_iterator(
                     variables_numerator_end  (),
                     Make_quotient( Compute_quotient(), d)); }
    
    
    // access to current solution
    ET  solution_numerator  ( ) const;
    ET  solution_denominator( ) const { return d*d; }
    
    Quotient<ET>  solution( ) const
        { return Quotient<ET>( solution_numerator(), d*d); }
    
    
    // access to lambda
    Lambda_numerator_iterator
    lambda_numerator_begin( ) const { return lambda.begin(); }
    
    Lambda_numerator_iterator
    lambda_numerator_end  ( ) const { return lambda.end(); }
    
    Lambda_value_iterator
    lambda_value_begin( ) const
        { return Lambda_value_iterator(
                     lambda_numerator_begin(),
                     Make_quotient( Compute_quotient(), d)); }
    
    Lambda_value_iterator
    lambda_value_end  ( ) const
        { return Lambda_value_iterator(
                     lambda_numerator_end  (),
                     Make_quotient( Compute_quotient(), d)); }
    
    // access to auxiliary problem
    A_artificial_iterator  a_artificial_begin( ) const { return art_A.begin();}
    A_artificial_iterator  a_artificial_end  ( ) const { return art_A.end  ();}
    
    C_auxiliary_iterator   c_auxiliary_begin ( ) const { return aux_c.begin();}
    C_auxiliary_iterator   c_auxiliary_end   ( ) const { return aux_c.end  ();}
    
    
    // access to variables of dual LP
    ET  dual_variable( int i) const;
    
    
    // operations
    Status  pivot( )
        { CGAL_optimisation_precondition( phase() > 0);
          CGAL_optimisation_precondition( phase() < 3);
          pivot_step();
          return status(); }
    
    Status  solve( )
        { CGAL_optimisation_precondition( phase() > 0);
          while ( phase() < 3) { pivot_step(); }
          return status(); }
    
    
    // altering the pricing strategy
    const Pricing_strategy&  pricing_strategy() const;
    /*
      { return *strategyP; }
    */
    
    void  set_pricing_strategy( Pricing_strategy& pricing_strategy);
    /*
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl
                      << "-----------------------" << std::endl
                      << "Pricing Strategy Change" << std::endl
                      << "-----------------------" << std::endl;
            }
             
    
            strategyP = &pricing_strategy;
            strategyP->set( *this, vout2);
        }
    */
    
    // transition (to phase II)
    void  transition( );
    
    

  private:
    
    // some constants
    const NT                 nt_0, nt_1, nt_minus_1;
    const ET                 et_0, et_1, et_2;
                                              

    
    // verbose output streams
    CGAL::Verbose_ostream    vout1;     // used for some verbose output
    CGAL::Verbose_ostream    vout2;     // used for more verbose output
    CGAL::Verbose_ostream    vout3;     // used for very verbose output
    CGAL::Verbose_ostream    vout;      // used for any  verbose output
                                                                       
    
    
    // given QP
    int                      qp_n;      // number of variables
    int                      qp_m;      // number of constraints
    
    A_iterator               qp_A;      // constraint matrix
    B_iterator               qp_b;      // right-hand-side vector
    C_iterator               qp_c;      // objective vector
    D_iterator               qp_D;      // objective matrix
    
    // HACK
    unsigned int             max_basis; 
    
    // current status
    Indices                  B;         // basis
    
    Basis_inverse<ET,Is_lp>  inv_M_B;   // inverse of basis matrix
    
    Values                   lambda;
    Values                   x_B;       // solution restricted to basis
    
    int                      m_phase;   // phase of the Simplex method
    Status                   m_status;  // status of last pivot step
    int                      m_iterations;// number of pivot steps
    
    
    const ET&                d;         // reference to `inv_M_B.denominator()'
    
    
    Indices                  in_B;      // position in basis, -1 if non-basic
    
    
    // pricing strategy
    //Pricing_strategy_default strategy_default;
    Pricing_strategy*        strategyP;
    
    
    int                      art_basic; // number of basic artificial variables
    
    
    // auxiliary problem
    A_artificial             art_A;     // artificial part of constraint matrix
    
    C_auxiliary              aux_c;     // auxiliary objective vector
    
    
    // additional variables
    Values                   b;         // exact version of `qp_b'
    Values                   minus_c_B; // exact version of `-qp_c'
                                        // restricted to basis
    
    
    bool                     is_phase_I;// flag indicating phase I
    
    
    int                      j;         // index of entering variable `x_j'
    
    
    Values                   A_j;       // exact version of j-th column of A
    Values                   two_D_Bj;  // exact version of twice the j-th
                                        // column of D restricted to basis B
    
    
    Values                   q_lambda;  // used in the ratio test
    Values                   q_x;       // --------- " ----------
    
    
    int                      i;         // index of leaving variable `x_i'
    ET                       x_i;       // numerator of leaving variable `x_i'
    ET                       q_i;       // corresponding `q_i'
    
    
    ET                       mu_j;      //   numerator of `t_j'
    ET                       nu;        // denominator of `t_j'
    
    

    
    // ------------------------------------------------------------------------
    
    // pivot step
    void  pivot_step( );
    
    
    // set-up functions
    void  set_up_auxiliary_problem( );
    
    void  set_up_basis( );
    
    void  set_up_initial_solution( );
    
    void  set_up_additional_variables( );
    
    
    // pricing
    // -------
    void  pricing( )
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl
                      << "Pricing" << std::endl
                      << "-------" << std::endl;
            }
             
    
            if ( is_phase_I && ( art_basic == 0)) {
                j = -1;                 // all artificial variables non-basic
    
                
                CGAL_optimisation_debug {
                    vout2 << "artificial variables are non-basic" << std::endl;
                }
                 
            } else {
    
                // call pricing strategy
                j = strategyP->pricing();
    
                
                CGAL_optimisation_debug {
                    if ( j < 0) {
                        vout2 << "entering variable: none" << std::endl;
                    } else {
                        vout1 << "  ";
                        vout << "entering"; vout2 << " variable"; vout << ": ";
                        vout  << j; vout2 << std::endl;
                    }
                }
                 
            }
        }
    
    /*
    // ratio test
    // ----------
    void  ratio_test( )
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl
                      << "Ratio Test" << std::endl
                      << "----------" << std::endl;
            }
             
    
            // compute `q_lambda' and `q_x'
            
            CGAL_optimisation_debug {
                
                vout2 << "     A_j: ";
                if ( vout2.verbose()) {
                    std::copy( A_j.begin(), A_j.begin()+qp_m,
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
            
            compute_q( Is_lp());
            
            
            CGAL_optimisation_debug {
                
                vout2 << "     q_x: ";
                if ( vout2.verbose()) {
                    std::copy( q_x.begin(), q_x.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
    
            // check `t_i's
            
            x_i = et_1;                 // trick: initialize
            q_i = et_0;                 // minimum with +oo
            
            
            Value_iterator  x_it = x_B.begin();
            Value_iterator  q_it = q_x.begin();
            for ( unsigned int k = 0; k < B.size(); ++k, ++x_it, ++q_it) {
                if ( ( *q_it > et_0) && ( ( *x_it * q_i) < ( x_i * *q_it))) {
                    i = k; x_i = *x_it; q_i = *q_it;
                }
            }
            
            
    
            // check `t_j'
            
            check_t_j( Is_lp());
                                
    
            
            CGAL_optimisation_debug {
                vout2 << std::endl;
                for ( unsigned int k = 0; k < B.size(); ++k) {
                    vout2 << "t_" << k << ": " << x_B[ k] << '/' << q_x[ k]
                          << ( ( q_i > et_0) && ( i == (int)k) ? " *" : "")
                          << std::endl;
                }
                if ( ! ( CGAL::check_tag( Is_lp()) || is_phase_I)) {
                    vout2 << "t_j: " << mu_j << '/' << nu
                          << ( ( q_i > et_0) && ( i < 0) ? " *" : "")
                          << std::endl;
                }
                vout2 << std::endl;
                if ( q_i > et_0) {
                    if ( i < 0) {
                        vout2 << "leaving variable: none" << std::endl;
                    } else {
                        vout1 << ", ";
                        vout  << "leaving"; vout2 << " variable"; vout << ": ";
                        vout  << B[ i];
                        vout2 << " (= B[ " << i << "])" << std::endl;
                    }
                }
            }
             
        }

    
    // initialization of ratio-test/update loop
    void init_ratio_test_update_loop( )
        {
            // store exact version of `A_j' (implicit conversion to ET)
            if ( j < qp_n) {
                // original variable
                std::copy(  qp_A[ j     ],  qp_A[ j     ]+qp_m, A_j.begin());
            } else {
                // artificial variable
                std::copy( art_A[ j-qp_n], art_A[ j-qp_n]+qp_m, A_j.begin());
            }
    
            // store exact version of `2 D_{B,j}'
            store_2_D_Bj( Is_lp());
        }

    // storing of exact version of `2 D_{B,j}'
    void  store_2_D_Bj( Tag_false)      // QP
        {
            if ( j < qp_n) {
                // original variable
                Access_D_Bj  access_D_Bj( qp_D[ j], nt_0, 0, qp_n);
                std::transform( D_Bj_iterator( B.begin(), access_D_Bj),
                                D_Bj_iterator( B.end  (), access_D_Bj),
                                two_D_Bj.begin(),
                                std::bind1st( std::multiplies<ET>(), et_2));
            } else {
                // artificial variable
                std::fill_n( two_D_Bj.begin(), B.size(), et_0);
            }
        }
    
    void  store_2_D_Bj( Tag_true)       // LP
        {
            // nop
        }
    
    
    // computation of `q_lambda' and `q_x'
    void  compute_q( Tag_false)         // QP
        {
            
            CGAL_optimisation_debug {
                
                vout2 << "  2 D_Bj: ";
                if ( vout2.verbose()) {
                    std::copy( two_D_Bj.begin(), two_D_Bj.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
    
            inv_M_B.multiply( A_j.begin(), two_D_Bj.begin(),
                              q_lambda.begin(), q_x.begin());
    
            
            CGAL_optimisation_debug {
                vout2 << std::endl;
                
                vout2 << "q_lambda: ";
                if ( vout2.verbose()) {
                    std::copy( q_lambda.begin(), q_lambda.begin()+qp_m,
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
        }
    
    void  compute_q( Tag_true)          // LP
        {
            inv_M_B.multiply_x( A_j.begin(), q_x.begin());
        }
    
    
    // computation and checking of `t_j'
    void  check_t_j( Tag_false)         // QP
        {
            // compute `nu'
            nu = std::inner_product( q_x.begin(), q_x.end(),
                                     two_D_Bj.begin(),
                 std::inner_product( q_lambda.begin(), q_lambda.end(),
                                     A_j.begin(),
                                     ( j < qp_n) ? -et_2*d*ET( qp_D[ j][ j])
                                                 : et_0));
    
            if ( ! is_phase_I) {
    
                // compute `mu_j'
                mu_j = std::inner_product( x_B.begin(), x_B.end(),
                                           two_D_Bj.begin(),
                       std::inner_product( lambda.begin(), lambda.end(),
                                           A_j.begin(),
                                           d * ( is_phase_I ? ET( aux_c[ j])
                                                            : ET( qp_c[ j]))));
    
                // check `t_j'
                if ( ( nu < et_0) && ( ( mu_j * q_i) > ( x_i * nu))) {
                    i = -1; q_i = et_1;
                }
            }
        }
    
    void  check_t_j( Tag_true)          // LP
        {
            // nop
        }

    
    // update
    // ------
    void  update( )
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl
                      << "Update" << std::endl
                      << "------";
            }
             
    
            // update basis and basis inverse
            
            update_basis( Is_lp());
                                   
    
            // compute current solution
            
            if ( is_phase_I) {
                inv_M_B.multiply_l( minus_c_B.begin(), lambda.begin());
                inv_M_B.multiply_x(         b.begin(),    x_B.begin());
            } else {
                inv_M_B.multiply( b.begin(), minus_c_B.begin(),
                                  lambda.begin(), x_B.begin());
            }
            
            
            CGAL_optimisation_debug {
                vout2 << std::endl;
                
                vout2 << "  -c_B: ";
                if ( vout2.verbose()) {
                    std::copy( minus_c_B.begin(), minus_c_B.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
                
                vout2 << "lambda: ";
                if ( vout2.verbose()) {
                    std::copy( lambda.begin(), lambda.begin()+qp_m,
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
                
                vout2 << "   x_B: ";
                if ( vout2.verbose()) {
                    std::copy( x_B.begin(), x_B.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
        }

    
    // update of basis and basis inverse
    void  update_basis( Tag_false)      // QP
    {
        // append variable to basis
        if ( ( i < 0) || ( B.size() == (unsigned int)qp_m)) {
    
            
            CGAL_optimisation_debug {
                vout2 << std::endl << "--> non-basic variable "
                      << j << " enters basis" << std::endl;
            }
             
    
            // update basis
            unsigned int  k = B.size();
            B.push_back( j);            // `j' enters basis
            in_B[ j] = k;
            if ( is_phase_I && ( j >= qp_n)) ++art_basic;
    
            
            CGAL_optimisation_debug {
                vout2 << "new ";
                
                vout2 << "basis: ";
                if ( vout2.verbose()) {
                    std::copy( B.begin(), B.end(),
                               std::ostream_iterator<int>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
                vout3 << "new basis-inverse:" << std::endl;
            }
             
    
            // update basis inverse
            inv_M_B.append( q_lambda.begin(), q_x.begin(), nu);
    
            // update status
            x_B.push_back( et_0);
            q_x.push_back( et_0);
            if ( k >= two_D_Bj.size()) {
                 two_D_Bj.push_back( et_0);
                minus_c_B.push_back( -( is_phase_I ? ET( aux_c[ j])
					           : ET(  qp_c[ j])));
            } else {
                minus_c_B[ k] = is_phase_I ? -ET( aux_c[ j]) : -ET( qp_c[ j]);
            }
    
            // mark variable as entered
            j = -1;
        }
    
        // remove variable from basis
        if ( i >= 0) {
    
            
            CGAL_optimisation_debug {
                vout2 << std::endl << "<-- basic variable "
                      << B[ i] << " leaves basis" << std::endl;
            }
             
    
            // make `i' the last index in `B'
            int  k = B.size()-1;
            int  l = B[ i];
            if ( i < k) {
                inv_M_B.swap( i, k);
                std::swap( in_B[ B[ i]], in_B[ B[ k]]);
                        B[ i] =         B[ k];
                minus_c_B[ i] = minus_c_B[ k];
                 two_D_Bj[ i] =  two_D_Bj[ k];
            }
    
            // update basis
            B.pop_back();               // `i' leaves basis
            in_B[ l] = -1;
            if ( is_phase_I && ( l >= qp_n)) --art_basic;
    
            
            CGAL_optimisation_debug {
                vout2 << "new ";
                
                vout2 << "basis: ";
                if ( vout2.verbose()) {
                    std::copy( B.begin(), B.end(),
                               std::ostream_iterator<int>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
                vout3 << "new basis-inverse:" << std::endl;
            }
             
    
            // update basis inverse
            inv_M_B.remove( k);
    
            // update status
            x_B.pop_back();
            q_x.pop_back();
    
            // notify pricing strategy
            strategyP->leaving_basis( l);
        }
    }

    
    void  update_basis( Tag_true)       // LP
    {
        
        CGAL_optimisation_debug {
            vout2 << std::endl << "<--> non-basic variable " << j
                  << " replaces basic variable " << B[ i] << std::endl;
        }
         
    
        // update basis
        int  l = B[ i];
        in_B[ B[ i] ] = -1;             // `i' leaves basis
        in_B[    j  ] = i;              // `j' enters basis
        B[ i] = j;                      // `j' replaces `i' in basis
        if ( is_phase_I && ( j >= qp_n)) ++art_basic;
        if ( is_phase_I && ( l >= qp_n)) --art_basic;
    
        
        CGAL_optimisation_debug {
            vout2 << "new ";
            
            vout2 << "basis: ";
            if ( vout2.verbose()) {
                std::copy( B.begin(), B.end(),
                           std::ostream_iterator<int>( vout2.out(), " "));
                vout2.out() << std::endl;
            }
             
            vout3 << "new basis-inverse:" << std::endl;
        }
         
    
        // update basis inverse
        inv_M_B.replace( i, q_x.begin());
    
        // update status
        minus_c_B[ i] = -ET( is_phase_I ? aux_c[ j] : qp_c[ j]);
    
        // mark variable as entered
        j = -1;
    
        // notify pricing strategy
        strategyP->leaving_basis( l);
    }
    */    

    // ratio test
    // ----------
    void  ratio_test( )
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl
                      << "Ratio Test" << std::endl
                      << "----------" << std::endl;
            }
             
    
            // initialize
            
            init_ratio_test();
                              
    
            // compute `q_lambda' and `q_x'
            
            compute_q( Is_lp());
            
            
            CGAL_optimisation_debug {
                
                vout2 << "     q_x: ";
                if ( vout2.verbose()) {
                    std::copy( q_x.begin(), q_x.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
    
            // check `t_i's
            
            x_i = et_1;                 // trick: initialize
            q_i = et_0;                 // minimum with +oo
            
            
            Value_iterator  x_it = x_B.begin();
            Value_iterator  q_it = q_x.begin();
            for ( unsigned int k = 0; k < B.size(); ++k, ++x_it, ++q_it) {
                if ( ( *q_it > et_0) && ( ( *x_it * q_i) < ( x_i * *q_it))) {
                    i = k; x_i = *x_it; q_i = *q_it;
                }
            }
            
            
    
            // check `t_j'
            
            check_t_j( Is_lp());
                                
    
            
            CGAL_optimisation_debug {
                vout2 << std::endl;
                for ( unsigned int k = 0; k < B.size(); ++k) {
                    vout2 << "t_" << k << ": " << x_B[ k] << '/' << q_x[ k]
                          << ( ( q_i > et_0) && ( i == (int)k) ? " *" : "")
                          << std::endl;
                }
                if ( ! ( CGAL::check_tag( Is_lp()) || is_phase_I)) {
                    vout2 << "t_j: " << mu_j << '/' << nu
                          << ( ( q_i > et_0) && ( i < 0) ? " *" : "")
                          << std::endl;
                }
                vout2 << std::endl;
                if ( q_i > et_0) {
                    if ( i < 0) {
                        vout2 << "leaving variable: none" << std::endl;
                    } else {
                        vout1 << ", ";
                        vout  << "leaving"; vout2 << " variable"; vout << ": ";
                        vout  << B[ i];
                        vout2 << " (= B[ " << i << "])" << std::endl;
                    }
                }
            }
             
        }
    
    
    // initialization of ratio-test
    void init_ratio_test( )
        {
            // store exact version of `A_j' (implicit conversion to ET)
            if ( j < qp_n) {
                // original variable
                std::copy(  qp_A[ j     ],  qp_A[ j     ]+qp_m, A_j.begin());
            } else {
                // artificial variable
                std::copy( art_A[ j-qp_n], art_A[ j-qp_n]+qp_m, A_j.begin());
            }
    
            
            CGAL_optimisation_debug {
                
                vout2 << "     A_j: ";
                if ( vout2.verbose()) {
                    std::copy( A_j.begin(), A_j.begin()+qp_m,
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
    
            // store exact version of `2 D_{B,j}'
            store_2_D_Bj( Is_lp());
        }
    
    // storing of exact version of `2 D_{B,j}'
    void  store_2_D_Bj( Tag_false)      // QP
        {
            if ( j < qp_n) {
                // original variable
                Access_D_Bj  access_D_Bj( qp_D[ j], nt_0, 0, qp_n);
                std::transform( D_Bj_iterator( B.begin(), access_D_Bj),
                                D_Bj_iterator( B.end  (), access_D_Bj),
                                two_D_Bj.begin(),
                                std::bind1st( std::multiplies<ET>(), et_2));
            } else {
                // artificial variable
                std::fill_n( two_D_Bj.begin(), B.size(), et_0);
            }
    
            
            CGAL_optimisation_debug {
                
                vout2 << "  2 D_Bj: ";
                if ( vout2.verbose()) {
                    std::copy( two_D_Bj.begin(), two_D_Bj.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
        }
    
    void  store_2_D_Bj( Tag_true)       // LP
        {
            // nop
        }
    
    
    // computation of `q_lambda' and `q_x'
    void  compute_q( Tag_false)         // QP
        {
            inv_M_B.multiply( A_j.begin(), two_D_Bj.begin(),
                              q_lambda.begin(), q_x.begin());
    
            
            CGAL_optimisation_debug {
                vout2 << std::endl;
                
                vout2 << "q_lambda: ";
                if ( vout2.verbose()) {
                    std::copy( q_lambda.begin(), q_lambda.begin()+qp_m,
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
        }
    
    void  compute_q( Tag_true)          // LP
        {
            inv_M_B.multiply_x( A_j.begin(), q_x.begin());
        }
    
    
    // computation and checking of `t_j'
    void  check_t_j( Tag_false)         // QP
        {
            // compute `nu'
            nu = std::inner_product( q_x.begin(), q_x.end(),
                                     two_D_Bj.begin(),
                 std::inner_product( q_lambda.begin(), q_lambda.end(),
                                     A_j.begin(),
                                     ( j < qp_n) ? -et_2*d*ET( qp_D[ j][ j])
                                                 : et_0));
    
            if ( ! is_phase_I) {
    
                // compute `mu_j'
                mu_j = std::inner_product( x_B.begin(), x_B.end(),
                                           two_D_Bj.begin(),
                       std::inner_product( lambda.begin(), lambda.end(),
                                           A_j.begin(),
                                           d * ( is_phase_I ? ET( aux_c[ j])
                                                            : ET( qp_c[ j]))));
    
                // check `t_j'
                if ( ( nu < et_0) && ( ( mu_j * q_i) > ( x_i * nu))) {
                    i = -1; q_i = et_1;
                }
            }
        }
    
    void  check_t_j( Tag_true)          // LP
        {
            // nop
        }
    
    
    // update
    // ------
    void  update( )
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl
                      << "Update" << std::endl
                      << "------";
            }
             
    
            // update basis and basis inverse
            
            update_basis( Is_lp());
                                   
    
            // compute current solution
            
            compute_current_solution();
                                       
        }
    
    
    // append variable to basis
    void  append_variable( )
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl << "--> non-basic variable "
                      << j << " enters basis" << std::endl;
            }
             
    
            // update basis
            unsigned int  k = B.size();
            B.push_back( j);            // `j' enters basis
            in_B[ j] = k;
            if ( is_phase_I && ( j >= qp_n)) ++art_basic;
    
            
            CGAL_optimisation_debug {
                vout2 << "new ";
                
                vout2 << "basis: ";
                if ( vout2.verbose()) {
                    std::copy( B.begin(), B.end(),
                               std::ostream_iterator<int>(vout2.out()," "));
                    vout2.out() << std::endl;
                }
                 
                vout3 << "new basis-inverse:" << std::endl;
            }
             
    
            // update basis inverse
            inv_M_B.append( q_lambda.begin(), q_x.begin(), nu);
    
            // update status
            x_B.push_back( et_0);
            q_x.push_back( et_0);
            if ( k >= two_D_Bj.size()) {
                two_D_Bj.push_back( et_0);
                minus_c_B.push_back( -( is_phase_I ? ET( aux_c[ j])
					           : ET(  qp_c[ j])));
            } else {
                minus_c_B[ k] = is_phase_I ? -ET( aux_c[ j]) : -ET( qp_c[ j]);
            }

	    CGAL_optimisation_assertion( check_basis( Is_lp()));
        }
    
    // remove variable from basis
    void  remove_variable( )
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl << "<-- basic variable "
                      << B[ i] << " leaves basis" << std::endl;
            }
             
    
            // make `i' the last index in `B'
            int  k = B.size()-1;
            int  l = B[ i];
            if ( i < k) {
                inv_M_B.swap( i, k);
                std::swap( in_B[ B[ i]], in_B[ B[ k]]);
                        B[ i] =         B[ k];
                minus_c_B[ i] = minus_c_B[ k];
                 two_D_Bj[ i] =  two_D_Bj[ k];
            }
    
            // update basis
            B.pop_back();               // `i' leaves basis
            in_B[ l] = -1;
            if ( is_phase_I && ( l >= qp_n)) --art_basic;
    
            
            CGAL_optimisation_debug {
                vout2 << "new ";
                
                vout2 << "basis: ";
                if ( vout2.verbose()) {
                    std::copy( B.begin(), B.end(),
                               std::ostream_iterator<int>(vout2.out()," "));
                    vout2.out() << std::endl;
                }
                 
                vout3 << "new basis-inverse:" << std::endl;
            }
             
    
            // update basis inverse
            inv_M_B.remove( k);
    
            // update status
            x_B.pop_back();
            q_x.pop_back();
    
            // notify pricing strategy
            strategyP->leaving_basis( l);
    
            // new basis found?
            if ( B.size() == (unsigned int)qp_m) { i = -1; }

	    CGAL_optimisation_assertion( check_basis( Is_lp()));
        }
    
    // update of basis and basis inverse
    void  update_basis( Tag_false)      // QP
        {
	  if ( B.size() < max_basis) {

            // append variable to basis
            append_variable();
    
            // remove variable from basis (if any)
            if ( i >= 0) remove_variable();

	  } else {
	    remove_variable();
            inv_M_B.multiply( A_j.begin(), two_D_Bj.begin(),
                              q_lambda.begin(), q_x.begin());
            nu = std::inner_product( q_x.begin(), q_x.end(),
                                     two_D_Bj.begin(),
                 std::inner_product( q_lambda.begin(), q_lambda.end(),
                                     A_j.begin(),
                                     ( j < qp_n) ? -et_2*d*ET( qp_D[ j][ j])
                                                 : et_0));
            append_variable();
	  }
        }
    
    
    void  update_basis( Tag_true)       // LP
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl << "<--> non-basic variable " << j
                      << " replaces basic variable " << B[ i] << std::endl;
            }
             
    
            // update basis
            int  l = B[ i];
            in_B[ l ] = -1;             // `i' leaves basis
            in_B[ j ] = i;              // `j' enters basis
            B[ i] = j;                  // `j' replaces `i' in basis
            if ( is_phase_I) {
                if ( j >= qp_n) ++art_basic;
                if ( l >= qp_n) --art_basic;
            }
            
            CGAL_optimisation_debug {
                vout2 << "new ";
                
                vout2 << "basis: ";
                if ( vout2.verbose()) {
                    std::copy( B.begin(), B.end(),
                               std::ostream_iterator<int>(vout2.out()," "));
                    vout2.out() << std::endl;
                }
                 
                vout3 << "new basis-inverse:" << std::endl;
            }
             
    
            // update basis inverse
            inv_M_B.replace( i, q_x.begin());
    
            // update status
            minus_c_B[ i] = -ET( is_phase_I ? aux_c[ j] : qp_c[ j]);
    
            // notify pricing strategy
            strategyP->leaving_basis( l);
    
            // new basis found
            i = -1;
        }
    
    
    // computation of current solution
    void   compute_current_solution( )
        {
            if ( is_phase_I) {
                inv_M_B.multiply_l( minus_c_B.begin(), lambda.begin());
                inv_M_B.multiply_x(         b.begin(),    x_B.begin());
            } else {
                inv_M_B.multiply( b.begin(), minus_c_B.begin(),
                                  lambda.begin(), x_B.begin());
            }
    
            
            CGAL_optimisation_debug {
                vout2 << std::endl;
                
                vout2 << "  -c_B: ";
                if ( vout2.verbose()) {
                    std::copy( minus_c_B.begin(), minus_c_B.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
                
                vout2 << "lambda: ";
                if ( vout2.verbose()) {
                    std::copy( lambda.begin(), lambda.begin()+qp_m,
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
                
                vout2 << "   x_B: ";
                if ( vout2.verbose()) {
                    std::copy( x_B.begin(), x_B.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
        }
    
    
    // iterated ratio test and update
    // ------------------------------
    void  iterated_ratio_test_update( )
        {
            loop_until_basis_found( Is_lp());
        }
    
    // loop until new basis is found
    void  loop_until_basis_found( Tag_false)    // QP
        {
            while ( 
                    i >= 0
                          ) {
    
                // ratio test (iterated)
                ratio_test_iterated();
    
                // check for unboundedness
                if ( 
                     q_i == et_0
                                ) {
                    m_phase  = 3;
                    m_status = UNBOUNDED;
                    
                    CGAL_optimisation_debug {
                        vout1 << "  ";
                        vout << "problem is UNBOUNDED" << std::endl;
                    }
                     
                    return;
                }
    
                // update (iterated)
                update_iterated();
            }
        }
    
    void  loop_until_basis_found( Tag_true)     // LP
        {
            // nop
        }
    
    
    // ratio test (iterated)
    void  ratio_test_iterated( )
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl
                      << "Ratio Test (iterated)" << std::endl
                      << "---------------------" << std::endl;
            }
             
    
            // get `q_x'
            unsigned int  l = in_B[ j];
            std::copy( inv_M_B.column_begin( qp_m+l)+qp_m,
                       inv_M_B.column_end  ( qp_m+l)     , q_x.begin());
    
            
            CGAL_optimisation_debug {
                
		vout2 << "q_x[ " << l << "]: ";
                if ( vout2.verbose()) {
                    std::copy( q_x.begin(), q_x.begin()+B.size(),
                               std::ostream_iterator<ET>( vout2.out(), " "));
                    vout2.out() << std::endl;
                }
                 
            }
             
    
            // store frequently used values
            ET            x_j = x_B[ l];
            ET            q_j = q_x[ l];
            int  sign_q_i_q_j = CGAL_NTS sign( q_j);
    
            // initialize minimum
            x_i = ( sign_q_i_q_j < 0) ? et_1 : -et_1;   // trick: initialize
            q_i = et_0;                                 // minimum with +oo
    
            // check `t_i's
            
            Value_iterator  x_it = x_B.begin();
            Value_iterator  q_it = q_x.begin();
            for ( unsigned int k = 0; k < B.size(); ++k, ++x_it, ++q_it) {
                if (   ( k != l)
            
                  // t_i > 0 ?
                    && (/*   (( *q_it > et_0) && ( *x_it * q_j < x_j * *q_it))
			||*/ (( *q_it < et_0) && ( *x_it * q_j > x_j * *q_it)))
            
                  // t_i < t_min ?
                    && (   (   ( CGAL_NTS sign( *q_it) * sign_q_i_q_j > 0)
                            && ( *x_it * q_i > x_i * *q_it))
                        || (   ( CGAL_NTS sign( *q_it) * sign_q_i_q_j < 0)
                            && ( *x_it * q_i < x_i * *q_it)))) {
            
                    // store new minimum
                    i = k; x_i = *x_it; q_i = *q_it;
                    sign_q_i_q_j = CGAL_NTS sign( q_j) * CGAL_NTS sign( q_i);
                }
            }
             
    
            // check `t_j'
            
            if ( ( x_j > et_0) && ( CGAL_NTS sign( x_i) * sign_q_i_q_j <= 0)) {
                i = -1; q_i = et_1;
            }
             
    
            
            CGAL_optimisation_debug {
                vout2 << std::endl;
                for ( unsigned int k = 0; k < B.size(); ++k) {
                    if ( k != l) {
                        vout2 << "t_" << k << ": "
                              << x_j * q_x[k] - x_B[k] * q_j << '/' << q_x[k]
                              << ( ( q_i != et_0) && (i == (int)k) ? " *" : "")
                              << std::endl;
            
                    }
                }
                vout2 << "t_j: " << x_j << '/' << et_1
                      << ( ( q_i > et_0) && ( i < 0) ? " *" : "")
                      << std::endl << std::endl;
                if ( q_i != et_0) {
                    if ( i < 0) {
                        vout2 << "leaving variable: none" << std::endl;
                    } else {
                        vout1 << ", ";
                        vout  << "leaving"; vout2 << " variable"; vout << ": ";
                        vout  << B[ i];
                        vout2 << " (= B[ " << i << "])" << std::endl;
                    }
                }
            }
             
        }
    
    
    // update (iterated)
    void  update_iterated( )
        {
            if ( i >= 0) {
    
                
                CGAL_optimisation_debug {
                    vout2 << std::endl
                          << "Update (iterated)" << std::endl
                          << "-----------------";
                }
                 
    
                // update basis and basis inverse
                remove_variable();
    
		CGAL_optimisation_assertion( check_basis( Is_lp()));

                // compute current solution
                compute_current_solution();
            }
        }

    bool check_basis( Tag_true)       // LP
        {
            return true;
        }

    bool check_basis( Tag_false)
	{
	    if ( is_phase_I) return true;

	    Values  result( qp_m+B.size());
	    Values  col_l( qp_m, et_0), col_x( B.size());
	    unsigned int i, j;

	    // first part
	    for ( i = 0; i < (unsigned int)qp_m; ++i) {

		// get source column
		for ( j = 0; j < B.size(); ++j) col_x[ j] = qp_A[ B[ j]][ i];

		// compute target column
		inv_M_B.multiply( col_l.begin(), col_x.begin(),
				  result.begin(), result.begin()+qp_m);

		// check result
		/*
		std::copy( result.begin(), result.end(),
			   std::ostream_iterator<ET>( std::cerr, " "));
		std::cerr << endl;
		*/
		for ( j = 0; j < qp_m+B.size(); ++j) {
		    if ( ( ( j == i) && result[ j] != d) ||
			 ( ( j != i) && result[ j] != et_0)) return false;
		}
	    }

	    // second part
	    for ( i = 0; i < B.size(); ++i) {
		j = B[ i];

		// get source column
                std::copy( qp_A[ j], qp_A[ j]+qp_m, col_l.begin());
                Access_D_Bj  access_D_Bj( qp_D[ j], nt_0, 0, qp_n);
                std::transform( D_Bj_iterator( B.begin(), access_D_Bj),
                                D_Bj_iterator( B.end  (), access_D_Bj),
                                col_x.begin(),
                                std::bind1st( std::multiplies<ET>(), et_2));

		// compute target column
		inv_M_B.multiply( col_l.begin(), col_x.begin(),
				  result.begin(), result.begin()+qp_m);

		// check result
		/*
		std::copy( result.begin(), result.end(),
			   std::ostream_iterator<ET>( std::cerr, " "));
		std::cerr << endl;
		*/
		for ( j = 0; j < qp_m+B.size(); ++j) {
		    if ( ( ( j == i+qp_m) && result[ j] != d) ||
			 ( ( j != i+qp_m) && result[ j] != et_0)) return false;
		}
	    }
	    return true;
	}
    
};
  
CGAL_END_NAMESPACE

#include <CGAL/_QP_solver/Pricing_strategy_base.h>
#include <CGAL/_QP_solver/Full_exact_pricing.h>

CGAL_BEGIN_NAMESPACE

template < class Rep_ >
QP_solver<Rep_>::
    QP_solver( int verbose, std::ostream& stream)
        : nt_0( 0), nt_1( 1), nt_minus_1( -nt_1),
          et_0( 0), et_1( 1), et_2( 2),
          vout1( verbose == 1, stream),
          vout2( verbose >= 2, stream),
          vout3( verbose == 3, stream),
          vout ( verbose >  0, stream),
          qp_n( 0), qp_m( 0), inv_M_B( vout3), m_phase( 0),
          d( inv_M_B.denominator()),
          
          strategyP( new Pricing_strategy_default)
                                       
        {
            
            CGAL_optimisation_debug {
                vout2 << "======================================" << std::endl
                      << "The CGAL Solver for Quadratic Programs" << std::endl
                      << "======================================" << std::endl;
            }
             
        }

template < class Rep_ >
const QP_solver<Rep_>::Pricing_strategy&
QP_solver<Rep_>::
    pricing_strategy() const { return *strategyP; }
    
template < class Rep_ >
void
QP_solver<Rep_>::
    set_pricing_strategy( QP_solver<Rep_>::Pricing_strategy& pricing_strategy)
        {
            
            CGAL_optimisation_debug {
                vout2 << std::endl
                      << "-----------------------" << std::endl
                      << "Pricing Strategy Change" << std::endl
                      << "-----------------------" << std::endl;
            }
             
    
            strategyP = &pricing_strategy;
            strategyP->set( *this, vout2);
        }

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/_QP_solver/QP_solver.C>
#endif

#endif // CGAL_QP_SOLVER_H

// ===== EOF ==================================================================
