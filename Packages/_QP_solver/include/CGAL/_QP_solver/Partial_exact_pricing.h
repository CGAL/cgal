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
// file          : include/CGAL/_QP_solver/Partial_exact_pricing.h
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.4
// revision_date : 2000/08/17
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Pricing Strategy with partial exact pricing
// ============================================================================
                                                                               

#ifndef CGAL_PARTIAL_EXACT_PRICING_H
#define CGAL_PARTIAL_EXACT_PRICING_H

// includes
#include <CGAL/_QP_solver/Pricing_strategy_base.h>
#include <CGAL/_QP_solver/Join_random_access_iterator.h>
#include <CGAL/_QP_solver/Access_by_index.h>
#include <vector>
#include <numeric>


CGAL_BEGIN_NAMESPACE
                    

// Class declaration
// =================
template < class Rep >
class Partial_exact_pricing;
                            

// Class interface
// ===============
template < class _Rep >
class Partial_exact_pricing
    : public CGAL::Pricing_strategy_base<_Rep> {
  public:
    // self
    typedef  _Rep                        Rep;
    typedef  Partial_exact_pricing<Rep>  Self;
    typedef  Pricing_strategy_base<Rep>  Base;

    // types from the base class
    typedef  typename Base::NT          NT;
    typedef  typename Base::ET          ET;

    typedef  typename Base::A_iterator  A_iterator;
    typedef  typename Base::B_iterator  B_iterator;
    typedef  typename Base::C_iterator  C_iterator;
    typedef  typename Base::D_iterator  D_iterator;

    typedef  typename Base::A_artificial_iterator
                                        A_artificial_iterator;
    typedef  typename Base::C_auxiliary_iterator
                                        C_auxiliary_iterator;

    typedef  typename Base::Basic_variable_index_iterator
                                        Basic_variable_index_iterator;

    typedef  typename Base::Is_lp       Is_lp;

    typedef  typename Base::Solver      Solver;

    typedef  typename Base::Tag_true    Tag_true;
    typedef  typename Base::Tag_false   Tag_false;

  private:
      // some constants
      ET  et_0, et_2;

      // data members
      std::vector<int>   N;         // non-basis
      int                s;         // size of active set

  public:
    
    // creation
    Partial_exact_pricing( ) : et_0( 0), et_2( 2) { }
    
    
    // initialization
    void  set( )
    {
        CGAL_optimisation_debug {
            vout() << "partial exact pricing" << std::endl;
        }
    }
    
    void  init( )
    {
        const Solver& solve = solver();
        int  n = solve.number_of_variables();
        int  m = solve.number_of_constraints();
        s = min( 2*m, n);
        N.erase( N.begin(), N.end());
        N.reserve( n);
        for ( int i = 0; i < n; ++i) N.push_back( i);
    }
    
    // operations
    int  pricing( )
    {
        typedef  CGAL::Access_by_index< CGAL_TYPENAME_MSVC_NULL
                     std::iterator_traits<D_iterator>::value_type,
                     false,false>       Access_D_Bj;
        typedef  CGAL::Join_random_access_iterator_1<
                     Basic_variable_index_iterator,
                     Access_D_Bj >      D_Bj_iterator;
    
        const Solver& solve = solver();
        int  n = solve.number_of_variables();
        ET   d = solve.variables_common_denominator();
    
        int   j,  min_k  = -1, min_j = -1;
        ET    mu, min_mu =  0;
        bool  is_phase_I = ( solve.phase() == 1);
    
        // loop over all active non-basic variables
        for ( int k = 0; k < s; ++k) {
    
            j = N[ k];
    
            // compute mu_j
            if ( is_phase_I) {      // phase I
                if ( j < n) {          // original variable
                    mu = std::inner_product(
                        solve.lambda_numerator_begin(),
                        solve.lambda_numerator_end(),
                        solve.a_begin()[ j],
                        d * solve.c_auxiliary_begin()[ j]);
                } else {               // artificial variable
                    mu = std::inner_product(
                        solve.lambda_numerator_begin(),
                        solve.lambda_numerator_end(),
                        solve.a_artificial_begin()[ j-n],
                        d * solve.c_auxiliary_begin()[ j]);
                }
            } else {                // phase II
                mu = std::inner_product(
                    solve.lambda_numerator_begin(),
                    solve.lambda_numerator_end(),
                    solve.a_begin()[ j],
                    d * solve.c_begin()[ j]);
                // is QP?
                if ( ! CGAL::check_tag( Is_lp())) {
                    mu += et_2 * std::inner_product(
                        solve.basic_variables_numerator_begin(),
                        solve.basic_variables_numerator_end(),
                        D_Bj_iterator( solve.basic_variables_index_begin(),
                                       Access_D_Bj( solve.d_begin()[ j])),
                        et_0);
                }
            }
    
            CGAL_optimisation_debug {
                vout() << "mu_" << j << ": " << mu << std::endl;
            }
    
            // new minimum?
            if ( ( mu < min_mu) ||
                 ( ( min_j >= n) && ( j < n) && ( mu == min_mu))) {
                min_k  = k;
                min_j  = j;
                min_mu = mu;
            }
        }
    
        if ( min_k < 0) {
    
        // --------------------------------------------------------------------
        vout() << "no entering variable found so far, test remaining variables"
        // --------------------------------------------------------------------
               << std::endl;
    
            // loop over all remaining non-basic variables
            for ( int k = s; k < (int)N.size(); ++k) {
    
                j = N[ k];
    
                // compute mu_j
                if ( is_phase_I) {      // phase I
                    if ( j < n) {          // original variable
                        mu = std::inner_product(
                            solve.lambda_numerator_begin(),
                            solve.lambda_numerator_end(),
                            solve.a_begin()[ j],
                            d * solve.c_auxiliary_begin()[ j]);
                    } else {               // artificial variable
                        mu = std::inner_product(
                            solve.lambda_numerator_begin(),
                            solve.lambda_numerator_end(),
                            solve.a_artificial_begin()[ j-n],
                            d * solve.c_auxiliary_begin()[ j]);
                    }
                } else {                // phase II
                    mu = std::inner_product(
                        solve.lambda_numerator_begin(),
                        solve.lambda_numerator_end(),
                        solve.a_begin()[ j],
                        d * solve.c_begin()[ j]);
                    // is QP?
                    if ( ! CGAL::check_tag( Is_lp())) {
                        mu += et_2 * std::inner_product(
                            solve.basic_variables_numerator_begin(),
                            solve.basic_variables_numerator_end(),
                            D_Bj_iterator( solve.basic_variables_index_begin(),
                                           Access_D_Bj( solve.d_begin()[ j])),
                            et_0);
                    }
                }
    
                CGAL_optimisation_debug {
                    vout() << "mu_" << j << ": " << mu << std::endl;
                }
    
                // improving variable?
                if ( mu < et_0) {
                    std::swap( N[ k], N[ s]);
    
                    // new minimum?
                    if ( ( mu < min_mu) ||
                         ( ( min_j >= n) && ( j < n) && ( mu == min_mu))) {
                        min_k  = s;
                        min_j  = j;
                        min_mu = mu;
                    }
    
                    ++s;
                }
            }
        }
        vout() << std::endl;
    
        // return index of entering variable
        if ( min_k >= 0) {
            j = N[ min_k];
            --s;
            N[ min_k] = N[ s];
            N[ s] = N.back();
            N.pop_back();
            return j;
        }
        return -1;
    }
    
    void  leaving_basis( int i)
    {
        if ( s == (int)N.size()) {
            N.push_back( i);
        } else {
            N.push_back( N[ s]);
            N[ s] = i;
        }
        ++s;
    }
    
    void  transition( )
    {
        const Solver& solve = solver();
        int  n = solve.number_of_variables();
        int  m = solve.number_of_constraints();
    
        // remove artificial variables from N
        int i = 0;
        for ( int j = n-m; j < n; ++j) {
            if ( N[ j] < n) {
                while ( N[ i] < n) { ++i; }
                N[ i] = N[ j];
            }
        }
        N.erase( N.end()-m, N.end());
        s = min( (int)(m*CGAL::NTS::sqrt<double>(n)), n-m);
    }
    
    
};
  

CGAL_END_NAMESPACE
                  

#endif // CGAL_PARTIAL_EXACT_PRICING_H

// ===== EOF ==================================================================
