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
// file          : include/CGAL/_QP_solver/Partial_filtered_pricing.h
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.4
// revision_date : 2000/08/17
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Pricing Strategy with partial filtered pricing
// ============================================================================
                                                                               

#ifndef CGAL_PARTIAL_FILTERED_PRICING_H
#define CGAL_PARTIAL_FILTERED_PRICING_H

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
class Partial_filtered_pricing;
                               

// Class interface
// ===============
template < class _Rep >
class Partial_filtered_pricing
    : public CGAL::Pricing_strategy_base<_Rep> {
  public:
    // self
    typedef  _Rep                        Rep;
    typedef  Partial_filtered_pricing<Rep>  Self;
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
      NT  nt_0, nt_1, nt_2;
      ET  et_0,       et_2;

      // data members
      std::vector<int>   N;         // non-basis
      int                s;         // size of active set
      std::vector<NT>    row_max_A;
      std::vector<NT>    row_max_D;
      std::vector<bool>  row_valid;
      NT                 row_max_c;
      std::vector<NT>    col_max;

  public:
    
    // creation
    Partial_filtered_pricing( )
        : nt_0( 0), nt_1( 1), nt_2( 2), et_0( 0), et_2( 2) { }
    
    
    // initialization
    void  set( )
    {
        CGAL_optimisation_debug {
            vout() << "partial filtered pricing" << std::endl;
        }
    }
    
    void  init( )
    {
        int i, j;
    
        const Solver& solve = solver();
        int  n = solve.number_of_variables();
        int  m = solve.number_of_constraints();
        s = min( 2*m, n);
        N.erase( N.begin(), N.end());
        N.reserve( n);
        for ( i = 0; i < n; ++i) N.push_back( i);
    
        // compute maxima
        row_max_A = std::vector<NT>(   m, nt_1);
        col_max   = std::vector<NT>( n+m, nt_0);
        A_iterator  Aj = solve.a_begin();
        NT z;
        for ( j = 0; j < n; ++j, ++Aj) {
            for ( i = 0; i < m; ++i) {
                z = CGAL::NTS::abs( (*Aj)[ i]);
                if ( z > row_max_A[ i]) row_max_A[ i] = z;
                if ( z > col_max  [ j]) col_max  [ j] = z;
            }
        }
        for ( j = n; j < n+m; ++j) {
            col_max[ j] = nt_1;
        }
        row_max_c = nt_1;
    }
    
    void  transition( )
    {
        const Solver& solve = solver();
        int  n = solve.number_of_variables();
        int  m = solve.number_of_constraints();
    
        // remove artificial variables from N
        int j, i = 0;
        for ( j = n-m; j < n; ++j) {
            if ( N[ j] < n) {
                while ( N[ i] < n) { ++i; }
                N[ i] = N[ j];
            }
        }
        N.erase( N.end()-m, N.end());
        s = min( (int)(m*CGAL::NTS::sqrt<double>(n)), n-m);
    
        // update row/column maxima of `A'
        C_iterator  c_i = solve.c_begin();
        NT z;
        for ( i = 0; i < n; ++i, ++c_i) {
            z = CGAL::NTS::abs( *c_i);
            if ( z > col_max[ i]) col_max[ i] = z;
            if ( z > row_max_c  ) row_max_c   = z;
        }
    
        // compute row/column maxima of `D'
        if ( ! CGAL::check_tag( Is_lp())) {
            row_max_D = std::vector< NT >( n, nt_0);
            row_valid = std::vector<bool>( n, false);
        }
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
        int  m = solve.number_of_constraints();
        int  b = solve.number_of_basic_variables();
        ET   d = solve.variables_common_denominator();
        NT   nt_d = CGAL::to_double( d);
    
        int   i, j, k, min_k  = -1, min_j = -1;
        NT    nt_mu, nt_min_mu =  0;
        ET    mu, min_mu =  0;
        bool  is_phase_I = ( solve.phase() == 1);
    
        // get inexact versions of `lambda' and `x_B'
        std::vector<NT>  lambda, x_B;
        lambda.reserve( m);
        std::transform( solve.lambda_numerator_begin(),
                        solve.lambda_numerator_end(),
                        std::back_inserter( lambda), To_double());
        if ( ! ( CGAL::check_tag( Is_lp()) || is_phase_I)) {
            x_B.reserve( b);
            std::transform( solve.basic_variables_numerator_begin(),
                            solve.basic_variables_numerator_end(),
                            std::back_inserter( x_B), To_double());
        }
    
        // loop over all active non-basic variables
        for ( k = 0; k < s; ++k) {
    
            j = N[ k];
    
            // compute mu_j
            if ( is_phase_I) {      // phase I
                if ( j < n) {          // original variable
                    nt_mu = std::inner_product(
                        lambda.begin(), lambda.end(),
                        solve.a_begin()[ j],
                        nt_d * solve.c_auxiliary_begin()[ j]);
                } else {               // artificial variable
                    nt_mu = std::inner_product(
                        lambda.begin(), lambda.end(),
                        solve.a_artificial_begin()[ j-n],
                        nt_d * solve.c_auxiliary_begin()[ j]);
                }
            } else {                // phase II
                nt_mu = std::inner_product(
                    lambda.begin(), lambda.end(),
                    solve.a_begin()[ j],
                    nt_d * solve.c_begin()[ j]);
                // is QP?
                if ( ! CGAL::check_tag( Is_lp())) {
                    nt_mu += nt_2 * std::inner_product(
                        x_B.begin(), x_B.end(),
                        D_Bj_iterator( solve.basic_variables_index_begin(),
                                       Access_D_Bj( solve.d_begin()[ j])),
                        nt_0);
                }
            }
    
            CGAL_optimisation_debug {
                vout() << "nt_mu_" << j << ": " << nt_mu << std::endl;
            }
    
            // new minimum?
            if ( ( nt_mu < nt_min_mu) ||
                 ( ( min_j >= n) && ( j < n) && ( nt_mu == nt_min_mu))) {
                min_k  = k;
                min_j  = j;
                nt_min_mu = nt_mu;
            }
        }
    
        // exact check of entering variable
        if ( min_k >= 0) {
            j = N[ min_k];
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
            if ( mu >= et_0) {
                vout() << "entering variable defeated by exact check\n";
                min_k = -1; min_j = -1; nt_min_mu = nt_0;
            }
        }
    
        if ( min_k < 0) {
    
        // --------------------------------------------------------------------
        vout() << "no entering variable found so far, test remaining variables"
        // --------------------------------------------------------------------
               << std::endl;
    
            // loop over all remaining non-basic variables
            for ( k = s; k < (int)N.size(); ++k) {
    
                j = N[ k];
    
                // compute mu_j
                if ( is_phase_I) {      // phase I
                    if ( j < n) {          // original variable
                        nt_mu = std::inner_product(
                            lambda.begin(), lambda.end(),
                            solve.a_begin()[ j],
                            nt_d * solve.c_auxiliary_begin()[ j]);
                    } else {               // artificial variable
                        nt_mu = std::inner_product(
                            lambda.begin(), lambda.end(),
                            solve.a_artificial_begin()[ j-n],
                            nt_d * solve.c_auxiliary_begin()[ j]);
                    }
                } else {                // phase II
                    nt_mu = std::inner_product(
                        lambda.begin(), lambda.end(),
                        solve.a_begin()[ j],
                        nt_d * solve.c_begin()[ j]);
                    // is QP?
                    if ( ! CGAL::check_tag( Is_lp())) {
                        nt_mu += nt_2 * std::inner_product(
                            x_B.begin(), x_B.end(),
                            D_Bj_iterator( solve.basic_variables_index_begin(),
                                           Access_D_Bj( solve.d_begin()[ j])),
                            nt_0);
                    }
                }
    
                CGAL_optimisation_debug {
                    vout() << "nt_mu_" << j << ": " << nt_mu << std::endl;
                }
    
                // improving variable?
                if ( nt_mu < nt_0) {
                    std::swap( N[ k], N[ s]);
    
                    // new minimum?
                    if ( ( nt_mu < nt_min_mu) ||
                         ( ( min_j >= n) && ( j < n) &&
                           ( nt_mu == nt_min_mu))) {
                        min_k  = s;
                        min_j  = j;
                        nt_min_mu = nt_mu;
                    }
    
                    ++s;
                }
            }
    
            // exact check of entering variable
            if ( min_k >= 0) {
                j = N[ min_k];
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
                if ( mu >= et_0) {
                    vout() << "entering variable defeated by exact check\n";
                    min_k = -1; nt_min_mu = nt_0;
                }
            }
        }
        if ( min_k < 0) {
    
            // ----------------------------------------------------------------
            vout()
            << "no entering variable found so far, revert to exact arithmetic"
            // ----------------------------------------------------------------
            << std::endl;
    
            // compute first error bound
            k = m+b+1;
            NT q = ldexp( 1.015625*k*(k+1), -53);
            NT max_1 = nt_d * row_max_c;
            NT max_2 = nt_d;
            NT z;
            for ( i = 0; i < m; ++i) {
                z = CGAL::NTS::abs( lambda[ i]) * row_max_A[ i];
                if ( z > max_1) max_1 = z;
                z = CGAL::NTS::abs( lambda[ i]);
                if ( z > max_2) max_2 = z;
            }
            if ( ! CGAL::check_tag( Is_lp())) {
                typename std::iterator_traits<D_iterator>::value_type  row_D;
                for ( i = 0; i < b; ++i) {
                    k = solve.basic_variables_index_begin()[ i];
                    row_D = solve.d_begin()[ k];
                    if ( ! row_valid[ k]) {
                        NT  max = nt_0;
                        for ( j = 0; j < n; ++j) {
                            z = CGAL::NTS::abs( row_D[ j]);
                            if ( z >     max    )     max     = z;
                            if ( z > col_max[ j]) col_max[ j] = z;
                        }
                        row_max_D[ k] = max;
                    }
                    z = CGAL::NTS::abs( x_B[ i]) * row_max_D[ k];
                    if ( z > max_1) max_1 = z;
                    z = CGAL::NTS::abs( x_B[ i]);
                    if ( z > max_2) max_2 = z;
                }
            }
            NT bound_1 = max_1 * q, max_u_q = max_2 * q, bound_2;
            CGAL_optimisation_debug {
                vout() << "[ first bound: " << bound_1 << " ]" << std::endl;
            }
    
            // loop again over all non-basic variables to verify optimality
            k = 0;
            while ( k < (int)N.size() && min_k < 0) {
    
                j = N[ k];
    
                // compute mu_j (inexact)
                if ( is_phase_I) {      // phase I
                    if ( j < n) {          // original variable
                        nt_mu = std::inner_product(
                            lambda.begin(), lambda.end(),
                            solve.a_begin()[ j],
                            nt_d * solve.c_auxiliary_begin()[ j]);
                    } else {               // artificial variable
                        nt_mu = std::inner_product(
                            lambda.begin(), lambda.end(),
                            solve.a_artificial_begin()[ j-n],
                            nt_d * solve.c_auxiliary_begin()[ j]);
                    }
                } else {                // phase II
                    nt_mu = std::inner_product(
                        lambda.begin(), lambda.end(),
                        solve.a_begin()[ j],
                        nt_d * solve.c_begin()[ j]);
                    // is QP?
                    if ( ! CGAL::check_tag( Is_lp())) {
                        nt_mu += nt_2 * std::inner_product(
                            x_B.begin(), x_B.end(),
                            D_Bj_iterator( solve.basic_variables_index_begin(),
                                           Access_D_Bj( solve.d_begin()[ j])),
                            nt_0);
                    }
                }
    
                CGAL_optimisation_debug {
                    vout() << "nt_mu_" << j << ": " << nt_mu;
                }
    
                // check against first bound
                if ( nt_mu >= bound_1) {
                    CGAL_optimisation_debug {
                        vout() << " [ certified by first bound ]" << std::endl;
                    }
                } else {
                    // compute second bound
                    bound_2 = col_max[ j] * max_u_q;
                    if ( nt_mu >= bound_2) {
                        CGAL_optimisation_debug {
                            vout() << " [ certified by second bound: "
                                   << bound_2 << " ]" << std::endl;
                        }
                    } else {
    
                        // compute mu_j (exact)
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
                            vout() << " [ exact computation needed: "
                                   << mu << " ]" << std::endl;
                        }
    
                        if ( mu < et_0) min_k = k;
                    }
                }
                ++k;
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
    
    
};
  

CGAL_END_NAMESPACE
                  

#endif // CGAL_PARTIAL_FILTERED_PRICING_H

// ===== EOF ==================================================================
