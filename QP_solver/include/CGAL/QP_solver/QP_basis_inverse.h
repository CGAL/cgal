// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp 
//                 Kaspar Fischer 
                                                                               
#ifndef CGAL_QP_SOLVER_QP_BASIS_INVERSE_H
#define CGAL_QP_SOLVER_QP_BASIS_INVERSE_H

#include <CGAL/license/QP_solver.h>


#include <CGAL/QP_solver/basic.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <vector>

namespace CGAL {
                    
// =================
// class declaration
// =================
template < class ET_, class Is_LP_ >
class QP_basis_inverse;

// ===============
// class interface
// ===============
template < class ET_, class Is_LP_ >
class QP_basis_inverse {
  public:
    // self
    typedef  ET_                        ET;
    typedef  Is_LP_                     Is_LP;
    typedef  QP_basis_inverse<ET,Is_LP>
                                        Self;

  private:
    
    // private types
    typedef std::vector<ET>             Row;
    typedef std::vector<Row>            Matrix;

  public:

    /*
     * Note: Some member functions below are suffixed with '_'.
     * They are member templates and their declaration is "hidden",
     * because they are also implemented in the class interface.
     * This is a workaround for M$-VC++, which otherwise fails to
     * instantiate them correctly.
     */

    // creation and initialization
    // ---------------------------
    QP_basis_inverse( CGAL::Verbose_ostream&  verbose_ostream);

    // set-up
    void  set( int n, int m, int nr_equalities);
    
    // init
    template < class InputIterator >                            // phase I
    void  init_( unsigned int art_size, InputIterator art_first);

    /*
    template < class InputIterator >                            // phase II
    void  init_( ...);
    */

    // transition to phase II
    template < class InputIterator >                            // QP case
    void  transition_( InputIterator twice_D_it);

    void  transition( );                                        // LP case

    // access
    // ------
    const ET&  denominator( ) const { return d; }

    const ET&  entry( unsigned int row, unsigned int column) const
        { return entry( row, column, Is_LP()); }

    // multiplication functions
    // ------------------------
    // matrix-vector multiplication ( y = M v )
    template < class ForwardIterator, class OutputIterator >  inline
    void  multiply( ForwardIterator v_l_it, ForwardIterator v_x_it,
                     OutputIterator y_l_it,  OutputIterator y_x_it) const
        { multiply( v_l_it, v_x_it, y_l_it, y_x_it, Is_LP(), Tag_true()); }
    
    // special matrix-vector multiplication functions for LPs
    template < class ForwardIterator, class OutputIterator >  inline
    void  multiply_l( ForwardIterator v_x_it, OutputIterator y_l_it) const
        { CGAL_qpe_assertion( is_LP || is_phaseI);
          multiply__l( v_x_it, y_l_it); }
    
    template < class ForwardIterator, class OutputIterator >  inline
    void  multiply_x( ForwardIterator v_l_it, OutputIterator y_x_it) const
        { CGAL_qpe_assertion( is_LP || is_phaseI);
	  multiply__x( v_l_it, y_x_it); }
    
    // vector-matrix multiplication ( x^T = u^T M )
    template < class ForwardIterator, class OutputIterator >  inline
    void  multiply_transposed( ForwardIterator u_l_it, ForwardIterator u_x_it,
                                OutputIterator x_l_it,  OutputIterator x_x_it)
                                                                          const
        { multiply( u_l_it, u_x_it, x_l_it, x_x_it); } // M_B^{-1} is symmetric
    
    // special vector-matrix multiplication functions for LPs
    template < class ForwardIterator, class OutputIterator >  inline
    void  multiply_transposed_l( ForwardIterator u_x_it,
                                  OutputIterator x_l_it) const
        { multiply_l( u_x_it, x_l_it); }        // Note: M_B^{-1} is symmetric
    
    template < class ForwardIterator, class OutputIterator >  inline
    void  multiply_transposed_x( ForwardIterator u_l_it,
                                  OutputIterator x_x_it) const
        { multiply_x( u_l_it, x_x_it); }        // Note: M_B^{-1} is symmetric
    
    // vector-vector multiplication ( y = u^T v )
    template < class InputIterator1, class InputIterator2 >  inline
    typename std::iterator_traits<InputIterator1>::value_type
    inner_product( InputIterator1 u_l_it, InputIterator1 u_x_it,
		   InputIterator2 v_l_it, InputIterator2 v_x_it) const
        { return inner_product_l( u_l_it, v_l_it)
	       + inner_product_x( u_x_it, v_x_it); }
    
    template < class InputIterator1, class InputIterator2 >  inline
    typename std::iterator_traits<InputIterator1>::value_type
    inner_product_l( InputIterator1 u_l_it, InputIterator2 v_l_it) const
        { return inner_product( u_l_it, v_l_it, s); }
    
    template < class InputIterator1, class InputIterator2 >  inline
    typename std::iterator_traits<InputIterator1>::value_type
    inner_product_x( InputIterator1 u_x_it, InputIterator2 v_x_it) const
        { return inner_product( u_x_it, v_x_it, b); }
	
    
    // update functions
    // ----------------
    // entering of original variable (update type U1)
    template < class ForwardIterator >
    void  enter_original_( ForwardIterator y_l_it,
			   ForwardIterator y_x_it, const ET& z);
    
    // leaving of original variable (update type U2)
    void  leave_original( );
    
    // entering of slack variable (update type U3)
    void  enter_slack( );

    // leaving of slack variable (update type U4)
    template < class ForwardIterator >
    void  leave_slack_( ForwardIterator u_x_it);
    
    // replacing of original by original variable (update type U5)
    template < class ForwardIterator >
    void  enter_original_leave_original_( ForwardIterator y, unsigned int k);
    
    // replacing of slack by slack variable (update type U6)
    template < class ForwardIterator >
    void  enter_slack_leave_slack_( ForwardIterator u, unsigned int k);
    
    // replacing of slack by original variable (update type U7)
    template < class ForwardIterator1, class ForwardIterator2 >
    void  enter_original_leave_slack_( ForwardIterator1 y, ForwardIterator2 u);
    
    // replacing of original by slack variable (update type U8)
    void  enter_slack_leave_original( );
    
    
    // replacing of original by original variable with precondition in QP-case
    // for phaseII                               (update type UZ_1)
    template < class ForwardIterator >
    void  z_replace_original_by_original(ForwardIterator y_l_it,
                                         ForwardIterator y_x_it,
                                         const ET& s_delta, const ET& s_nu,
					                     unsigned int k_i);


    // replacing of original by slack variable with precondition in QP-case
    // for phaseII                               (update type UZ_2)
    void  z_replace_original_by_slack( );


    // replacing of slack by original variable with precondition in QP-case
    // for phaseII                               (update type UZ_3)
    template < class ForwardIterator >
    void  z_replace_slack_by_original(ForwardIterator y_l_it,
                                      ForwardIterator y_x_it,
				      ForwardIterator u_x_it, const ET& hat_kappa,
				      const ET& hat_nu);


    // replacing of slack by slack variable with precondition in QP-case
    // for phaseII                               (update type UZ_4)
    template < class ForwardIterator >
    void  z_replace_slack_by_slack(ForwardIterator u_x_it, unsigned int k_j);
    

    // copying of reduced basis inverse row in (upper) C-half
    template < class OutIt >
    void  copy_row_in_C(OutIt y_l_it, OutIt y_x_it, unsigned int k);
    
    // copying of reduced basis inverse row in (lower) B_O-half
    template < class OutIt >
    void  copy_row_in_B_O(OutIt y_l_it, OutIt y_x_it, unsigned int k);


    // swap functions
    void  swap_variable( unsigned int j)                // ``to the end'' of R
        { CGAL_qpe_assertion( j < b);
	  swap_variable( j, Is_LP()); }
    void  swap_constraint( unsigned int i)              // ``to the end'' of P
        { CGAL_qpe_assertion( i < s);
	  swap_constraint( i, Is_LP()); }

  private:
    
    // constants
    const ET                 et0, et1, et2;
                                        
    // data members
    Matrix                   M;         // basis inverse, stored row-wise
    ET                       d;         // denominator

    unsigned int             l;         // minimum of `n' and `m'
    unsigned int             s;         // size of `E \cup S_N'
    unsigned int             b;         // size of `B_O'

    bool                     is_phaseI; // flag indicating phase I
    bool                     is_phaseII;// flag indicating phase II
    const bool               is_LP;     // flag indicating a linear    program
    const bool               is_QP;     // flag indicating a quadratic program

    CGAL::Verbose_ostream&   vout;      // used for diagnostic output

    Row                      x_l, tmp_l;  // used in the 
    Row                      x_x, tmp_x;  // update functions

    // private member functions
    // ------------------------
    // set-up
    void  set( Tag_false);        // QP case
    void  set( Tag_true );        // LP case

    // init
    template < class InIt >                                     // QP case
    void  init_( unsigned int art_size, InIt art_first, Tag_false);
    template < class InIt >                                     // LP case
    void  init_( unsigned int art_size, InIt art_first, Tag_true );

    // access
    const ET&  entry( unsigned int row,
		      unsigned int column, Tag_false) const;    // QP case
    const ET&  entry( unsigned int row,
		      unsigned int column, Tag_true ) const;    // LP case

    // matrix-vector multiplication
    template < class ForIt, class OutIt, class Use1stArg >      // QP case
    void  multiply_( ForIt v_l_it, ForIt v_x_it,
		     OutIt y_l_it, OutIt y_x_it, Tag_false, Use1stArg) const;
    template < class ForIt, class OutIt, class  DummyArg >      // LP case
    void  multiply_( ForIt v_l_it, ForIt v_x_it,
		     OutIt y_l_it, OutIt y_x_it, Tag_true,   DummyArg) const;

    // special matrix-vector multiplication functions for LPs
    template < class ForIt, class OutIt >
    void  multiply__l_( ForIt v_x_it, OutIt y_l_it) const;
    template < class ForIt, class OutIt >
    void  multiply__x_( ForIt v_l_it, OutIt y_x_it) const;

    // in-place update
    template < class ForIt >                                    // QP case
    void  update_inplace_QP_( ForIt  y_l_it, ForIt  y_x_it,
			      const ET&  d_new, const ET&  d_old);
    template < class ForIt1, class ForIt2 >                     // LP case
    void  update_inplace_LP_( ForIt1 x_x_it, ForIt2 y_x_it,
			      const ET&  d_new, const ET&  d_old);
			      
    template < class ForIt >                                  // QP case only
    void  z_update_inplace( ForIt psi1_l_it, ForIt psi1_x_it,
                            ForIt psi2_l_it, ForIt psi2_x_it,
			    const ET& omega0, const ET& omega1,
			    const ET& omega2, const ET& omega3); 

    void  update_entry( ET& entry,   const ET& d_new,
			const ET& y, const ET& d_old) const;

    // swap functions
    void  swap_variable  ( unsigned int, Tag_true );            // LP case
    void  swap_variable  ( unsigned int, Tag_false);            // QP case
    void  swap_constraint( unsigned int, Tag_true );            // LP case
    void  swap_constraint( unsigned int, Tag_false);            // QP case	
    
    // diagnostic output
    void  print( );

// ----------------------------------------------------------------------------

// ===============================
// class implementation (template)
// ===============================

  public:

    // creation and initialization
    // ---------------------------
    // init
    template < class InputIterator >
    void
    init( unsigned int art_size, InputIterator art_first)
    {
	CGAL_qpe_assertion_msg( art_size <= l, \
	    "There are more equality constraints than original variables!");

        init( art_size, art_first, Is_LP());
	d = et1;
        CGAL_qpe_assertion( s == art_size);
        CGAL_qpe_assertion( b == art_size);

	is_phaseI  = true;
	is_phaseII = false;

        if ( x_x.size() < art_size) {
            x_x.insert( x_x.end(), art_size-x_x.size(), et0);
        }
	
        if ( tmp_x.size() < art_size) {
            tmp_x.insert( tmp_x.end(), art_size-tmp_x.size(), et0);
        }
        
        CGAL_qpe_debug {
            if ( vout.verbose()) print();
        }
    }

    // transition to phase II
    template < class InputIterator >                            // QP case
    void  transition( InputIterator twice_D_it)
    {
	typename Matrix::iterator  m_it1, m_it2, p_begin, r_begin;
	typename Row   ::iterator  x_it;
	unsigned int      row, col;

	// fill missing rows
	for (row= 0; row< s; ++ row) {
	    CGAL_qpe_assertion(M[row].size()==0);
	    M[row].insert(M[row].end(), row+1, et0);
	}

	// compute new basis inverse [ upper-left part: -(Q^T * 2 D_B * Q) ]
	// -----------------------------------------------------------------
	// compute 'Q^T * 2 D_B' ( Q = A_B^-1 )	
	p_begin = M.begin();
	r_begin = M.begin();
	if (b > 0) r_begin += l+s-1;   	// initialize only if for loops 
					// are entered
	for ( col = 0; col < b; ++col, ++twice_D_it) {
            ++p_begin;

	    // get column of D (D symmetric)
	    std::copy( *twice_D_it, *twice_D_it+b, x_l.begin());

	    // compute 'Q^T * 2 D_Bj'
	    multiply__l( x_l.begin(), x_x.begin());

	    // store resulting column in 'P' and 'R'
	    x_it  = x_x.begin();
	    m_it2 = r_begin;
	    for ( row = l+col; row >= l; --row, --m_it2, ++x_it) {
		(*m_it2)[ row] = *x_it;
	    }
	    m_it1 = p_begin;
	    for ( row = col+1; row <  s; ++row, ++m_it1, ++x_it) {
		(*m_it1)[ col] = *x_it;
	    }
	}

	// compute '(Q^T * 2 D_B) * Q' ( Q = A_B^-1 )
	m_it1 = M.begin();
	m_it2 = r_begin;
	for ( row = 0; row < s; ++row, ++m_it1, --m_it2) {

	    // get row of '(Q^T * 2 D_B)' (stored in 'P' and 'R')
	    std::copy(m_it1->begin()  ,m_it1->begin()+row,    x_l.begin());
	    std::copy(m_it2->begin()+l,m_it2->begin()+(l+b-row),
	    	x_l.begin()+row);

	    // compute '(Q^T * 2 D_B)_i * Q'
	    multiply__l( x_l.begin(), x_x.begin());

	    // negate and store result in 'P'
	    std::transform( x_x.begin(), x_x.begin()+row+1,
			    m_it1->begin(), std::negate<ET>());

	    // clean up in 'R'
	    std::fill_n( m_it2->begin()+l, b-row, et0);
	}

	// scale A_B^-1
	m_it1 = M.begin()+l;
	for ( row = 0; row < s; ++row, ++m_it1) {
	    std::transform( m_it1->begin(), m_it1->begin()+s, m_it1->begin(),
			    boost::bind2nd( std::multiplies<ET>(), d));
	}

	// new denominator: |det(A_B)|^2
	d *= d;

	// update status
	transition();
    }

    // update functions
    // ----------------
    // entering of original variable (update type U1)
    template < class ForwardIterator >
    void
    enter_original( ForwardIterator y_l_it,
                    ForwardIterator y_x_it, const ET& z)
    {
        // assert QP case
        Assert_compile_time_tag( Tag_false(), Is_LP());
    
        // update matrix in-place
        // ----------------------
        // handle sign of new denominator
        CGAL_qpe_assertion( z != et0);
        bool  z_neg = ( z < et0);

        // update matrix
        update_inplace_QP( y_l_it, y_x_it, z, ( z_neg ? -d : d));
                                                                      
        // append new row and column ("after R")
        // -------------------------------------
        typename Row::iterator  row_it;
	ForwardIterator           y_it;
        unsigned int            col, k = l+(++b);
    
//      // allocate new row, if necessary
//      // BG: replaced this by the ensure_physical_row construct below
//      if ( k <= M.size()) {
//           row_it = M[ k-1].begin();
//      } else {
//           row_it = M.insert( M.end(), Row( k, et0))->begin();
//           x_x.push_back( et0);
// 	     tmp_x.push_back( et0);
//      }
	ensure_physical_row(k-1);
	row_it = M[ k-1].begin();
    
        // store entries in new row
        for (   col = 0,                              y_it = y_l_it;
                col < s;
              ++col,     ++row_it,                  ++y_it         ) {
            *row_it = ( z_neg ? *y_it : -( *y_it));
        }
        for (   col = l,     row_it += l-s,           y_it = y_x_it;
                col < k-1;
              ++col,       ++row_it,                ++y_it         ) {
            *row_it = ( z_neg ? *y_it : -( *y_it));
        }
        *row_it = ( z_neg ? -d : d);
    
        // store new denominator
	// ---------------------
        d = ( z_neg ? -z : z);
        CGAL_qpe_assertion( d > et0);
    
        CGAL_qpe_debug {
            if ( vout.verbose()) print();
        }
    }
    
    // leaving of slack variable (update type U4)
    template < class ForwardIterator >
    void
    leave_slack( ForwardIterator u_x_it)
    {
        // assert QP case
        Assert_compile_time_tag( Tag_false(), Is_LP());
    
        // update matrix in-place
        // ----------------------
        // compute new row/column of basis inverse
        multiply( u_x_it,                               // dummy (not used)
		  u_x_it, x_l.begin(), x_x.begin(),
		  Tag_false(),                          // QP
		  Tag_false());                         // ignore 1st argument
        ET    z     = -inner_product_x( x_x.begin(), u_x_it);
        bool  z_neg = ( z < et0);
        CGAL_qpe_assertion( z != et0);
    
        // update matrix
        update_inplace_QP( x_l.begin(), x_x.begin(), z, ( z_neg ? -d : d));
                                                                      
        // insert new row and column ("after P")
        // -------------------------------------
        typename Row   ::iterator  row_it, x_it;
        typename Matrix::iterator  col_it;
        unsigned int               col, k = l+b;
    
        // store entries in new row
	if (M[s].size()==0) {
	   // row has to be filled first
	   M[s].insert(M[s].end(), s+1, et0);
	}
        for (   col = 0,   row_it = M[ s].begin(),        x_it = x_l.begin();
                col < s;
              ++col,     ++row_it,                      ++x_it              ) {
            *row_it = ( z_neg ? *x_it : -( *x_it));
        }
        *row_it = ( z_neg ? -d : d);
    
        for (   col = l,   col_it = M.begin()+l,          x_it = x_x.begin();
                col < k;
              ++col,     ++col_it,                      ++x_it              ) {
            (*col_it)[ s] = ( z_neg ? *x_it : -( *x_it));
        }
        ++s;
    
        // store new denominator
	// ---------------------
        d = ( z_neg ? -z : z);
        CGAL_qpe_assertion( d > et0);
    
        CGAL_qpe_debug {
            if ( vout.verbose()) print();
	}
    }

    // replacing of original variable (update type U5) [ replace column ]
    template < class RandomAccessIterator >
    void
    enter_original_leave_original( RandomAccessIterator y_x_it, unsigned int k)
    {
        // assert LP case or phase I
	CGAL_qpe_assertion( is_LP || is_phaseI);
	CGAL_qpe_assertion( k < b);

        // update matrix in place
        // ----------------------
        typename Matrix::iterator  matrix_it;
        typename Row   ::iterator     row_it, row_k_it, row_k;

        // handle sign of new denominator
        ET    z     = y_x_it[ k];
        bool  z_neg = ( z < et0);
        if ( z_neg) d = -d;

	// QP (in phase I)?
	matrix_it = M.begin();
	if ( is_QP) matrix_it += l;
	row_k = ( matrix_it+k)->begin();

        // rows: 0..s-1 without k
        unsigned int  row, col;
	ET            minus_y;
        for (   row = 0;
                row < s;
              ++row,     ++matrix_it, ++y_x_it) {
	    if ( row != k) {

		// columns: 0..b-1
		minus_y = -( *y_x_it);
		for (   col = 0, row_it = matrix_it->begin(), row_k_it = row_k;
			col < b;
		      ++col,   ++row_it,                    ++row_k_it       ){
        
		    // update in place
		    update_entry( *row_it, z, minus_y * *row_k_it, d);
		}
	    }
        }

	// rows: k (flip signs, if necessary)
	if ( z_neg) {
	    for (   col = 0,   row_it = row_k;
		    col < b;
		  ++col,     ++row_it        ) {
        
		*row_it = -( *row_it);
	    }
	}

        // store new denominator
        // ---------------------
        d = ( z_neg ? -z : z);
        CGAL_qpe_assertion( d > et0);

	// diagnostic output
        CGAL_qpe_debug {
            if ( vout.verbose()) print();
	}
    }
    
    // replacing of slack variable (update type U6) [ replace row ]
    template < class ForwardIterator >
    void
    enter_slack_leave_slack( ForwardIterator u_x_it, unsigned int k)
    {
        // assert LP case or phase I
	CGAL_qpe_assertion( is_LP || is_phaseI);
	CGAL_qpe_assertion( k < s);

        // compute new row of basis inverse
        multiply__l( u_x_it, x_x.begin());

        // update matrix in place
        // ----------------------
        typename Matrix::iterator  matrix_it;
        typename Row   ::iterator     row_it, x_it;

        // handle sign of new denominator
        ET    z     = x_x[ k];
        bool  z_neg = ( z < et0);
        if ( z_neg) d = -d;

	// QP (in phase I)?
	matrix_it = M.begin();
	if ( is_QP) matrix_it += l;

        // rows: 0..s-1
        unsigned int          row, col;
	ET            minus_m_row;
        for (   row = 0;
                row < s;
              ++row,     ++matrix_it) {

	    // columns: 0..b-1
	    minus_m_row = -( *matrix_it)[ k];
	    for (   col = 0,   row_it = matrix_it->begin(), x_it = x_x.begin();
		    col < b;
		  ++col,     ++row_it,                    ++x_it             ){

		if ( col != k) {                // all columns but k

		    // update in place
		    update_entry( *row_it, z, minus_m_row * *x_it, d);

		} else {                        // column k

		    // flip sign, if necessary
		    if ( z_neg) *row_it = -( *row_it);

		}
	    }
	}

        // store new denominator
        // ---------------------
        d = ( z_neg ? -z : z);
        CGAL_qpe_assertion( d > et0);

	// diagnostic output
        CGAL_qpe_debug {
            if ( vout.verbose()) print();
	}
    }
    
    // replacing of slack by original variable (update type U7) [ grow ]
    template < class ForwardIterator1, class ForwardIterator2 >
    void  enter_original_leave_slack( ForwardIterator1 y_x_it,
				      ForwardIterator2 u_x_it)
    {
        // assert LP case or phase I
	CGAL_qpe_assertion( is_LP || is_phaseI);

        // update matrix in-place
        // ----------------------
        // compute new row of basis inverse
        multiply__l( u_x_it, x_x.begin());
        ET    z     = d*u_x_it[ b] - inner_product_x( y_x_it, u_x_it);
        bool  z_neg = ( z < et0);
        CGAL_qpe_assertion( z != et0);
	
        // update matrix
	update_inplace_LP( x_x.begin(), y_x_it, z, ( z_neg ? -d : d));
                                                  
        // append new row and column
        // -------------------------
        typename Matrix::iterator  matrix_it;
        typename Row   ::iterator     row_it, x_it;
        unsigned int                  row, col;

	// QP (in phase I)?
	if ( is_QP) {
	    ensure_physical_row(l+b);
	    row_it = M[l+b].begin();
	    matrix_it = M.begin() + l;
	} else {
	    row_it = M[s].begin();
	    matrix_it = M.begin();
	} 	
	
	// store 'x_x' in new row 
	x_it = x_x.begin();
	for ( col = 0; col < b; ++col, ++row_it, ++x_it) {
		*row_it = ( z_neg ? *x_it : -( *x_it));
	}
        *row_it = ( z_neg ? -d : d);
	
	// store 'y_x' in new col
	for ( row = 0; row < s; ++row, ++matrix_it, ++y_x_it) {
		(*matrix_it)[b] = ( z_neg ? *y_x_it : -( *y_x_it));
	}
	++s; ++b;
    
        // store new denominator
	// ---------------------
        d = ( z_neg ? -z : z);
        CGAL_qpe_assertion( d > et0);
    
        CGAL_qpe_debug {
            if ( vout.verbose()) print();
        }
    }
  // due to basis_matrix_stays_regular fix, needs reconsideration
  //private:

    // private member functions
    // ------------------------
    // init (QP case)
    template < class InIt >  inline
    void
    init( unsigned int art_size, InIt art_first, Tag_false)
    {
	// only Q is used to store A_B^-1 in phase I
        for ( s = l, b = 0; b < art_size; ++s, ++b, ++art_first) {
	    ensure_physical_row(s);
            M[ s][ b] = ( art_first->second ? -et1 : et1);
        }

        s = art_size;
    }

    // init (LP case)
    template < class InIt >  inline
    void
    init( unsigned int art_size, InIt art_first, Tag_true)
    {
        for ( s = 0; s < art_size; ++s, ++art_first) {
	    std::fill_n( M[ s].begin(), art_size, et0);
	    M[ s][ s] = ( art_first->second ? -et1 : et1);
	}
	b = art_size;
    }
    
    // append row in Q if no allocated row available
    void ensure_physical_row (unsigned int row) {
    	unsigned int rows = static_cast<unsigned int>(M.size());
	CGAL_qpe_assertion(rows >= row);
	if (rows == row) {
            M.push_back(Row(row+1, et0));
	    
	    // do we have to grow x_x?
	    CGAL_qpe_assertion(x_x.size() >= row-l);
	    if (x_x.size() == row-l)
	       x_x.push_back(et0);
	    
	    // do we have to grow tmp_x?
	    CGAL_qpe_assertion(tmp_x.size() >= row-l);
	    if (tmp_x.size() == row-l)
	       tmp_x.push_back(et0);
	    
            CGAL_qpe_assertion(M[row].size()==row+1);
	    CGAL_qpe_debug {
	      if ( vout.verbose()) {
                    vout << "physical row " << (row) << " appended in Q\n";
	      }
            }
	}
    }
    
    // matrix-vector multiplication (QP case)
    template < class ForIt, class OutIt, class Use1stArg >
    void
    multiply( ForIt v_l_it, ForIt v_x_it,
              OutIt y_l_it, OutIt y_x_it, Tag_false,
	      Use1stArg use_1st_arg) const
    {
	// use 'LP' functions in phase I
	if ( is_phaseI) {
	    multiply__l( v_x_it, y_l_it);
	    multiply__x( v_l_it, y_x_it);
	    return;
	}

	// phase II
        typename Matrix::const_iterator  matrix_it;
        typename Row   ::const_iterator     row_it;     // left  of diagonal
        typename Matrix::const_iterator  column_it;     // right of diagonal
        ForIt                                 v_it;
    
        unsigned int  row, count, k = l+b;
        ET            sum;
    
        // compute  P v_l + Q^T v_x	
        for (   row = 0,   matrix_it = M.begin();
                row < s;
              ++row,                                                ++y_l_it) {
            sum = et0;

	    // P v_l
	    if ( check_tag( use_1st_arg)) {

		// P: left of diagonal (including)
		for (   row_it =  matrix_it->begin(),            v_it = v_l_it;
			row_it != matrix_it->end();
		      ++row_it,                                ++v_it) {
		    sum += *row_it * *v_it;
		}

		// P: right of diagonal (excluding)
		for (   count = row+1,   column_it  = ++matrix_it;
			count < s;
		      ++count,         ++column_it,                ++v_it) {
		    sum += (*column_it)[ row] * *v_it;
		}
	    }
    
            // Q^T:
            for (   count = 0,       column_it  = M.begin()+l,   v_it = v_x_it;
                    count < b;
                  ++count,         ++column_it,                ++v_it) {
                sum += (*column_it)[ row] * *v_it;
            }
    
            // store result
            *y_l_it = sum;
        }
    
        // compute  Q v_l + R v_x
        for (   row = l,   matrix_it = M.begin()+l;
                row < k;
              ++row,                                                ++y_x_it) {
            sum = et0;

	    // Q v_l
	    if ( check_tag( use_1st_arg)) {

		// Q:
		for (   count = 0,  row_it = matrix_it->begin(), v_it = v_l_it;
			count < s;
		      ++count,    ++row_it,                    ++v_it) {
		    sum += *row_it * *v_it;
		}
	    }
    
	    // R: left of diagonal (including)
            for (                row_it =  matrix_it->begin()+l, v_it = v_x_it;
                                 row_it != matrix_it->end();
                               ++row_it,                       ++v_it) {
                sum += *row_it * *v_it;
            }
    
            // R: right of diagonal (excluding)
            for (   count = row+1,   column_it = ++matrix_it;
                    count < k;
                  ++count,         ++column_it,                ++v_it) {
                sum += (*column_it)[ row] * *v_it;
            }
    
            // store result
            *y_x_it = sum;
        }
    }
    
    // matrix-vector multiplication (LP case)
    template < class ForIt, class OutIt, class Dummy > inline
    void
    multiply( ForIt v_l_it, ForIt v_x_it,
              OutIt y_l_it, OutIt y_x_it, Tag_true, Dummy) const
    {
        multiply__l( v_x_it, y_l_it);
        multiply__x( v_l_it, y_x_it);
    }
    
    // special matrix-vector multiplication functions for LPs
    template < class ForIt, class OutIt > inline
    void
    multiply__l( ForIt v_x_it, OutIt y_l_it) const
    {
        typename Matrix::const_iterator  matrix_it = M.begin();
        typename Matrix::const_iterator  column_it;
        ForIt                                 v_it;
    
        unsigned int  row, count;
        ET            sum;
    
	// QP?
	if ( is_QP) matrix_it += l;

        // compute  Q^T v_x
        for ( row = 0; row < s; ++row,                              ++y_l_it) {
            sum = et0;
    
            for (   count = 0,   column_it = matrix_it,   v_it = v_x_it;
                    count < b;
                  ++count,     ++column_it,             ++v_it) {
                sum += (*column_it)[ row] * *v_it;
            }
    
            *y_l_it = sum;
        }
    }
    
    template < class ForIt, class OutIt > inline
    void
    multiply__x( ForIt v_l_it, OutIt y_x_it) const
    {
        typename Matrix::const_iterator  matrix_it = M.begin();
        unsigned int  row;

	// QP?
	if ( is_QP) matrix_it += l;

        // compute  Q v_l
        for (   row = 0;
                row < b;
              ++row,     ++matrix_it, ++y_x_it) {

	    *y_x_it = inner_product( matrix_it->begin(), v_l_it, s);
	}
    }
    
    // vector-vector multiplication  
    template < class InIt1, class InIt2 > inline
    typename std::iterator_traits<InIt1>::value_type  
    inner_product( InIt1 u_it, InIt2 v_it, unsigned int n) const
    {
	typedef  typename std::iterator_traits<InIt1>::value_type  NT;
    
        // compute u^T v
	NT sum = NT( 0);
        for ( unsigned int count = 0; count < n; ++count, ++u_it, ++v_it) {
            sum += NT(*u_it) * NT(*v_it);
        }
    
        return sum;
    }
    
    // in-place update
    template < class ForIt >                                    // QP case
    void  update_inplace_QP( ForIt y_l_it, ForIt y_x_it,
			     const ET& d_new, const ET& d_old)
    {
        typename Matrix::      iterator  matrix_it;
        typename Row   ::      iterator     row_it;
        typename Row   ::const_iterator      y_it1, y_it2;
    
        unsigned int  row, col, k = l+b;
    
        // rows: 0..s-1  ( P )
        for (   row = 0,   y_it1 = y_l_it,   matrix_it = M.begin();
                row < s;
              ++row,     ++y_it1,          ++matrix_it            ) {
    
            // columns: 0..row  ( P )
            for (   row_it =  matrix_it->begin(),   y_it2 = y_l_it;
                    row_it != matrix_it->end();
                  ++row_it,                       ++y_it2         ) {
    
                update_entry( *row_it, d_new, *y_it1 * *y_it2, d_old);
            }
        }
    
        // rows: l..k-1  ( Q R )
        for (   row = l,   y_it1 = y_x_it,   matrix_it += l-s;
                row < k;
              ++row,     ++y_it1,          ++matrix_it       ) {
    
            // columns: 0..s-1  ( Q )
            for (   col = 0,   row_it =  matrix_it->begin(),   y_it2 = y_l_it;
                    col < s;
                  ++col,     ++row_it,                       ++y_it2         ){
    
                update_entry( *row_it, d_new, *y_it1 * *y_it2, d_old);
            }
    
            // columns: l..k-1  ( R )
            for (              row_it += l-s,                  y_it2 = y_x_it;
                               row_it != matrix_it->end();
                             ++row_it,                       ++y_it2         ){
    
                update_entry( *row_it, d_new, *y_it1 * *y_it2, d_old);
            }
        }
    }
    
    template < class ForIt1, class ForIt2 >                     // LP case
    void  update_inplace_LP( ForIt1 x_x_it, ForIt2 y_x_it,
			     const ET& d_new, const ET& d_old)
    {
        typename Matrix::      iterator  matrix_it;
        typename Row   ::      iterator     row_it;
	ForIt1                                x_it;
    
        unsigned int  row, col;
	ET            y;

	// QP (in phase I)?
	matrix_it = M.begin();
	if ( is_QP) matrix_it += l;

        // rows: 0..s-1  ( Q )
        for (   row = 0;
                row < s;
              ++row,     ++y_x_it, ++matrix_it) {
    
            // columns: 0..b-1  ( Q )
	    y = *y_x_it;
            for (   col = 0,   row_it =  matrix_it->begin(),   x_it = x_x_it;
		    col < b;
		  ++col,     ++row_it,                       ++x_it         ){
    
                update_entry( *row_it, d_new, y * *x_it, d_old);
            }
        }
    }
    
    
    template < class RandomAccessIterator >
    typename std::iterator_traits<RandomAccessIterator>::value_type 
    inv_M_B_row_dot_col( int row, RandomAccessIterator u_l_it) const
    {
	typename Row::const_iterator row_it;
	if ( is_QP) {
	    row_it = M[l + row].begin();
	} else {
	    row_it = M[row].begin();
	}
	return inner_product(row_it, u_l_it, b);	
    }

};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// creation
template < class ET_, class Is_LP_ >  inline
QP_basis_inverse<ET_,Is_LP_>::
QP_basis_inverse( CGAL::Verbose_ostream&  verbose_ostream)
    : et0( 0), et1( 1), et2( 2),
      is_LP( check_tag( Is_LP())), is_QP( ! is_LP),
      vout( verbose_ostream)
{ }

// transition (LP case)
template < class ET_, class Is_LP_ >  inline
void  QP_basis_inverse<ET_,Is_LP_>::
transition( )
{
    is_phaseI  = false;
    is_phaseII = true;

    CGAL_qpe_debug {
	if ( vout.verbose()) print();
    }
}

// set-up (QP case)
template < class ET_, class Is_LP_ >  inline
void  QP_basis_inverse<ET_,Is_LP_>::
set( Tag_false)
{
    M.reserve( l);
    // only allocate empty rows
    for ( unsigned int i = 0; i < l; ++i )
       M.push_back(Row(0, et0)); 
}
    
// set-up (LP case)
template < class ET_, class Is_LP_ >  inline
void  QP_basis_inverse<ET_,Is_LP_>::
set( Tag_true)
{
    M.reserve( l);
    for ( unsigned int i = 0; i < l; ++i) M.push_back( Row( l, et0));
}

// access (QP case)
template < class ET_, class Is_LP_ >  inline
const ET_&  QP_basis_inverse<ET_,Is_LP_>::
entry( unsigned int r, unsigned int c, Tag_false) const
{
    CGAL_qpe_assertion( ( r < s) || ( ( r >= l) && ( r < l+b)));
    CGAL_qpe_assertion( ( c < s) || ( ( c >= l) && ( c < l+b)));
    return ( c < r ? M[ r][ c] : M[ c][ r]);
}

// access (LP case)
template < class ET_, class Is_LP_ >  inline
const ET_&  QP_basis_inverse<ET_,Is_LP_>::
entry( unsigned int r, unsigned int c, Tag_true) const
{
    CGAL_qpe_assertion( r < s);
    CGAL_qpe_assertion( c < b);
    return M[ r][ c];
}

// in-place update
template < class ET_, class Is_LP_ >  inline
void  QP_basis_inverse<ET_,Is_LP_>::
update_entry( ET& entry, const ET& d_new, const ET& y, const ET& d_old) const
{
    entry *= d_new;
    entry += y;
    entry = CGAL::integral_division(entry, d_old);
}

} //namespace CGAL

#include <CGAL/QP_solver/QP_basis_inverse_impl.h>

#endif // CGAL_QP_SOLVER_QP_BASIS_INVERSE_H

// ===== EOF ==================================================================
