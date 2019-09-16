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

namespace CGAL {

// =============================
// class implementation (cont'd)
// =============================

// set-up
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
void  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
set( )
{
    // reserve memory for NT versions of current solution
    //int  l = (std::min)( this->solver().number_of_variables(),
    //		           this->solver().number_of_constraints());
    int l = this->solver().get_l();
    lambda_NT.resize( l, nt0);
    set( l, Is_linear());
}

// initialization; BG: who calls this???
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
void  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
init( )
{
    // get properties of quadratic program
    n = this->solver().number_of_variables();
    int  m = this->solver().number_of_constraints();
    int  s = this->solver().number_of_slack_variables();

    // clear row and column maxima, if necessary
    if ( ! row_max_A.empty()) row_max_A.clear();
    if ( ! handled_A.empty()) handled_A.clear();
    if ( ! col_max  .empty()) col_max  .clear();
    init( Is_linear());

    // initialize row and column maxima
    // assuming nonneg coefficients
    row_max_c = nt0;
    C_auxiliary_iterator v_it;
    for ( v_it = this->solver().c_auxiliary_value_iterator_begin();
              v_it != this->solver().c_auxiliary_value_iterator_end(); ++v_it) {
        if (*v_it > row_max_c) row_max_c = (*v_it);
    }
        
    handled_A.insert( handled_A.end(), m, false);
    row_max_A.insert( row_max_A.end(), m, nt0);

    col_max.insert( col_max.end(), n,   nt0);               // original
    col_max.insert( col_max.end(), s, nt1);               // slack
    							  //auxiliary
    for ( v_it = this->solver().c_auxiliary_value_iterator_begin();
          v_it != this->solver().c_auxiliary_value_iterator_end(); ++v_it) {
	      col_max.insert( col_max.end(), (*v_it));
    }
}

// operations
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
void  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
init_NT( )
{
    // ToDo: scale 'x_B_O', 'lambda', and 'd' if necessary

    // get inexact version of 'lambda'
    std::transform( this->solver().get_lambda_begin(),
		    this->solver().get_lambda_end(),
		    lambda_NT.begin(), et2nt_obj);

    // get inexact version of 'x_B_O'
    init_NT( Is_linear());

    // get inexact version of 'd'
    d_NT = et2nt_obj( this->solver().variables_common_denominator());
}

// Here is how it works (Sven's thesis, page 99):
// ----------------------------------------------
// R_0   := max_j |c_j|   (largest objective function coefficient)
// R_i^A := max_j |A_ij|  (largest entry in row i of A)
// R_i^D := max_j 2|D_ij| (largest entry in row i of D)
// C_j   := max_j (|c_j|, max_i |A_ij|, max_i |D_ij|)
//
// U := max(
//          d_NT * R_0, 
//          max_{i basic} |lambda_NT[i]| * R_i^A,
//          max_{j basic} |x_NT[j]| * R_j^D
//         )
// W := max(
//          d_NT, 
//          max_{i basic} |lambda_NT[i]|,
//          max_{j basic} |x_NT[j]|
//         )
// Note: all max_{j basic} are over basic ORIGINAL variables only
//
// q := (1+1/64) * 
//      (#basic constraints + #basic original variables + 1) *
//      (#basic constraints + #basic original variables + 2) * 2^(-53)
//
// the actual error bound for mu_j_NT is then
//
// min (U * q, W * q * C_j)
//
// the code maintains and manipulates:
// 
// row_max_c    = R_0   (set in init() for phase I, transition() for phase II)
// row_max_A[i] = R_i^A (set in update_maxima())
// row_max_D[i] = R_i^D (set in update_maxima(Tag_false))
// col_max[j]   = C_j   (set in update_maxima() + Tag_false variant)
// bound1       = U * q (computed in update_maxima() + Tag_true/Tag_false )
// bound2_wq    = W * q (computed in update_maxima() + Tag_true/Tag_false )
//
// Changes for upper-bounded case:
// -------------------------------
// in the scalar product considered on p.100, we have another term:
// d * w_j, where w_j is the correction term 2 x_N^T D_Nj that shows
// up additionally in the pricing and that the QP_solver maintains in 
// the vector w. The first bound on p.101 yields U, while the second
// bound yields W. This means that
// - C_j = max_y|y_i| has an additional term |w_NT[j]| in the maximum,
// - U has an additional term d_NT * R_w in the maximum, where
// 
//     R_w = max_j |w_j|
//
// The problem is that R_w is not a static quantity like max_j|c_j|,
// so it has to be handled per j in certify_mu_j_NT:
// --> for q: (c+b)*(c+b+1) -> (c+b+1)*(c+b+2) 
// --> for bound1: take d_NT * |w_NT[j]| into account for U
// --> for bound2: take |w_NT[j]| into account for C_j


template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
void  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
update_maxima( )
{
    // get properties of quadratic program
    int  c = this->solver().number_of_basic_constraints();
    int  b = this->solver().number_of_basic_original_variables();

    // initialize error bounds
    bound1    = d_NT * row_max_c;
    bound2_wq = d_NT;

    // update row and column maxima of 'A'
    A_iterator  a_it = this->solver().a_begin();
    R_iterator  r_it = this->solver().row_type_begin();
    int         row, col;
    NT          row_max, z;

    Basic_constraint_index_iterator  it;
    Values_NT_iterator v_it = lambda_NT.begin();
    for ( it =  this->solver().basic_constraint_indices_begin();
	  it != this->solver().basic_constraint_indices_end(); ++it, ++v_it) {
	row = *it;

	// row not handled yet?
	if ( ! handled_A[ row]) {

	    // slack variable (or artificial in phase I) involved?
	    row_max = (    ( r_it[ row] != CGAL::EQUAL)
			|| ( this->solver().phase() == 1) ? nt1 : nt0);

	    // scan row and update maxima
	    for ( col = 0; col < n; ++col) {
		z = CGAL::abs( *((*(a_it + col))+row));
		if ( z > row_max      ) row_max       = z;
		if ( z > col_max[ col]) col_max[ col] = z;
	    }
	    row_max_A[ row] = row_max;
	    handled_A[ row] = true;
	}
	// update bounds
	z = CGAL::abs( *v_it);
	if ( z > bound2_wq) bound2_wq = z;
	z *= row_max_A[ row];
	if ( z > bound1) bound1 = z;
    }

    // update row and column maxima of 'D', if necessary
    if (this->solver().phase() == 1) {
    	update_maxima( Tag_true());
    } else {
        update_maxima( Is_linear());
    }

    // finalize error bounds
    // ToDo: use std::numeric_limits for 'machine epsilon' to
    // support types different from double
    set_q(c, b);
    bound1    *= q;
    bound2_wq *= q;

    CGAL_qpe_debug {
	this->vout() << std::endl
	       << "first bound for certification: " << bound1 << std::endl;
    }
}

template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >                // QP case
void  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
update_maxima( Tag_false)
{
    // update row and column maxima of 'D'
    D_row_iterator  d_row_it;
    int             row, col;
    NT              row_max, z;

    Basic_variable_index_iterator  it;
    Values_NT_iterator v_it = x_B_O_NT.begin();
    for ( it =  this->solver().basic_original_variable_indices_begin();
	  it != this->solver().basic_original_variable_indices_end();
	  ++it, ++v_it) {
	row = *it;

	// row not handled yet?
	if ( ! handled_D[ row]) {
	    d_row_it = this->solver().d_begin()[ row];

	    // scan row and update maxima
	    row_max = nt0;
	    for ( col = 0; col < n; ++col, ++d_row_it) {
		z = CGAL::abs( *d_row_it);
		if ( z > row_max      ) row_max       = z;
		if ( z > col_max[ col]) col_max[ col] = z;
	    }
	    row_max_D[ row] = row_max;
	    handled_D[ row] = true;
	}
	// update bounds
	z = CGAL::abs( (*v_it));
	if ( z > bound2_wq) bound2_wq = z;
	z *= row_max_D[ row];
	if ( z > bound1) bound1 = z;
    }
}

template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >                // LP case
void  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
update_maxima( Tag_true)
{

    // update basic original variables maxima
    NT              z;

    Values_NT_iterator v_it;
    for ( v_it = x_B_O_NT.begin(); v_it != x_B_O_NT.end(); ++v_it) {

	// update bounds
	z = CGAL::abs( (*v_it));
	if ( z > bound2_wq) bound2_wq = z;
    }
}

template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
bool  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
certify_mu_j_NT( int j, Tag_true) const  // standard form
{
    // compute 'mu_j' with inexact arithmetic again
    NT  mu = mu_j_NT( j);

    CGAL_qpe_debug {
	this->vout() << "mu_" << j << " [NT]: " << mu;
    }

    // check against first bound
    if ( mu >= bound1) {
	CGAL_qpe_debug {
	    this->vout() << "  ok [1]" << std::endl;
	}
	return true;
    }

    // compute and check against second bound
    NT  bound2 = bound2_wq * col_max[ j];
    if ( mu >= bound2) {
	CGAL_qpe_debug {
	    this->vout() << "  ok [2: " << bound2 << "]" << std::endl;
	}
	return true;
    }

    // compute and check exact 'mu_j'
    ET  mu_et = this->mu_j( j);

    CGAL_qpe_debug {
	this->vout() << "  " << ( mu_et >= this->et0 ? "ok" : "MISSED")
	       << " [exact: " << mu_et << "]" << std::endl;
    }
    return ( mu_et >= this->et0);
}

template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
bool  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
certify_mu_j_NT( int j, Tag_false) const // case with bounds
{
  // compute 'mu_j' with inexact arithmetic again; this call sets w_j_NT
  NT mu = CGAL::abs(mu_j_NT( j)); // may have both signs in bounded case
  NT w_j_NT_abs = CGAL::abs(w_j_NT);

  CGAL_qpe_debug {
      this->vout() << "|mu_|" << j << " [NT]: " << mu;
  }

  // check against first bound, taking w_j_NT into account
  // q has been set in previous call to update_maxima()
  // bound1_with_w = (U + d_NT*|w_j|) * q
  NT bound1_with_w = bound1 + (d_NT * w_j_NT_abs * q); 
  if ( mu >= bound1_with_w) {
    CGAL_qpe_debug {
      this->vout() << "  ok [1]" << std::endl;
    }
    return true;
  }

  // compute and check against second bound, taking w_j_NT into account
  // bound2_with_w = max (W, |w_j|)
  NT  bound2_with_w = bound2_wq * (std::max)(col_max[ j], w_j_NT_abs);
  if ( mu >= bound2_with_w) {
    CGAL_qpe_debug {
      this->vout() << "  ok [2: " << bound2_with_w << "]" << std::endl;
    }
    return true;
  }

  // the bounds are insufficient: compute and check exact 'mu_j'
  ET  mu_et = this->mu_j( j);
  bool improving = is_improving(j, mu_et, this->et0);

  CGAL_qpe_debug {
    this->vout() << "  " << (!improving ? "ok" : "MISSED")
		 << " [exact: " << mu_et << "]" << std::endl;
  }
  return ( !improving);
}

// transition
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
void  QP__filtered_base<Q,ET,Tags,NT_,ET2NT_>::
transition( )
{

    // update row and column maxima with original objective vector 'c'
    C_iterator  c_it = this->solver().c_begin();
    NT  z;
    row_max_c = nt0;
    for ( int i = 0; i < n; ++i, ++c_it) {
	z = CGAL::abs( *c_it);
	if ( z > row_max_c  ) row_max_c   = z;
	if ( z > col_max[ i]) col_max[ i] = z;
    }

    // initialize row maxima of 'D', if necessary
    transition( n, Is_linear());
}

} //namespace CGAL

// ===== EOF ==================================================================
