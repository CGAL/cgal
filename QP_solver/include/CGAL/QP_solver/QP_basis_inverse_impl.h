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

// creation and initialization
// ---------------------------
// set-up
template < class ET_, class Is_LP_ >
void  QP_basis_inverse<ET_,Is_LP_>::
set( int n, int m, int nr_equalities)
{
    CGAL_qpe_assertion( n >= 0);
    CGAL_qpe_assertion( m >= 0);
    b = s = 0;
    // l is the maximum size of the basis in phase I
    l = (std::min)( n+nr_equalities+1, m);
    if ( ! M.empty()) M.clear();
    set( Is_LP());
    
    if ( ! x_l.empty()) x_l.clear();
    if ( ! x_x.empty()) x_x.clear();
       
    x_l.insert( x_l.end(), l, et0);
    x_x.insert( x_x.end(), l, et0); // has to grow later QP-case
    
    if ( ! tmp_l.empty()) tmp_l.clear();
    if ( ! tmp_x.empty()) tmp_x.clear();

    tmp_l.insert( tmp_l.end(), l, et0);
    tmp_x.insert( tmp_x.end(), l, et0); // has to grow later QP-case

}

// update functions
// ----------------
// leaving of original variable (update type U2)
template < class ET_, class Is_LP_ >
void  QP_basis_inverse<ET_,Is_LP_>::
leave_original( )
{
    // assert QP case
    Assert_compile_time_tag( Tag_false(), Is_LP());

    // determine new denominator (`z')
    --b;
    ET    z     = M[ l+b][ l+b];
    bool  z_neg = ( z < et0);
    CGAL_qpe_assertion( z != et0);

    // update matrix in place
    update_inplace_QP( M[ l+b].begin(), M[ l+b].begin()+l,
		       -z, ( z_neg ? d : -d));
                                                                 
    // store new denominator
    d = ( z_neg ? -z : z);
    CGAL_qpe_assertion( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
}

// entering of slack variable (update type U3)
template < class ET_, class Is_LP_ >
void  QP_basis_inverse<ET_,Is_LP_>::
enter_slack( )
{
    // assert QP case
    Assert_compile_time_tag( Tag_false(), Is_LP());

    // determine new denominator (`z')
    --s;
    ET    z     = M[ s][ s];
    bool  z_neg = ( z < et0);
    CGAL_qpe_assertion( z != et0);

    // update matrix in place
    typename Matrix::iterator  col_it;
    typename Row   ::iterator    x_it;
    unsigned int               col;
    for (   col = 0,   col_it = M.begin()+l,   x_it = x_x.begin();
            col < b;
          ++col,     ++col_it,               ++x_it              ) {
        *x_it = (*col_it)[ s];
    }
    update_inplace_QP( M[ s].begin(), x_x.begin(), -z, ( z_neg ? d : -d));

    // store new denominator
    d = ( z_neg ? -z : z);
    CGAL_qpe_assertion( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
}

// replacing of original by slack variable (update type U8)
template < class ET_, class Is_LP_ >
void  QP_basis_inverse<ET_,Is_LP_>::
enter_slack_leave_original( )
{
    // assert LP case or phase I
    CGAL_qpe_assertion( is_LP || is_phaseI);

    // update matrix in-place
    // ----------------------
    typename Matrix::iterator  matrix_it;
    typename Row   ::iterator       x_it;
    unsigned int                     row;

    // QP (in phase I)?
    matrix_it = M.begin();
    if ( is_QP) matrix_it += l;

    // get last column of basis inverse (store it in 'x_x')
    --s; --b;
    for (   row = 0,   x_it = x_x.begin();
	    row < s;
	  ++row,     ++x_it,               ++matrix_it) {
	*x_it = (*matrix_it)[ b];
    }
    ET    z     = (*matrix_it)[ b];
    bool  z_neg = ( z < et0);
    CGAL_qpe_assertion( z != et0);

    // update matrix
    update_inplace_LP( matrix_it->begin(), x_x.begin(), -z, ( z_neg ? d : -d));

    // store new denominator
    // ---------------------
    d = ( z_neg ? -z : z);
    CGAL_qpe_assertion( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
}


// replacing of original by original variable with precondition in QP-case
// for phaseII                               (update type UZ_1)
template < class ET_, class Is_LP_ >
template < class ForwardIterator >
void  QP_basis_inverse<ET_,Is_LP_>::
z_replace_original_by_original(ForwardIterator y_l_it,
                               ForwardIterator y_x_it, const ET& s_delta,
                               const ET& s_nu, unsigned int k_i)
{

    // assert QP case and phaseII
    CGAL_qpe_assertion(is_QP && is_phaseII);


    // prepare \hat{k}_{1} -scalar
    ET  hat_k_1 = *(y_x_it + k_i);

    // prepare \hat{\rho} -vector in x_l, x_x
    copy_row_in_B_O(x_l.begin(), x_x.begin(), k_i);
    
    // prepare \hat{v} -vector in tmp_l, tmp_x
    
    // tmp_l -part
    std::transform(y_l_it, (y_l_it+s), x_l.begin(), tmp_l.begin(),
        compose2_2(std::plus<ET>(), Identity<ET>(),
        boost::bind1st(std::multiplies<ET>(), s_delta)));
    
    // tmp_x -part    
    std::transform(y_x_it, (y_x_it+b), x_x.begin(), tmp_x.begin(),
        compose2_2(std::plus<ET>(), Identity<ET>(),
        boost::bind1st(std::multiplies<ET>(), s_delta)));
    tmp_x[k_i] -= d;
    
    // prepare \hat{k}_{2} -scalar
    ET  hat_k_2 = s_nu - (et2 * s_delta * hat_k_1);
    
    CGAL_qpe_assertion( d != et0);
        
    // update matrix in place
    z_update_inplace(x_l.begin(), x_x.begin(), tmp_l.begin(), tmp_x.begin(),
                      hat_k_1 * hat_k_1, -hat_k_2, -hat_k_1, d*d);
    
    // store new denominator
    d = CGAL::integral_division(hat_k_1 * hat_k_1, d);

    CGAL_qpe_assertion( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
    
}


// replacing of original by slack variable with precondition in QP-case
// for phaseII                               (update type UZ_2)
template < class ET_, class Is_LP_ >
void  QP_basis_inverse<ET_,Is_LP_>::
z_replace_original_by_slack( )
{

    // assert QP case and phaseII
    CGAL_qpe_assertion(is_QP && is_phaseII);

    // adapt s and b
    --s; --b;

    // prepare \hat{\rho} -vector in x_l, x_x
    copy_row_in_B_O(x_l.begin(), x_x.begin(), b);
    
    // prepare \hat{\varrho} -vector in tmp_l, tmp_x
    copy_row_in_C(tmp_l.begin(), tmp_x.begin(), s);
    
    // prepare \hat{\kappa} -scalar
    ET  hat_kappa = M[l+b][s];
    
    // prepare \hat{\xi} -scalar
    ET hat_xi = M[s][s];
        
    CGAL_qpe_assertion( d != et0);
    
    // update matrix in place
    z_update_inplace(x_l.begin(), x_x.begin(), tmp_l.begin(), tmp_x.begin(),
                           hat_kappa * hat_kappa, hat_xi, -hat_kappa, d * d);
		     
    // store new denominator
    d = CGAL::integral_division(hat_kappa * hat_kappa, d);

    CGAL_qpe_assertion( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }

}


// replacing of slack by original variable with precondition in QP-case
// for phaseII                               (update type UZ_3)
template < class ET_, class Is_LP_ >
template < class ForwardIterator >
void  QP_basis_inverse<ET_,Is_LP_>::
z_replace_slack_by_original(ForwardIterator y_l_it,
                            ForwardIterator y_x_it,
			                ForwardIterator u_x_it, const ET& hat_kappa,
		                    const ET& hat_nu)
{
    // assert QP case and phaseII
    CGAL_qpe_assertion(is_QP && is_phaseII);
    
    // get copies of y_l_it and y_x_it for later use
    ForwardIterator y_l_it_copy = y_l_it;
    ForwardIterator y_x_it_copy = y_x_it;

    CGAL_qpe_assertion( d != et0);
    
    // prepare \hat{\phi}
     
    // prepare \hat{\varphi} -vector in x_l, x_x
    multiply(u_x_it, u_x_it, x_l.begin(), x_x.begin(), Tag_false(),
             Tag_false());
	     
    // prepare \hat{\kappa} -scalar
    
    // prepare \hat{\nu} -scalar
   
    // update matrix in place
    z_update_inplace(x_l.begin(), x_x.begin(), y_l_it, y_x_it,
                     hat_kappa * hat_kappa, -hat_nu, hat_kappa, d * d);    
    
    // append new rows and columns
    // ---------------------------
    typename Row   ::iterator  row_it, x_l_it, x_x_it;
    typename Matrix::iterator  matrix_it;
    unsigned int               count;
    
    // insert new row and column at the end of block P
    CGAL_qpe_assertion(M.size()>=s+1);
    if (M[s].size()==0) {
	// row has to be filled first
        M[s].insert(M[s].end(), s+1, et0);
    }
     
    
    // P-block: left of diagonal (including element on diagonal)
    y_l_it = y_l_it_copy;
    for (  row_it = M[s].begin(), x_l_it = x_l.begin();
           row_it != M[s].end() - 1;
	 ++row_it,  ++x_l_it,  ++y_l_it                ) {
        *row_it = 
	  CGAL::integral_division((hat_nu * *x_l_it)-(hat_kappa * *y_l_it), d);  
    } 
    *row_it = -hat_nu;
     
    // Q-block
    y_x_it = y_x_it_copy;
    for (  matrix_it = M.begin()+l, count = 0, x_x_it = x_x.begin();
           count < b;
	 ++matrix_it,  ++count, ++x_x_it, ++y_x_it                  ) {
        (*matrix_it)[s] = 
	  CGAL::integral_division((hat_nu * *x_x_it) - (hat_kappa * *y_x_it), d);
    }
          
    // insert new row and column at the end of blocks Q and R
    ensure_physical_row(l+b);
    
    // Q-block
    for (  row_it = M[l+b].begin(), count = 0, x_l_it = x_l.begin();
           count < s;
	 ++row_it,  ++count,  ++x_l_it                              ) {
        *row_it = CGAL::integral_division(-hat_kappa * *x_l_it, d);
    }
    *row_it = hat_kappa;
    
    // R-block
    for (  row_it = M[l+b].begin()+l, count = 0, x_x_it = x_x.begin();
           count < b;
	 ++row_it,  ++count,  ++x_x_it                                ) {
        *row_it = CGAL::integral_division(-hat_kappa * *x_x_it, d);
    }
    *row_it = et0;
    
    //adapt s and b
    ++s; ++b; 

    // store new denominator
    d = CGAL::integral_division(hat_kappa * hat_kappa, d);

    CGAL_qpe_assertion( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
		     
}


// replacing of slack by slack variable with precondition in QP-case
// for phaseII                               (update type UZ_4)
template < class ET_, class Is_LP_ >
template < class ForwardIterator >
void  QP_basis_inverse<ET_,Is_LP_>::
z_replace_slack_by_slack(ForwardIterator u_x_it, unsigned int k_j)
{

    // assert QP case and phaseII
    CGAL_qpe_assertion(is_QP && is_phaseII);

    // prepare \hat{v} -vector in x_l, x_x
    multiply(u_x_it, u_x_it, x_l.begin(), x_x.begin(),Tag_false(),
             Tag_false());
    x_l[k_j] -= d;
    
    // prepare \hat{\varrho} -vector in tmp_l, tmp_x
    copy_row_in_C(tmp_l.begin(), tmp_x.begin(), k_j);
    
    // prepare \hat{k}_{1} -scalar
    ET  hat_k_1 = inner_product_x(tmp_x.begin(), u_x_it);
    
    // prepare \hat{k}_{3} -scalar
    ET  hat_k_3 = -M[k_j][k_j];
    
    CGAL_qpe_assertion( d != et0);    
    
    // update matrix in place
    z_update_inplace(x_l.begin(), x_x.begin(), tmp_l.begin(), tmp_x.begin(),
                     hat_k_1 * hat_k_1, -hat_k_3, -hat_k_1, d * d);
		     
    // store new denominator
    d = CGAL::integral_division(hat_k_1 * hat_k_1, d);

    CGAL_qpe_assertion( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
      
}


// copying of reduced basis inverse row in (upper) C-half
template < class ET_, class Is_LP_ >
template < class OutIt >
void  QP_basis_inverse<ET_,Is_LP_>::
copy_row_in_C(OutIt y_l_it, OutIt y_x_it, unsigned int r)
{
    typename Matrix::const_iterator  matrix_it;
    typename Row   ::const_iterator     row_it;
    unsigned int  count;
    
    // P-block: left of diagonal (including element on diagonal)
    matrix_it = M.begin()+r;
    for (  row_it = matrix_it->begin();
           row_it != matrix_it->end(); 
	 ++row_it, ++y_l_it            ) {
        *y_l_it = *row_it;    
    }
    
    // P-block: right of diagonal (excluding element on diagonal)
    for (  matrix_it = M.begin()+r+1, count = r+1;
           count < s; 
	 ++matrix_it,  ++count,  ++y_l_it         ) {
        *y_l_it = (*matrix_it)[r];
    }
    
    // Q-block
    for (  matrix_it = M.begin()+l, count = 0;
           count < b;
	 ++matrix_it,  ++count,  ++y_x_it     ) {
	*y_x_it = (*matrix_it)[r]; 
    } 
}


// copying of reduced basis inverse row in (lower) B_O-half
template < class ET_, class Is_LP_ >
template < class OutIt >
void  QP_basis_inverse<ET_,Is_LP_>::
copy_row_in_B_O(OutIt y_l_it, OutIt y_x_it, unsigned int r)
{
    typename Matrix::const_iterator  matrix_it;
    typename Row   ::const_iterator     row_it;
    unsigned int  count;
    
    // Q-block
    matrix_it = M.begin()+l+r;
    for (  row_it = matrix_it->begin(), count = 0;
           count < s;
	 ++row_it,  ++count,  ++y_l_it           ) {
        *y_l_it = *row_it;
    }
    
    // R-block: left of diagonal (including element on diagonal)
    for (  row_it = matrix_it->begin()+l; 
           row_it != matrix_it->end();
	 ++row_it,  ++y_x_it            ) {
        *y_x_it = *row_it;
    }
    
    // R-block: right of diagonal (excluding element on diagonal)
    for (  matrix_it = M.begin()+l+r+1, count = r+1;
           count < b;
	 ++matrix_it,  ++count,  ++y_x_it           ) {
        *y_x_it = (*matrix_it)[l+r];
    }

}

template < class ET_, class Is_LP_ >
template < class ForIt >
void  QP_basis_inverse<ET_,Is_LP_>::
z_update_inplace( ForIt psi1_l_it, ForIt psi1_x_it,
                  ForIt psi2_l_it, ForIt psi2_x_it,
	             const ET& omega0, const ET& omega1,
		         const ET& omega2, const ET& omega3)
{
    typename Matrix::      iterator  matrix_it;
    typename Row   ::      iterator     row_it;
    typename Row   ::const_iterator      y_it1_r, y_it1_c, y_it2_r, y_it2_c;
	
    unsigned int  row, col, k = l+b;
    ET           u_elem;

    // rows: 0..s-1  ( P )
    for (  row = 0, matrix_it = M.begin(),
           y_it1_r = psi1_l_it,  y_it2_r = psi2_l_it;
	   row < s;
         ++row, ++matrix_it, ++y_it1_r, ++y_it2_r  ) {
	      
        // columns: 0..row  ( P )
        for (   row_it =  matrix_it->begin(),
	        y_it1_c = psi1_l_it,  y_it2_c = psi2_l_it;
                row_it != matrix_it->end();
              ++row_it,  ++y_it1_c,  ++y_it2_c            ) {
                
            u_elem = (*y_it1_r * *y_it2_c) + (*y_it2_r * *y_it1_c);
	    u_elem *= omega2;
	    u_elem += omega1 * *y_it1_r * *y_it1_c;
            update_entry( *row_it, omega0, u_elem, omega3);
        } 
    }
	
    // rows: l..k-1  ( Q R )
    for (  row = l, matrix_it = M.begin()+l,
	   y_it1_r = psi1_x_it,  y_it2_r = psi2_x_it;
	   row != k;
	 ++row,  ++matrix_it,  ++y_it1_r,  ++y_it2_r ) {
	    
        // columns: 0..s-1  ( Q )
        for (   col = 0,   row_it =  matrix_it->begin(),
	        y_it1_c = psi1_l_it,  y_it2_c = psi2_l_it;
                col < s;
              ++col, ++row_it,  ++y_it1_c,  ++y_it2_c     ){
    
            u_elem = (*y_it1_r * *y_it2_c) + (*y_it2_r * *y_it1_c);
	       u_elem *= omega2;
	       u_elem += omega1 * *y_it1_r * *y_it1_c; 
	       update_entry( *row_it, omega0, u_elem, omega3);
        }
    
        // columns: l..k-1  ( R )
        for (  row_it = matrix_it->begin()+l,
	       y_it1_c = psi1_x_it,  y_it2_c = psi2_x_it;
               row_it != matrix_it->end();
             ++row_it,  ++y_it1_c,  ++y_it2_c            ){
		 
            u_elem = (*y_it1_r * *y_it2_c) + (*y_it2_r * *y_it1_c);
            u_elem *= omega2;
	        u_elem += omega1 * *y_it1_r * *y_it1_c;     
            update_entry( *row_it, omega0, u_elem, omega3);
        }
	    
    } 
	
} 


// swap functions
// --------------
// swap variable ``to the end'' of R
template < class ET_, class Is_LP_ >                            // LP case
void  QP_basis_inverse<ET_,Is_LP_>::
swap_variable( unsigned int j, Tag_true)
{
    unsigned int  k = b-1;
    if ( j == k) return;

    // swap rows
    // ---------
    typename Row::iterator   row_j_it = M[ j].begin();
    typename Row::iterator   row_k_it = M[ k].begin();
    unsigned int  count;

    // swap entries 0..b-1 (row <-> row) [in Q]
    for (   count = 0;
            count < b;
          ++count,     ++row_j_it, ++row_k_it) {
        std::iter_swap( row_j_it, row_k_it);
    }
}

template < class ET_, class Is_LP_ >                            // QP case
void  QP_basis_inverse<ET_,Is_LP_>::
swap_variable( unsigned int j, Tag_false)
{
    unsigned int  i = l+j, k = l+b-1;
    if ( i == k) return;

    // swap rows and columns
    // ---------------------
    typename    Row::iterator   row_i_it = M[ i].begin();
    typename    Row::iterator   row_k_it = M[ k].begin();
    typename Matrix::iterator  matrix_it = M.begin()+(i+1);
    unsigned int  count;

    // swap entries 0..s-1 (row <-> row) [in Q]
    for (   count = 0;
            count < s;
          ++count,     ++row_i_it, ++row_k_it) {
        std::iter_swap( row_i_it, row_k_it);
    }

    if ( is_phaseII) {

	// swap entries l..i-1 (row <-> row) [in R]
	for (   count = l,   row_i_it += l-s,   row_k_it += l-s;
		count < i;
	      ++count,     ++row_i_it,        ++row_k_it       ) {
	    std::iter_swap( row_i_it, row_k_it);
	}

	// swap entries i+1..k-1 (column <-> row) [in R]
	for ( ++count,                        ++row_k_it;
		count < k;
	      ++count,     ++matrix_it,       ++row_k_it) {
	    std::swap( ( *matrix_it)[ i], *row_k_it);
	}

	// swap entries i,i with k,k (entry <-> entry) [in R]
	std::iter_swap( row_i_it, row_k_it);
    }
}

// swap constraint ``to the end'' of P
template < class ET_, class Is_LP_ >                            // LP case
void  QP_basis_inverse<ET_,Is_LP_>::
swap_constraint( unsigned int i, Tag_true)
{
    unsigned int  k = s-1;
    if ( i == k) return;

    // swap columns
    // ------------
    typename Matrix::iterator  matrix_it = M.begin();
    unsigned int  count;

    // swap entries 0..s-1 (column <-> column) [in Q]
    for (   count = 0;
            count < s;
          ++count,     ++matrix_it) {
        std::swap( ( *matrix_it)[ i], ( *matrix_it)[ k]);
    }
}

template < class ET_, class Is_LP_ >                            // QP case
void  QP_basis_inverse<ET_,Is_LP_>::
swap_constraint( unsigned int i, Tag_false)
{
 
    if ( i == s-1) return;

    // swap rows and columns
    // ---------------------
    typename    Row::iterator   row_i_it = M[ i].begin();
    typename    Row::iterator   row_k_it = M[ s-1].begin();
    typename Matrix::iterator  matrix_it = M.begin()+i;
    unsigned int  count;

    if ( is_phaseI) {

	// skip empty P
	matrix_it =M.begin() + l;

    } else {

	// swap entries 0..i-1 (row <-> row) [in P]
	for (   count =  0;
		count < i;
	      ++count,      ++row_i_it, ++row_k_it) {
	    std::iter_swap( row_i_it, row_k_it);
	}

	// swap entries i+1..s-2 (column <-> row) [in P]
	for ( count = i + 1, ++matrix_it, ++row_k_it;
		count < s-1;
	      ++count,     ++matrix_it, ++row_k_it) {
	    std::swap( ( *matrix_it)[ i], *row_k_it);
	}
	// the remaining two entries to be swapped on the main diagonal
	std::swap(M[i][i], M[s-1][s-1]);

	// advance to Q
	matrix_it = M.begin() + l;
    }

    // swap entries l..l+b (column <-> column) [in Q]
    for (   count = 0;
            count < b;
          ++count,     ++matrix_it) {
        std::swap( ( *matrix_it)[ i], ( *matrix_it)[ s-1]);
    }
}

// diagnostic output
// -----------------
template < class ET_, class Is_LP_ >
void  QP_basis_inverse<ET_,Is_LP_>::
print( )
{
    // P
    if ( is_LP || is_phaseII) {
	for ( unsigned int row = 0; row < s; ++row) {
	    std::copy( M[ row].begin(),
		       M[ row].begin() + ( is_LP ? s : row+1),
		       std::ostream_iterator<ET>( vout.out(), " "));
	    vout.out() << std::endl;
	}
	if ( is_QP) vout.out() << "= = = = = = = = = =" << std::endl;
    }

    // Q & R
    if ( is_QP) {
	for ( unsigned int row = l; row < l+b; ++row) {
	    std::copy( M[ row].begin(), M[ row].begin()+s,
		       std::ostream_iterator<ET>( vout.out(), " "));
	    if ( is_phaseII) {
		vout.out() << "|| ";
		std::copy( M[ row].begin()+l, M[ row].end(),
			   std::ostream_iterator<ET>( vout.out(), " "));
	    }
	    vout.out() << std::endl;
	}
    }
    vout.out() << "denominator = " << d << std::endl;
}

} //namespace CGAL

// ===== EOF ==================================================================
