// ============================================================================
//
// Copyright (c) 1997-2004 The CGAL Consortium
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
// file          : include/CGAL/QP_engine/QPE__filtered_base.C
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Base Class for Filtered Pricing of the QPE Solver
// ============================================================================

CGAL_BEGIN_NAMESPACE

// =============================
// class implementation (cont'd)
// =============================

// set-up
template < class Rep_, class NT_, class ET2NT_ >
void  QPE__filtered_base<Rep_,NT_,ET2NT_>::
set( )
{
    // reserve memory for NT versions of current solution
    int  l = std::min( solver().number_of_variables(),
		       solver().number_of_constraints());
    lambda_NT.resize( l, nt0);
    set( l, Is_linear());
}

// initialization
template < class Rep_, class NT_, class ET2NT_ >
void  QPE__filtered_base<Rep_,NT_,ET2NT_>::
init( )
{
    // get properties of quadratic program
    int  n = solver().number_of_variables();
    int  m = solver().number_of_constraints();
    int  s = solver().number_of_slack_variables();
    int  a = solver().number_of_artificial_variables();

    // clear row and column maxima, if necessary
    if ( ! row_max_A.empty()) row_max_A.clear();
    if ( ! handled_A.empty()) handled_A.clear();
    if ( ! col_max  .empty()) col_max  .clear();
    init( Is_linear());

    // initialize row and column maxima
    row_max_c = nt1;

    handled_A.insert( handled_A.end(), m, false);
    row_max_A.insert( row_max_A.end(), m, nt0);

    col_max.insert( col_max.end(), n,   nt0);               // original
    col_max.insert( col_max.end(), s+a, nt1);               // slack+artificial
}

// operations
template < class Rep_, class NT_, class ET2NT_ >
void  QPE__filtered_base<Rep_,NT_,ET2NT_>::
init_NT( )
{
    // ToDo: scale 'x_B_O', 'lambda', and 'd' if necessary

    // get inexact version of 'lambda'
    std::transform( solver().lambda_numerator_begin(),
		    solver().lambda_numerator_end(),
		    lambda_NT.begin(), et2nt_obj);

    // get inexact version of 'x_B_O'
    init_NT( Is_linear());

    // get inexact version of 'd'
    d_NT = et2nt_obj( solver().variables_common_denominator());
}

template < class Rep_, class NT_, class ET2NT_ >
void  QPE__filtered_base<Rep_,NT_,ET2NT_>::
update_maxima( )
{
    // get properties of quadratic program
    int  n = solver().number_of_variables();
    int  c = solver().number_of_basic_constraints();
    int  b = solver().number_of_basic_original_variables();

    // initialize error bounds
    bound1    = d_NT * row_max_c;
    bound2_wq = d_NT;

    // update row and column maxima of 'A'
    A_iterator  a_it = solver().a_begin();
    R_iterator  r_it = solver().row_type_begin();
    int         row, col;
    NT          row_max, z;

    Basic_constraint_index_iterator  it;
    for ( it =  solver().basic_constraints_index_begin();
	  it != solver().basic_constraints_index_end(); ++it) {
	row = *it;

	// row not handled yet?
	if ( ! handled_A[ row]) {

	    // slack variable (or artificial in phase I) involved?
	    row_max = (    ( r_it[ row] != Rep::EQUAL)
			|| ( solver().phase() == 1) ? nt1 : nt0);

	    // scan row and update maxima
	    for ( col = 0; col < n; ++col) {
		z = CGAL::abs( a_it[ col][ row]);
		if ( z > row_max      ) row_max       = z;
		if ( z > col_max[ col]) col_max[ col] = z;
	    }
	    row_max_A[ row] = row_max;
	    handled_A[ row] = true;

	    // update bounds
	    z = CGAL::abs( lambda_NT[ row]);
	    if ( z > bound2_wq) bound2_wq = z;
	    z *= row_max_A[ row];
	    if ( z > bound1) bound1 = z;
	}
    }

    // update row and column maxima of 'D', if necessary
    update_maxima( Is_linear());

    // finalize error bounds
    // ToDo: use std::numeric_limits for 'machine epsilon'
    NT  q = std::ldexp( 1.015625 * ( c+b+1) * ( c+b+2), -53);
    bound1    *= q;
    bound2_wq *= q;

    CGAL_qpe_debug {
	vout() << std::endl
	       << "first bound for certification: " << bound1 << std::endl;
    }
}

template < class Rep_, class NT_, class ET2NT_ >                // QP case
void  QPE__filtered_base<Rep_,NT_,ET2NT_>::
update_maxima( Tag_false)
{
    int  n = solver().number_of_variables();

    // update row and column maxima of 'D'
    D_row_iterator  d_row_it;
    int             row, col;
    NT              row_max, z;

    Basic_variable_index_iterator  it;
    for ( it =  solver().basic_original_variables_index_begin();
	  it != solver().basic_original_variables_index_end(); ++it) {
	row = *it;

	// row not handled yet?
	if ( ! handled_D[ row]) {
	    d_row_it = solver().d_begin()[ row];

	    // scan row and update maxima
	    row_max = nt0;
	    for ( col = 0; col < n; ++col, ++d_row_it) {
		z = CGAL::abs( *d_row_it);
		if ( z > row_max      ) row_max       = z;
		if ( z > col_max[ col]) col_max[ col] = z;
	    }
	    row_max_D[ row] = row_max;
	    handled_D[ row] = true;

	    // update bounds
	    z = CGAL::abs( x_B_O_NT[ row]);
	    if ( z > bound2_wq) bound2_wq = z;
	    z *= row_max_D[ row];
	    if ( z > bound1) bound1 = z;
	}
    }
}

template < class Rep_, class NT_, class ET2NT_ >
bool  QPE__filtered_base<Rep_,NT_,ET2NT_>::
certify_mu_j_NT( int j) const
{
    // compute 'mu_j' with inexact arithmetic again
    NT  mu = mu_j_NT( j);

    CGAL_qpe_debug {
	vout() << "mu_" << j << " [NT]: " << mu;
    }

    // check against first bound
    if ( mu >= bound1) {
	CGAL_qpe_debug {
	    vout() << "  ok [1]" << std::endl;
	}
	return true;
    }

    // compute and check against second bound
    NT  bound2 = bound2_wq * col_max[ j];
    if ( mu >= bound2) {
	CGAL_qpe_debug {
	    vout() << "  ok [2: " << bound2 << "]" << std::endl;
	}
	return true;
    }

    // compute and check exact 'mu_j'
    ET  mu_et = mu_j( j);

    CGAL_qpe_debug {
	vout() << "  " << ( mu_et >= et0 ? "ok" : "MISSED")
	       << " [exact: " << mu_et << "]" << std::endl;
    }
    return ( mu_et >= et0);
}

// transition
template < class Rep_, class NT_, class ET2NT_ >
void  QPE__filtered_base<Rep_,NT_,ET2NT_>::
transition( )
{
    // get properties of quadratic program
    int  n = solver().number_of_variables();

    // update row and column maxima with original objective vector 'c'
    typename Rep::C_iterator  c_it = solver().c_begin();
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

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
