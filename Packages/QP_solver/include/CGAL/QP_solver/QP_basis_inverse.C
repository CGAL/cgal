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
// file          : include/CGAL/QP_engine/QPE_basis_inverse.C
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Basis Inverse for the Quadratic Programming Engine
// ============================================================================

CGAL_BEGIN_NAMESPACE

// =============================
// class implementation (cont'd)
// =============================

// creation and initialization
// ---------------------------
// set-up
template < class ET_, class Is_LP_ >
void  QPE_basis_inverse<ET_,Is_LP_>::
set( int n, int m, int nr_equalities)
{
    CGAL_qpe_precondition( n > 0);
    CGAL_qpe_precondition( m > 0);
    b = s = 0;
    // l is the maximum size of the basis in phase I
    l = std::min( n+nr_equalities+1, m);
    if ( ! M.empty()) M.clear();
    set( Is_LP());
    
    if ( ! x_l.empty()) x_l.clear();
    if ( ! x_x.empty()) x_x.clear();
   
    x_l.insert( x_l.end(), l, et0);
    x_x.insert( x_x.end(), nr_equalities+1, et0); // has to grow later
}

// update functions
// ----------------
// leaving of original variable (update type U2)
template < class ET_, class Is_LP_ >
void  QPE_basis_inverse<ET_,Is_LP_>::
leave_original( )
{
    // assert QP case
    Assert_compile_time_tag( Tag_false(), Is_LP());

    // determine new denominator (`z')
    --b;
    ET    z     = M[ l+b][ l+b];
    bool  z_neg = ( z < et0);
    CGAL_qpe_precondition( z != et0);

    // update matrix in place
    update_inplace_QP( M[ l+b].begin(), M[ l+b].begin()+l,
		       -z, ( z_neg ? d : -d));
                                                                 
    // store new denominator
    d = ( z_neg ? -z : z);
    CGAL_qpe_postcondition( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
}

// entering of slack variable (update type U3)
template < class ET_, class Is_LP_ >
void  QPE_basis_inverse<ET_,Is_LP_>::
enter_slack( )
{
    // assert QP case
    Assert_compile_time_tag( Tag_false(), Is_LP());

    // determine new denominator (`z')
    --s;
    ET    z     = M[ s][ s];
    bool  z_neg = ( z < et0);
    CGAL_qpe_precondition( z != et0);

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
    CGAL_qpe_postcondition( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
}

// replacing of original by slack variable (update type U8)
template < class ET_, class Is_LP_ >
void  QPE_basis_inverse<ET_,Is_LP_>::
enter_slack_leave_original( )
{
    // assert LP case or phase I
    CGAL_qpe_precondition( is_LP || is_phaseI);

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
    CGAL_qpe_precondition( z != et0);

    // update matrix
    update_inplace_LP( matrix_it->begin(), x_x.begin(), -z, ( z_neg ? d : -d));

    // store new denominator
    // ---------------------
    d = ( z_neg ? -z : z);
    CGAL_qpe_postcondition( d > et0);

    CGAL_qpe_debug {
        if ( vout.verbose()) print();
    }
}

// swap functions
// --------------
// swap variable ``to the end'' of R
template < class ET_, class Is_LP_ >                            // LP case
void  QPE_basis_inverse<ET_,Is_LP_>::
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
void  QPE_basis_inverse<ET_,Is_LP_>::
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
void  QPE_basis_inverse<ET_,Is_LP_>::
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
void  QPE_basis_inverse<ET_,Is_LP_>::
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
	matrix_it += l;

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
void  QPE_basis_inverse<ET_,Is_LP_>::
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

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
