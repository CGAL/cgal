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
// file          : include/CGAL/QPE_solver.C
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Quadratic Programming Engine - Solver
// ============================================================================

CGAL_BEGIN_NAMESPACE

// =============================
// class implementation (cont'd)
// =============================

// creation & initialization
// -------------------------
// creation
template < class Rep_ >
QPE_solver<Rep_>::
QPE_solver( )
    : et0( 0), et1( 1), et2( 2),
      strategyP( static_cast< Pricing_strategy*>( 0)),
      inv_M_B( vout4),
      d( inv_M_B.denominator()),
      m_phase( -1), is_phaseI( false), is_phaseII( false),
      is_LP( check_tag( Is_linear())), is_QP( ! is_LP),
      no_ineq( check_tag( Has_no_inequalities())), has_ineq( ! no_ineq)
{ }

// set-up of QP
template < class Rep_ >
void
QPE_solver<Rep_>::
set( int n, int m,
     typename QPE_solver<Rep_>::A_iterator A_it,
     typename QPE_solver<Rep_>::B_iterator b_it,
     typename QPE_solver<Rep_>::C_iterator c_it,
     typename QPE_solver<Rep_>::D_iterator D_it,
     typename QPE_solver<Rep_>::Row_type_iterator Row_type_it)
{
    // store QP
    CGAL_qpe_precondition( n > 0);
    CGAL_qpe_precondition( m > 0);
    qp_n = n; qp_m = m;
    qp_A = A_it; qp_b = b_it; qp_c = c_it; qp_D = D_it;
    qp_r = Row_type_it;
    l = std::min( n, m);

    // set up slack variables and auxiliary problem
    // --------------------------------------------
    // clear slack and artificial part of `A', if necessary
    if ( ! slack_A.empty()) slack_A.clear();
    if ( !   art_A.empty())   art_A.clear();
    if ( !   art_s.empty())   art_s.clear();

    // reserve memory for slack and artificial part of `A'
    if ( has_ineq) {
	unsigned int  art_size = std::count( qp_r, qp_r+qp_m, Rep::EQUAL);
	slack_A.reserve( qp_m - art_size);
	  art_A.reserve(        art_size);
          art_s.insert ( art_s.end(), qp_m, A_entry( 0));
    } else {
	  art_A.reserve( qp_m);
    }

    // initialize slack and artificial part of `A'
    B_entry  b0( 0), b_max( b0);
    C_entry  c1( 1);
    int              i_max = -1;
    for ( int i = 0; i < qp_m; ++i) {

        if ( has_ineq && ( qp_r[ i] != Rep::EQUAL)) {   // slack variable

	    if ( qp_r[ i] == Rep::LESS_EQUAL) {                 // '<='
		if ( qp_b[ i] < b0) {

                    // special entry '< -0'
		    art_s[ i] = -c1;
		    if ( -qp_b[ i] > b_max) {
			i_max = slack_A.size();
			b_max = -qp_b[ i];
		    }
		}

		// store slack column
		slack_A.push_back( std::make_pair( i, false));

	    } else {                                            // '>='
		if ( qp_b[ i] > b0) {

                    // special entry '> +0'
		    art_s[ i] = c1;
		    if (  qp_b[ i] > b_max) {
			i_max = slack_A.size();
			b_max = qp_b[ i];
		    }
		}

		// store slack column
		slack_A.push_back( std::make_pair( i, true));
	    }

        } else {                                        // artificial variable

            art_A.push_back( std::make_pair( i, qp_b[ i] < b0));

	}
    }

    // special artificial column needed?
    if ( i_max >= 0) {
	art_s_i = -i_max;
    } else {
	art_s_i = -1;
	art_s.clear();
    }

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout.verbose()) {
	    if ( vout2.verbose()) {
		vout2.out() << "======" << std::endl
			    << "Set-Up" << std::endl
			    << "======" << std::endl;
	    }
	    vout.out() << "[ " << ( is_LP ? "LP" : "QP")
		       << ", " << n << " variables, " << m << " constraints";
	    if ( vout2.verbose() && ( ! slack_A.empty())) {
		vout2.out() << " (" << slack_A.size() << " inequalities)";
	    }
	    vout.out() << " ]" << std::endl;
	    if ( vout2.verbose()) {
		if ( is_QP) {
		    vout2.out() << "flag: D "
				<< ( check_tag( Is_symmetric()) ? "" : "not ")
				<< "symmetric" << std::endl;
		}
		vout2.out() << "flag: " << ( has_ineq ? "has" : "no")
			    << " inequality constraints" << std::endl;
		if ( vout4.verbose()) print_program();
	    }
	}
    }

    // set up pricing strategy
    if ( strategyP != static_cast< Pricing_strategy*>( 0)) {
	strategyP->set( *this, vout2);
    }

    // set up basis inverse
    inv_M_B.set( qp_n, qp_m);

    // set phase
    m_phase    = 0;
    is_phaseI  = false;
    is_phaseII = false;
}

// initialization (phase I)
template < class Rep_ >
void
QPE_solver<Rep_>::
init( )
{
    CGAL_qpe_debug {
        vout2 << std::endl
              << "==============" << std::endl
              << "Initialization" << std::endl
              << "==============" << std::endl;
    }

    // set status
    m_phase    = 1;
    m_status   = UPDATE;
    m_pivots   = 0;
    is_phaseI  = true;
    is_phaseII = false;

    // initial basis and basis inverse
    init_basis();

    // initial solution
    init_solution();

    // initialize additional data members
    init_additional_data_members();

    // initialize pricing strategy
    CGAL_qpe_precondition( strategyP != static_cast< Pricing_strategy*>( 0));
    strategyP->init( 0);

    // basic feasible solution already available?
    if ( art_basic == 0) {

	// transition to phase II
	CGAL_qpe_debug {
	    if ( vout2.verbose()) {
		vout2.out() << std::endl
			    << "no artificial variables at all "
			    << "--> skip phase I"
			    << std::endl;
	    }
	}
	transition();
    }
}

// initial basis and basis inverse
template < class Rep_ >
void
QPE_solver<Rep_>::
init_basis( )
{
    int  i, s_i = -1;
    int  j, s   = slack_A.size();

    // has special artificial column?
    if ( ! art_s.empty()) {
	    s_i = -art_s_i;
	art_s_i = qp_n+s+art_A.size();
	art_A.push_back( std::make_pair( s_i, ! slack_A[ s_i].second));
    }

    // initialize indices of basic variables
    if ( ! in_B.empty()) in_B.clear();
    in_B.reserve( qp_n+s+art_A.size());
    in_B.insert( in_B.end(), qp_n, -1);        // no original variable is basic

    init_basis__slack_variables( s_i, Has_no_inequalities());

    if ( ! B_O.empty()) B_O.clear();
    B_O.reserve( qp_n);                         // all artificial variables are
    for ( i = 0, j = qp_n+s; i < (int)art_A.size(); ++i, ++j) {        // basic
	   B_O.push_back( j);
	in_B  .push_back( i);
    }
    art_basic = art_A.size();

    // initialize indices of 'basic' and 'nonbasic' constraints
    if ( ! C.empty()) C.clear();
    init_basis__constraints( s_i, Has_no_inequalities());

    // diagnostic output
    CGAL_qpe_debug {
        if ( vout.verbose()) print_basis();
    }

    // initialize basis inverse (explain: 'art_s' not needed here)
    inv_M_B.init( art_A.size(), art_A.begin());
}

template < class Rep_ >                                        // has ineq.
void  QPE_solver<Rep_>::
init_basis__slack_variables( int s_i, Tag_false)
{
    int  s = slack_A.size();

    // reserve memory
    if ( ! B_S.empty()) B_S.clear();
    B_S.reserve( s);

    // (almost) all slack variables are basic
    // Note: slack var. corresponding to special artificial var. is nonbasic
    for ( i = 0, j = qp_n; i < s; ++i, ++j) {
	if ( i != s_i) {
	    in_B  .push_back( B_S.size());
	       B_S.push_back( j);     
	} else {
	    in_B  .push_back( -1);
	}
    }
}

template < class Rep_ >                                        // has ineq.
void  QPE_solver<Rep_>::
init_basis__constraints( int s_i, Tag_false)
{
    // reserve memory
    if ( ! in_C  .empty()) in_C  .clear();
    if ( !    S_B.empty())    S_B.clear();
      C.reserve( l);
    S_B.reserve( slack_A.size());

    // store constraints' indices
    in_C.insert( in_C.end(), qp_m, -1);
    if ( s_i >= 0) s_i = slack_A[ s_i].first;
    for ( i = 0, j = 0; i < qp_m; ++i) {
	if ( qp_r[ i] == Rep::EQUAL) {          // equal. constraints are basic
	       C.push_back( i);
	    in_C[ i] = j;
	    ++j;
	} else {                                // ineq. constrai. are nonbasic
	    if ( i != s_i) S_B.push_back( i);
	}
    }
    if ( s_i >= 0) {                            // special ineq. con. is basic
	   C.push_back( s_i);
	in_C[ s_i] = j;
    }
}

// initial solution
template < class Rep_ >
void
QPE_solver<Rep_>::
init_solution( )
{
    // initialize exact version of `qp_b' restricted to basic constraints C
    // (implicit conversion to ET)
    if ( ! b_C.empty()) b_C.clear();
    init_solution__b_C( Has_no_inequalities());

    // initialize exact version of `-aux_c' restricted to basic variables B_O
    if ( ! minus_c_B.empty()) minus_c_B.clear();
    minus_c_B.insert( minus_c_B.end(), l, -et1);
    if ( art_s_i > 0) minus_c_B[ art_A.size()-1] *= ET( qp_n+qp_m);

    // allocate memory for current solution
    if ( ! lambda.empty()) lambda.clear();
    if ( ! x_B_O .empty()) x_B_O .clear();
    if ( ! x_B_S .empty()) x_B_S .clear();
    lambda.insert( lambda.end(), l, et0);
    x_B_O .insert( x_B_O .end(), l, et0);
    x_B_S .insert( x_B_S .end(), slack_A.size(), et0);

    // compute initial solution
    compute_solution();

    // diagnostic output
    CGAL_qpe_debug {
        if ( vout.verbose()) print_solution();
    }
}

// initialize additional data members
template < class Rep_ >
void
QPE_solver<Rep_>::
init_additional_data_members( )
{
    if ( ! A_Cj.empty()) A_Cj.clear();
    A_Cj.insert( A_Cj.end(), l, et0);
    if ( ! two_D_Bj.empty()) two_D_Bj.clear();
    two_D_Bj.insert( two_D_Bj.end(), l, et0);

    if ( ! q_lambda.empty()) q_lambda.clear();
    q_lambda.insert( q_lambda.end(), l, et0);
    if ( ! q_x_O.empty()) q_x_O.clear();
    q_x_O.insert( q_x_O.end(), l, et0);
    if ( ! q_x_S.empty()) q_x_S.clear();
    q_x_S.insert( q_x_S.end(), slack_A.size(), et0);
    
    if ( ! tmp_l.empty()) tmp_l.clear();
    tmp_l.insert( tmp_l.end(), l, et0);
    if ( ! tmp_l_2.empty()) tmp_l_2.clear();
    tmp_l_2.insert( tmp_l_2.end(), l, et0);
    if ( ! tmp_x.empty()) tmp_x.clear();
    tmp_x.insert( tmp_x.end(), l, et0);
    if ( ! tmp_x_2.empty()) tmp_x_2.clear();
    tmp_x_2.insert( tmp_x_2.end(), l, et0);
}

// transition (to phase II)
// ------------------------
template < class Rep >
void
QPE_solver<Rep>::
transition( )
{
    CGAL_qpe_debug {
	if ( vout.verbose()) {
	    vout2 << std::endl
		  << "----------------------" << std::endl
		  << 'T';
	    vout1 << "[ t"; vout  << "ransition to phase II"; vout1 << " ]";
	    vout  << std::endl;
	    vout2 << "----------------------";
	}
    }

    // update status
    m_phase    = 2;
    is_phaseI  = false;
    is_phaseII = true;

    // remove artificial variables
    in_B.erase( in_B.begin()+qp_n+slack_A.size(), in_B.end());

    // update basis inverse
    CGAL_qpe_debug {
	vout4 << std::endl << "basis-inverse:" << std::endl;
    }
    transition( Is_linear());
    check_basis_inverse();

    // initialize exact version of `-qp_c' (implicit conversion to ET)
    C_by_index_accessor  c_accessor( qp_c);
    std::transform( C_by_index_iterator( B_O.begin(), c_accessor),
                    C_by_index_iterator( B_O.end  (), c_accessor),
                    minus_c_B.begin(), std::negate<ET>());

    // compute initial solution of phase II
    compute_solution();

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout.verbose()) print_solution();
    }

    // notify pricing strategy
    strategyP->transition();
}

// access
// ------
// numerator of current solution
template < class Rep_ >
typename QPE_solver<Rep_>::ET
QPE_solver<Rep_>::
solution_numerator( ) const
{
    Index_const_iterator  i_it;
    Value_const_iterator  x_i_it, c_it;
    ET   s, z = et0;
    int  i;

    // foreach i
    x_i_it =       x_B_O.begin();
      c_it = minus_c_B  .begin();
    for ( i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++x_i_it, ++c_it){
        i = *i_it;

        // quadratic part
        s = et0;
        if ( is_QP && is_phaseII) {
       
            // foreach j < i
	    s += std::inner_product(x_B_O.begin(), x_i_it,
				D_pairwise_iterator(
				    B_O.begin(),
				    D_pairwise_accessor( qp_D, i)),
				et0);

            // D_{i,i} x_i
            s += ET( qp_D[ i][ i]) * *x_i_it;
        }
        // linear part
        s -= d * *c_it;

        // accumulate
        z += s * *x_i_it;
    }
    return z;
}

// pivot step
// ----------
template < class Rep_ >
void
QPE_solver<Rep_>::
pivot_step( )
{
    ++m_pivots;

    // diagnostic output
    CGAL_qpe_debug {
        vout2 << std::endl
              << "==========" << std::endl
              << "Pivot Step" << std::endl
              << "==========" << std::endl;
        vout  << "[ phase " << ( is_phaseI ? "I" : "II")
              << ", iteration " << m_pivots << " ]" << std::endl;
    }

    // pricing
    // -------
    pricing();

    // check for optimality
    if ( j < 0) {

        if ( is_phaseI) {                               // phase I

            // problem is infeasible
	    m_phase  = 3;
	    m_status = INFEASIBLE;
	    CGAL_qpe_debug {
		vout1 << "  ";
		vout << "problem is INFEASIBLE" << std::endl;
	    }

        } else {                                        // phase II

	    // optimal solution found
	    m_phase  = 3;
            m_status = OPTIMAL;
            CGAL_qpe_debug {
                vout1 << "  ";
		vout2 << std::endl;
                vout  << "solution is OPTIMAL" << std::endl;
            }
        }
        return;
    }

    // ratio test & update (step 1)
    // ----------------------------
    // initialize ratio test
    ratio_test_init();

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose() && is_QP && is_phaseII) {
	    vout2.out() << std::endl
			<< "----------------------------" << std::endl
			<< "Ratio Test & Update (Step 1)" << std::endl
			<< "----------------------------" << std::endl;
	}
    }

    // loop (step 1)
    do {

        // ratio test
        ratio_test_1();

        // check for unboundedness
        if ( q_i == et0) {
            m_phase  = 3;
            m_status = UNBOUNDED;
            CGAL_qpe_debug {
                vout1 << "  ";
                vout << "problem is UNBOUNDED" << std::endl;
            }
            return;
        }

        // update
        update_1();

    } while ( j >= 0);

    // ratio test & update (step 2)
    // ----------------------------
    if ( i >= 0) {

	// diagnostic output
	CGAL_qpe_debug {
	    vout2 << std::endl
		  << "----------------------------" << std::endl
		  << "Ratio Test & Update (Step 2)" << std::endl
		  << "----------------------------" << std::endl;
	}

	// compute index of entering variable
	j += in_B.size();

	// loop (step 2)
	while ( ( i >= 0) && basis_matrix_stays_regular()) {

	    // update
	    update_2( Is_linear());

	    // ratio test
	    ratio_test_2( Is_linear());
	}
    }

    // ratio test & update (step 3)
    // ----------------------------
    CGAL_qpe_assertion_msg( i < 0, "Step 3 should never be reached!");

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout1.verbose()) print_basis();
	if ( vout .verbose()) print_solution();
    }

    // transition to phase II (if possible)
    // ------------------------------------
    if ( is_phaseI && ( art_basic == 0)) {
	CGAL_qpe_debug {
	    if ( vout2.verbose()) {
		vout2.out() << std::endl
			    << "all artificial variables are nonbasic"
			    << std::endl;
	    }
	}
	transition();
    }
}

// pricing
template < class Rep_ >
void
QPE_solver<Rep_>::
pricing( )
{
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    vout2 << std::endl
		  << "-------" << std::endl
		  << "Pricing" << std::endl
		  << "-------" << std::endl;
	}
    }

    // call pricing strategy
    j = strategyP->pricing();

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout.verbose()) {
	    if ( j < 0) {
		vout2 << "entering variable: none" << std::endl;
	    } else {
		vout1 << "  ";
		vout  << "entering"; vout2 << " variable"; vout << ": ";
		vout  << j;
		vout2 << " (" << variable_type( j) << ')' << std::endl;
	    }
	}
    }
}

// initialization of ratio-test
template < class Rep_ >
void
QPE_solver<Rep_>::
ratio_test_init( )
{
    // store exact version of `A_Cj' (implicit conversion)
    ratio_test_init__A_Cj( A_Cj.begin(), j, Has_no_inequalities());

    // store exact version of `2 D_{B_O,j}'
    ratio_test_init__2_D_Bj( two_D_Bj.begin(), j, Is_linear());
}

template < class Rep_ >                                         // no ineq.
void  QPE_solver<Rep_>::
ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j_, Tag_true)
{
    // store exact version of `A_Cj' (implicit conversion)
    if ( j_ < qp_n) {                                   // original variable

	std::copy( qp_A[ j_], qp_A[ j_]+qp_m, A_Cj_it);

    } else {                                            // artificial variable

	unsigned int  k = j_;
	k -= qp_n;
	std::fill_n( A_Cj_it, qp_m, et0);
	A_Cj_it[ k] = ( art_A[ k].second ? -et1 : et1);
    }
}

template < class Rep_ >                                        // has ineq.
void  QPE_solver<Rep_>::
ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j_, Tag_false)
{
    // store exact version of `A_Cj' (implicit conversion)
    if ( j_ < qp_n) {                                   // original variable

	A_by_index_accessor  a_accessor( qp_A[ j_]);
	std::copy( A_by_index_iterator( C.begin(), a_accessor),
		   A_by_index_iterator( C.end  (), a_accessor),
		   A_Cj_it);

    } else {
	unsigned int  k = j_;
	k -= qp_n;
	std::fill_n( A_Cj_it, C.size(), et0);

	if ( k < slack_A.size()) {                      // slack variable

	    A_Cj_it[ in_C[ slack_A[ k].first]] = ( slack_A[ k].second ? -et1
						                      :  et1);

	} else {                                        // artificial variable
	    k -= slack_A.size();

	    if ( j_ != art_s_i) {                           // normal art.

		A_Cj_it[ in_C[ art_A[ k].first]] = ( art_A[ k].second ? -et1
						                      :  et1);

	    } else {                                        // special art.

		S_by_index_accessor  s_accessor( art_s.begin());
		std::copy( S_by_index_iterator( C.begin(), s_accessor),
			   S_by_index_iterator( C.end  (), s_accessor),
			   A_Cj_it);
	    }	
	}
    }
}

// ratio test (step 1)
template < class Rep_ >
void
QPE_solver<Rep_>::
ratio_test_1( )
{
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    vout2.out() << std::endl;
	    if ( is_LP || is_phaseI) {
		vout2.out() << "----------" << std::endl
			    << "Ratio Test" << std::endl
			    << "----------" << std::endl;
	    } else {
		vout2.out() << "Ratio Test (Step 1)" << std::endl
			    << "-------------------" << std::endl;
	    }
	    if ( vout3.verbose()) {
		vout3.out() << "    A_Cj: ";
		std::copy( A_Cj.begin(), A_Cj.begin()+C.size(),
			   std::ostream_iterator<ET>( vout3.out()," "));
		vout3.out() << std::endl;
		if ( is_QP && is_phaseII) {
		    vout3.out() << "  2 D_Bj: ";
		    std::copy( two_D_Bj.begin(), two_D_Bj.begin()+B_O.size(),
			       std::ostream_iterator<ET>( vout3.out()," "));
		    vout3.out() << std::endl;
		}
		vout3.out() << std::endl;
	    }
	}
    }

    // compute `q_lambda' and `q_x'
    ratio_test_1__q_x_O( Is_linear());
    ratio_test_1__q_x_S( Has_no_inequalities());

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout3.verbose()) {
	    if ( is_QP && is_phaseII) {
		vout3.out() << "q_lambda: ";
		std::copy( q_lambda.begin(), q_lambda.begin()+C.size(),
			   std::ostream_iterator<ET>( vout3.out()," "));
		vout3.out() << std::endl;
	    }
	    vout3.out() << "   q_x_O: ";
	    std::copy( q_x_O.begin(), q_x_O.begin()+B_O.size(),
		       std::ostream_iterator<ET>( vout3.out()," "));
	    vout3.out() << std::endl;

	    if ( has_ineq) {
		vout3.out() << "   q_x_S: ";
		std::copy( q_x_S.begin(), q_x_S.begin()+B_S.size(),
			   std::ostream_iterator<ET>( vout3.out()," "));
		vout3.out() << std::endl;
	    }
	    vout3.out() << std::endl;
	}
    }

    // check `t_i's
    x_i = et1;                                          // trick: initialize
    q_i = et0;                                          // minimum with +oo

    // what happens, if all original variables are nonbasic?
    ratio_test_1__t_i(   B_O.begin(),   B_O.end(),
		       x_B_O.begin(), q_x_O.begin(), Tag_false());
    ratio_test_1__t_i(   B_S.begin(),   B_S.end(),
		       x_B_S.begin(), q_x_S.begin(), Has_no_inequalities());

    // check `t_j'
    ratio_test_1__t_j( Is_linear());

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    for ( unsigned int k = 0; k < B_O.size(); ++k) {
		vout2.out() << "t_O_" << k << ": "
			    << x_B_O[ k] << '/' << q_x_O[ k]
			    << ( ( q_i > et0) && ( i == B_O[ k]) ? " *" : "")
			    << std::endl;
	    }
	    if ( has_ineq) {
		for ( unsigned int k = 0; k < B_S.size(); ++k) {
		    vout2.out() << "t_S_" << k << ": "
				<< x_B_S[ k] << '/' << q_x_S[ k]
				<< ( ( q_i > et0) && ( i == B_S[ k]) ? " *":"")
				<< std::endl;
		}
	    }
	    if ( is_QP && is_phaseII) {
		vout2.out() << std::endl
			    << "  t_j: " << mu << '/' << nu
			    << ( ( q_i > et0) && ( i < 0) ? " *" : "")
			    << std::endl;
	    }
	    vout2.out() << std::endl;
	}
	if ( q_i > et0) {
	    if ( i < 0) {
		vout2 << "leaving variable: none" << std::endl;
	    } else {
		vout1 << ", ";
		vout  << "leaving"; vout2 << " variable"; vout << ": ";
		vout  << i;
		if ( vout2.verbose()) {
		    if ( ( i < qp_n) || ( i >= (int)( qp_n+slack_A.size()))) {
			vout2.out() << " (= B_O[ " << in_B[ i] << "]: "
				    << variable_type( i) << ')' << std::endl;
		    } else {
			vout2.out() << " (= B_S[ " << in_B[ i] << "]: slack)"
				    << std::endl;
		    }
		}
	    }
	}
    }
}

template < class Rep_ >                                         // QP case
void
QPE_solver<Rep_>::
ratio_test_2( Tag_false)
{
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    vout2.out() << std::endl
			<< "Ratio Test (Step 2)" << std::endl
			<< "-------------------" << std::endl;
	}
    }

    // compute `p_lambda' and `p_x' (Note: `p_...' is stored in `q_...')
    ratio_test_2__p( Has_no_inequalities());
 
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout3.verbose()) {
	    vout3.out() << "p_lambda: ";
	    std::copy( q_lambda.begin(), q_lambda.begin()+C.size(),
		       std::ostream_iterator<ET>( vout3.out()," "));
	    vout3.out() << std::endl;
	    vout3.out() << "   p_x_O: ";
	    std::copy( q_x_O.begin(), q_x_O.begin()+B_O.size(),
		       std::ostream_iterator<ET>( vout3.out()," "));
	    vout3.out() << std::endl;
	    if ( has_ineq) {
		vout3.out() << "   p_x_S: ";
		std::copy( q_x_S.begin(), q_x_S.begin()+B_S.size(),
			   std::ostream_iterator<ET>( vout3.out()," "));
		vout3.out() << std::endl;
	    }
	    vout3.out() << std::endl;
	}
    }

    // check `mu_j_i's
    x_i = et0;                                          // initialize
    q_i = et1;                                          // minimum with 0

    Value_iterator  x_it = x_B_O.begin();
    Value_iterator  q_it = q_x_O.begin();
    Index_iterator  i_it;
    for ( i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++x_it, ++q_it) {
	if ( ( *q_it > et0) && ( ( *x_it * q_i) > ( x_i * *q_it))) {
	    i = *i_it; x_i = *x_it; q_i = *q_it;
	}
    }
    x_it = x_B_S.begin();
    q_it = q_x_S.begin();
    for ( i_it = B_S.begin(); i_it != B_S.end(); ++i_it, ++x_it, ++q_it) {
	if ( ( *q_it > et0) && ( ( *x_it * q_i) > ( x_i * *q_it))) {
	    i = *i_it; x_i = *x_it; q_i = *q_it;
	}
    }

    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    for ( unsigned int k = 0; k < B_O.size(); ++k) {
		vout2.out() << "mu_j_O_" << k << ": - "
			    << x_B_O[ k] << '/' << q_x_O[ k]
			    << ( ( q_i < et0) && ( i == B_O[ k]) ? " *" : "")
			    << std::endl;
	    }
	    for ( unsigned int k = 0; k < B_S.size(); ++k) {
		vout2.out() << "mu_j_S_" << k << ": - "
			    << x_B_S[ k] << '/' << q_x_S[ k]
			    << ( ( q_i < et0) && ( i == B_S[ k]) ? " *" : "")
			    << std::endl;
	    }
	    vout2.out() << std::endl;
	}
	if ( i < 0) {
	    vout2 << "leaving variable: none" << std::endl;
	} else {
	    vout1 << ", ";
	    vout  << "leaving"; vout2 << " variable"; vout << ": ";
	    vout  << i;
	    if ( vout2.verbose()) {
		if ( i < qp_n) {
		    vout2.out() << " (= B_O[ " << in_B[ i] << "]: original)"
				<< std::endl;
		} else {
		    vout2.out() << " (= B_S[ " << in_B[ i] << "]: slack)"
				<< std::endl;
		}
	    }
	}
    }
}

// update (step 1)
template < class Rep_ >
void
QPE_solver<Rep_>::
update_1( )
{
    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    vout2.out() << std::endl;
	    if ( is_LP || is_phaseI) {
		vout2.out() << "------" << std::endl
			    << "Update" << std::endl
			    << "------" << std::endl;
	    } else {
		vout2.out() << "Update (Step 1)" << std::endl
			    << "---------------" << std::endl;
	    }
	}
    }

    // update basis & basis inverse
    update_1( Is_linear());
    check_basis_inverse();

    // compute current solution
    compute_solution();
}

// update (step 2)
template < class Rep_ >                                         // QP case
void
QPE_solver<Rep_>::
update_2( Tag_false)
{
    CGAL_qpe_debug {
	vout2 << std::endl
	      << "Update (Step 2)" << std::endl
	      << "---------------" << std::endl;
    }

    // leave variable from basis
    leave_variable();
    check_basis_inverse();

    // compute current solution
    compute_solution();
}

// replace variable in basis
template < class Rep_ >
void
QPE_solver<Rep_>::
replace_variable( )
{
    CGAL_qpe_debug {
	vout2 <<   "<--> nonbasic (" << variable_type( j) << ") variable " << j
	      << " replaces basic (" << variable_type( i) << ") variable " << i
	      << std::endl << std::endl;
    }

    // replace variable
    replace_variable( Has_no_inequalities());

    // pivot step done
    i = j = -1;
}

template < class Rep_ >
void  QPE_solver<Rep_>::
replace_variable_original_original( )
{
    int  k = in_B[ i];

    // replace original variable [ in: j | out: i ]
    in_B  [ i] = -1;
    in_B  [ j] = k;
       B_O[ k] = j;

    minus_c_B[ k] = ( is_phaseI ? ( j < qp_n ? et0 : -et1) : -ET( qp_c[ j]));

    if ( is_phaseI) {
	if ( j >= qp_n) ++art_basic;
	if ( i >= qp_n) --art_basic;
    }

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }
	    
    // update basis inverse
    inv_M_B.enter_original_leave_original( q_x_O.begin(), k);
}

template < class Rep_ >
void  QPE_solver<Rep_>::
replace_variable_slack_slack( )
{
    int  k = in_B[ i];

    // replace slack variable [ in: j | out: i ]
    in_B  [ i] = -1;
    in_B  [ j] = k;
       B_S[ k] = j;
       S_B[ k] = slack_A[ j-qp_n].first;

    // replace inequality constraint [ in: i | out: j ]
    int old_row = S_B[ k];
    int new_row = slack_A[ i-qp_n].first;
    k = in_C[ old_row];

    in_C[ old_row] = -1;
    in_C[ new_row] = k;
       C[ k      ] = new_row;

     b_C[ k] = ET( qp_b[ new_row]);

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }

    // update basis inverse
    A_row_by_index_accessor  a_accessor( A_accessor( qp_A, 0, qp_n), new_row);
    std::copy( A_row_by_index_iterator( B_O.begin(), a_accessor),
	       A_row_by_index_iterator( B_O.end  (), a_accessor),
	       tmp_x.begin());
    if ( art_s_i > 0) {                                 // special artificial
	tmp_x[ in_B[ art_s_i]] = ET( art_s[ new_row]);
    }
    inv_M_B.enter_slack_leave_slack( tmp_x.begin(), k);
}

template < class Rep_ >
void  QPE_solver<Rep_>::
replace_variable_slack_original( )
{
    int  k = in_B[ i];

    // leave original variable [ out: i ]
    in_B  [ i         ] = -1; 
    in_B  [ B_O.back()] = k;
       B_O[ k] = B_O.back();
       B_O.pop_back();

    minus_c_B[ k] = minus_c_B[ B_O.size()];

    if ( is_phaseI && ( i >= qp_n)) --art_basic;

    // enter slack variable [ in: j ]
    int  old_row = slack_A[ j-qp_n].first;
    in_B  [ j] = B_S.size();
       B_S.push_back( j);
       S_B.push_back( old_row);

    // leave inequality constraint [ out: j ]
    int  l = in_C[ old_row];
     b_C[ l       ] = b_C[ C.size()-1];
       C[ l       ] = C.back();
    in_C[ C.back()] = l;
    in_C[ old_row ] = -1;
       C.pop_back();

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }

    // update basis inverse
    inv_M_B.swap_variable( k);
    inv_M_B.swap_constraint( l);
    inv_M_B.enter_slack_leave_original();
}

template < class Rep_ >
void  QPE_solver<Rep_>::
replace_variable_original_slack( )
{
    int  k = in_B[ i];

    // enter original variable [ in: j ]
    minus_c_B[ B_O.size()]
	= ( is_phaseI ? ( j < qp_n ? et0 : -et1) : -ET( qp_c[ j]));
    

    in_B  [ j] = B_O.size();
       B_O.push_back( j);

    if ( is_phaseI && ( j >= qp_n)) ++art_basic;

    // leave slack variable [ out: i ]
       B_S[ k         ] = B_S.back();
       S_B[ k         ] = S_B.back();
    in_B  [ B_S.back()] = k;
    in_B  [ i         ] = -1; 
       B_S.pop_back();
       S_B.pop_back();

    // enter inequality constraint [ in: i ]
    int new_row = slack_A[ i-qp_n].first;
     b_C[ C.size()] = ET( qp_b[ new_row]);
    in_C[ new_row ] = C.size();
       C.push_back( new_row);

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }

    // update basis inverse
    A_row_by_index_accessor  a_accessor( A_accessor( qp_A, 0, qp_n), new_row);
    std::copy( A_row_by_index_iterator( B_O.begin(), a_accessor),
	       A_row_by_index_iterator( B_O.end  (), a_accessor),
	       tmp_x.begin());
    if ( art_s_i > 0) {                                 // special art.
	tmp_x[ in_B[ art_s_i]] = ET( art_s[ new_row]);
    }
    inv_M_B.enter_original_leave_slack( q_x_O.begin(), tmp_x.begin());
    
}

// enter variable into basis
template < class Rep_ >
void
QPE_solver<Rep_>::
enter_variable( )
{
    CGAL_qpe_debug {
	vout2 << "--> nonbasic (" << variable_type( j) << ") variable "
	      << j << " enters basis" << std::endl << std::endl;
    }

    // update basis & basis inverse
    if ( no_ineq || ( j < qp_n)) {                      // original variable

	// enter original variable [ in: j ]
	if ( minus_c_B.size() == B_O.size()) {
	    minus_c_B.push_back( et0);
	        q_x_O.push_back( et0);
	      tmp_x  .push_back( et0);
	}
	minus_c_B[ B_O.size()] = -ET( qp_c[ j]);

	in_B  [ j] = B_O.size();
	   B_O.push_back( j);

	// diagnostic output
	CGAL_qpe_debug {
	    if ( vout2.verbose()) print_basis();
	}
	    
	// update basis inverse
	inv_M_B.enter_original( q_lambda.begin(), q_x_O.begin(), -nu);
	
    } else {                                            // slack variable

	// enter slack variable [ in: j ]
	in_B  [ j] = B_S.size();
	   B_S.push_back( j);
	   S_B.push_back( slack_A[ j-qp_n].first);

	// leave inequality constraint [ out: j ]
	int old_row = slack_A[ j-qp_n].first;
	int k = in_C[old_row];
	
	// reflect change of active constraints heading C in b_C
	b_C[ k] = b_C[C.size()-1];
		
	   C[ k] = C.back();
	in_C[ C.back()      ] = k;
	in_C[ old_row       ] = -1;
	   C.pop_back();
	
	// diagnostic output
	CGAL_qpe_debug {
	    if ( vout2.verbose()) print_basis();
	}

	// update basis inverse
	inv_M_B.swap_constraint( k);
	inv_M_B.enter_slack();
    }

    // variable entered
    j -= in_B.size();
}

// leave variable from basis
template < class Rep_ >
void
QPE_solver<Rep_>::
leave_variable( )
{
    CGAL_qpe_debug {
	vout2 << "<-- basic (" << variable_type( i) << ") variable "
	      << i << " leaves basis" << std::endl << std::endl;
    }

    // update basis & basis inverse
    int  k = in_B[ i];
    if ( no_ineq || ( i < qp_n)) {                      // original variable

	// leave original variable [ out: i ]
	in_B  [ i         ] = -1; 
	in_B  [ B_O.back()] = k;
	   B_O[ k] = B_O.back(); B_O.pop_back();

	minus_c_B [ k] = minus_c_B [ B_O.size()];
	  two_D_Bj[ k] =   two_D_Bj[ B_O.size()];

	// diagnostic output
	CGAL_qpe_debug {
	    if ( vout2.verbose()) print_basis();
	}

	// update basis inverse
	inv_M_B.swap_variable( k);
	inv_M_B.leave_original();

    } else {                                            // slack variable

	// leave slack variable [ out: i ]
	in_B  [ i         ] = -1; 
	in_B  [ B_S.back()] = k;
	   B_S[ k] = B_S.back(); B_S.pop_back();
	   S_B[ k] = S_B.back(); S_B.pop_back();

	// enter inequality constraint [ in: i ]
	int new_row = slack_A[ i-qp_n].first;

	A_Cj[ C.size()] = ( j < qp_n ? ET( qp_A[ j][ new_row]) : et0);

	 b_C[ C.size()] = ET( qp_b[ new_row]);
	in_C[ new_row ] = C.size();
	   C.push_back( new_row);

	// diagnostic output
	CGAL_qpe_debug {
	    if ( vout2.verbose()) print_basis();
	}

	// update basis inverse
	A_row_by_index_accessor  a_accessor( A_accessor( qp_A, 0, qp_n),
					     new_row);
	std::copy( A_row_by_index_iterator( B_O.begin(), a_accessor),
		   A_row_by_index_iterator( B_O.end  (), a_accessor),
		   tmp_x.begin());
	inv_M_B.leave_slack( tmp_x.begin());
    }

    // notify pricing strategy
    strategyP->leaving_basis( i);

    // variable left
    i = -1;
}

// compute solution
template < class Rep_ >
void
QPE_solver<Rep_>::
compute_solution()
{
    // compute current solution
    inv_M_B.multiply( b_C.begin(), minus_c_B.begin(),
                      lambda.begin(), x_B_O.begin());
    compute__x_B_S( Has_no_inequalities());
}

template < class Rep_ >
void  QPE_solver<Rep_>::
multiply__A_S_BxB_O( Value_iterator in, Value_iterator out) const
{
    // initialize
    std::fill_n( out, B_S.size(), et0);

    // foreach original column of A in B_O (artificial columns are zero in S_B
    A_column              a_col;                             // except special)
    Index_const_iterator  row_it, col_it;
    Value_iterator        out_it;
    ET                    in_value;
    for ( col_it = B_O.begin(); col_it != B_O.end(); ++col_it, ++in) {
	in_value = *in;
	out_it   = out;

	if ( *col_it < qp_n) {	                        // original variable
	    a_col = qp_A[ *col_it];

	    // foreach row of A in S_B
	    for ( row_it = S_B.begin(); row_it != S_B.end(); ++row_it,
		                                             ++out_it) {
		*out_it += ET( a_col[ *row_it]) * in_value;
	    }
	} else {
	    if ( *col_it == art_s_i) {                  // special artificial

		// foreach row of 'art_s'
		for ( row_it = S_B.begin(); row_it != S_B.end(); ++row_it,
		                                                 ++out_it) {
		    *out_it += ET( art_s[ *row_it]) * in_value;
		}
	    }
	}
    }
}

// check basis inverse
template < class Rep_ >
bool
QPE_solver<Rep_>::
check_basis_inverse()
{
    // diagnostic output
    CGAL_qpe_debug {
	vout4 << "check: " << std::flush;
    }
    bool ok;
    if (is_phaseI) {
    	ok = check_basis_inverse(Tag_true());
    } else {
        ok = check_basis_inverse( Is_linear());
    }

    // diagnostic output
    CGAL_qpe_debug {
	if ( ok) {
	    vout4 << "check ok";
	} else {
	    vout4 << "check failed";
	}
	vout4 << std::endl;
    }

    return ok;
}

template < class Rep_ >                                         // LP case
bool  QPE_solver<Rep_>::
check_basis_inverse( Tag_true)
{
    CGAL_qpe_debug {
	vout4 << std::endl;
    }
    bool res = true;
    unsigned int    row, rows =   C.size();
    unsigned int    col, cols = B_O.size();
    Index_iterator  i_it = B_O.begin();
    Value_iterator  q_it;
    for ( col = 0; col < cols; ++col, ++i_it) {
	ratio_test_init__A_Cj( tmp_l.begin(), *i_it, Has_no_inequalities());
	inv_M_B.multiply_x( tmp_l.begin(), q_x_O.begin());

	CGAL_qpe_debug {
	    if ( vout4.verbose()) {
		std::copy( tmp_l.begin(), tmp_l.begin()+rows,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << " ||  ";
		std::copy( q_x_O.begin(), q_x_O.begin()+cols,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << std::endl;
	    }
	}

	q_it = q_x_O.begin();
	for ( row = 0; row < rows; ++row, ++q_it) {
	    if ( *q_it != ( row == col ? d : et0)) {
		if ( ! vout4.verbose()) {
		    std::cerr << std::endl << "basis-inverse check: ";
		}
		std::cerr << "failed ( row=" << row << " | col=" << col << " )"
		          << std::endl;
		res = false;
	    }
	}
    }

    return res;
}

template < class Rep_ >                                         // QP case
bool  QPE_solver<Rep_>::
check_basis_inverse( Tag_false)
{
    bool res = true;
    unsigned int    row, rows =   C.size();
    unsigned int    col, cols = B_O.size();
    Value_iterator  v_it;
    Index_iterator  i_it;

    CGAL_qpe_debug {
	vout4 << std::endl;
    }

    // left part of M_B
    std::fill_n( tmp_l.begin(), rows, et0);
    for ( col = 0; col < rows; ++col) {

	// get column of A_B^T (i.e. row of A_B)
	row = ( has_ineq ? C[ col] : col);
	v_it = tmp_x.begin();
	for ( i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) {
	    *v_it = ( *i_it < qp_n ? qp_A[ *i_it][ row] :           // original
		      art_A[ *i_it - qp_n].first != (int)row ? et0 :// artific.
		      ( art_A[ *i_it - qp_n].second ? -et1 : et1));
	}
	
	if ( art_s_i >= 0) {              // special artificial variable basic?
	    int  k = qp_n+slack_A.size()+art_s_i;
	    if ( in_B[ k] >= 0) {
		tmp_x[ in_B[ k]] = art_s[ row];
	    }
	}
	
	inv_M_B.multiply( tmp_l.begin(), tmp_x.begin(),
			  q_lambda.begin(), q_x_O.begin());

	CGAL_qpe_debug {
	    if ( vout4.verbose()) {
		std::copy( tmp_l.begin(), tmp_l.begin()+rows,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << "| ";
		std::copy( tmp_x.begin(), tmp_x.begin()+cols,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << " ||  ";
		std::copy( q_lambda.begin(), q_lambda.begin()+rows,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << " |  ";
		std::copy( q_x_O.begin(), q_x_O.begin()+cols,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << std::endl;
	    }
	}

	v_it = q_lambda.begin();
	for ( row = 0; row < rows; ++row, ++v_it) {
	    if ( *v_it != ( row == col ? d : et0)) {
		if ( ! vout4.verbose()) {
		    std::cerr << std::endl << "basis-inverse check: ";
		}
		std::cerr << "failed ( row=" << row << " | col=" << col << " )"
		          << std::endl;
		//		return false;
		res = false;
	    }
	}
	v_it = std::find_if( q_x_O.begin(), q_x_O.begin()+cols,
			     std::bind2nd( std::not_equal_to<ET>(), et0));
	if ( v_it != q_x_O.begin()+cols) {
	    if ( ! vout4.verbose()) {
		std::cerr << std::endl << "basis-inverse check: ";
	    }
	    std::cerr << "failed ( row=" << rows+(v_it-q_x_O.begin())
		      << " | col=" << col << " )" << std::endl;
	    // ToDo: return false;
	    res = false;
	}
    }
    vout4 << "= = = = = = = = = =" << std::endl;

    // right part of M_B
    if ( is_phaseI) std::fill_n( tmp_x.begin(), B_O.size(), et0);
    i_it = B_O.begin();
    for ( col = 0; col < cols; ++col, ++i_it) {
	ratio_test_init__A_Cj  ( tmp_l.begin(), *i_it, Has_no_inequalities());
	ratio_test_init__2_D_Bj( tmp_x.begin(), *i_it, Tag_false());
	inv_M_B.multiply( tmp_l.begin(), tmp_x.begin(),
			  q_lambda.begin(), q_x_O.begin());

	CGAL_qpe_debug {
	    if ( vout4.verbose()) {
		std::copy( tmp_l.begin(), tmp_l.begin()+rows,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << "| ";
		std::copy( tmp_x.begin(), tmp_x.begin()+cols,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << " ||  ";
		std::copy( q_lambda.begin(), q_lambda.begin()+rows,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << " |  ";
		std::copy( q_x_O.begin(), q_x_O.begin()+cols,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << std::endl;
	    }
	}

	v_it = std::find_if( q_lambda.begin(), q_lambda.begin()+rows,
			     std::bind2nd( std::not_equal_to<ET>(), et0));
	if ( v_it != q_lambda.begin()+rows) {
	    if ( ! vout4.verbose()) {
		std::cerr << std::endl << "basis-inverse check: ";
	    }
	    std::cerr << "failed ( row=" << v_it-q_lambda.begin()
		      << " | col=" << col << " )" << std::endl;
	    //	    return false;
	    res = false;
	}
	v_it = q_x_O.begin();
	for ( row = 0; row < cols; ++row, ++v_it) {
	    if ( *v_it != ( row == col ? d : et0)) {
		if ( ! vout4.verbose()) {
		    std::cerr << std::endl << "basis-inverse check: ";
		}
		std::cerr << "failed ( row=" << row+rows << " | col="
			  << col << " )" << std::endl;
		//		return false;
		res = false;
	    }
	}
    }
    return res;
}

// setting the pricing strategy
template < class Rep_ >
void  QPE_solver<Rep_>::
set_pricing_strategy( Pricing_strategy& strategy)
{
    CGAL_qpe_precondition( phase() != 1);
    CGAL_qpe_precondition( phase() != 2);

    strategyP = &strategy;
    if ( phase() != -1) strategyP->set( *this, vout2);
}

// diagnostic output
// -----------------
template < class Rep_ >
void  QPE_solver<Rep_>::
set_verbosity( int verbose, std::ostream& stream)
{
    vout  = Verbose_ostream( verbose >  0, stream);
    vout1 = Verbose_ostream( verbose == 1, stream);
    vout2 = Verbose_ostream( verbose >= 2, stream);
    vout3 = Verbose_ostream( verbose >= 3, stream);
    vout4 = Verbose_ostream( verbose == 4, stream);
}

template < class Rep_ >
void  QPE_solver<Rep_>::
print_program( )
{
    int  row, i;

    // objective function
    vout4.out() << std::endl << "objective function:" << std::endl;
    if ( is_QP) {
	vout4.out() << "--> output of MATRIX must go here <--";
	vout4.out() << std::endl;
    }
    std::copy( qp_c, qp_c+qp_n,
	       std::ostream_iterator<C_entry>( vout4.out(), " "));
    vout4.out() << std::endl;
    vout4.out() << std::endl;

    // constraints
    vout4.out() << "constraints:" << std::endl;
    for ( row = 0; row < qp_m; ++row) {

	// original variables
	for ( i = 0; i < qp_n; ++i) vout4.out() << qp_A[ i][ row] << ' ';

	// slack variables
	if ( ! slack_A.empty()) {
	    vout4.out() << " |  ";
	    for ( i = 0; i < (int)slack_A.size(); ++i) {
		vout4.out() << ( slack_A[ i].first != row ? " 0" :
		               ( slack_A[ i].second ? "-1" : "+1")) << ' ';
	    }
	}

	// artificial variables
	if ( ! art_A.empty()) {
	    vout4.out() << " |  ";
	    for ( i = 0; i < (int)art_A.size(); ++i) {
		vout4.out() << ( art_A[ i].first != row ? " 0" :
		               ( art_A[ i].second ? "-1" : "+1")) << ' ';
	    }
	}
	if ( ! art_s.empty()) vout4.out() << " |  " << art_s[ row] << ' ';

	// rhs
	vout4.out() << " |  "
		    << ( qp_r[ row] == Rep::EQUAL      ? ' ' :
		       ( qp_r[ row] == Rep::LESS_EQUAL ? '<' : '>')) << "=  "
		    << qp_b[ row] << std::endl;
    }
}

template < class Rep_ >
void  QPE_solver<Rep_>::
print_basis( )
{
    vout1 << "  basis: ";
    vout2 << "basic variables" << ( has_ineq ? "  " : "") << ":  ";
    std::copy( B_O.begin(), B_O.end  (),
	       std::ostream_iterator<int>( vout.out(), " "));
    if ( vout2.verbose()) {
	if ( has_ineq && ( ! slack_A.empty())) {
	    vout2.out() << " |  ";
	    std::copy( B_S.begin(), B_S.end(),
		       std::ostream_iterator<int>( vout2.out(), " "));
	}
	if ( is_phaseI) {
	    vout2.out() << " (artificial: " << art_basic << ')';
	}
	if ( has_ineq) {
	    vout2.out() << std::endl
			<< "basic constraints:  ";
	    std::copy( C.begin(), C.begin()+qp_m-slack_A.size(),
		       std::ostream_iterator<int>( vout2.out(), " "));
	    if ( ! slack_A.empty()) {
		vout2.out() << " |  ";
		std::copy( C.begin()+qp_m-slack_A.size(), C.end(),
			   std::ostream_iterator<int>( vout2.out(), " "));
	    }
	}
	if ( vout3.verbose()) {
	    vout3.out() << std::endl
			<< std::endl
			<< "    in_B: ";
	    std::copy( in_B.begin(), in_B.end  (),
		       std::ostream_iterator<int>( vout3.out(), " "));
	    if ( has_ineq) {
		vout3.out() << std::endl
			    << "    in_C: ";
		std::copy( in_C.begin(), in_C.end  (),
			   std::ostream_iterator<int>( vout3.out(), " "));
	    }
	}
    }
    vout.out() << std::endl;
    vout4 << std::endl << "basis-inverse:" << std::endl;
}

template < class Rep_ >
void  QPE_solver<Rep_>::
print_solution( )
{
    if ( vout3.verbose()) {
	vout3.out() << std::endl
		    << "     b_C: ";
	std::copy( b_C.begin(), b_C.begin()+C.size(),
		   std::ostream_iterator<ET>( vout3.out()," "));
	vout3.out() << std::endl
		    << "  -c_B_O: ";
	std::copy( minus_c_B.begin(), minus_c_B.begin()+B_O.size(),
		   std::ostream_iterator<ET>( vout3.out()," "));
	vout3.out() << std::endl;
    }
    if ( vout2.verbose()) {
	vout2.out() << std::endl
		    << "  lambda: ";
	std::copy( lambda.begin(), lambda.begin()+C.size(),
		   std::ostream_iterator<ET>( vout2.out(), " "));
	vout2.out() << std::endl
		    << "   x_B_O: ";
	std::copy( x_B_O.begin(), x_B_O.begin()+B_O.size(),
		   std::ostream_iterator<ET>( vout2.out(), " "));
	vout2.out() << std::endl;
	if ( has_ineq) {
	    vout2.out() << "   x_B_S: ";
	    std::copy( x_B_S.begin(), x_B_S.begin()+B_S.size(),
		       std::ostream_iterator<ET>( vout2.out()," "));
	    vout2.out() << std::endl;
	}
	vout2.out() << std::endl;
    }
    Quotient<ET>  s = solution();
    vout1 << "  ";
    vout.out() << "solution: " << s << "  ~= " << to_double( s) << std::endl;
    vout2 << std::endl;
}

template < class Rep_ >
const char*  QPE_solver<Rep_>::
variable_type( int k) const
{
    return ( k <        qp_n                 ? "original"  :
	   ( k < (int)( qp_n+slack_A.size()) ? "slack"     :
	                                       "artificial"));
}


CGAL_END_NAMESPACE

// ===== EOF ==================================================================
