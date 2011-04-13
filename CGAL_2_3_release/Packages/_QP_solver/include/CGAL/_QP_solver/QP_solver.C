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
// file          : include/CGAL/_QP_solver/QP_solver.C
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.5
// revision_date : 2000/09/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Solver for Quadratic Programs
// ============================================================================
                                                                               
CGAL_BEGIN_NAMESPACE
                    
// Class Implementation (continued)
// ================================

// initialization
// --------------

// set-up of QP
template < class Rep_ >
void
QP_solver<Rep_>::
set( int n, int m, int max_b,
     typename QP_solver<Rep_>::A_iterator A_it,
     typename QP_solver<Rep_>::B_iterator b_it,
     typename QP_solver<Rep_>::C_iterator c_it,
     typename QP_solver<Rep_>::D_iterator D_it)
{
    
    CGAL_optimisation_debug {
        vout2 << std::endl
              << "------" << std::endl
              << "Set-Up" << std::endl
              << "------" << std::endl;
        vout  << "[ " << ( CGAL::check_tag( Is_lp()) ? "L" : "Q") << "P, "
                      << n << " variables, "
                      << m << " constraints ]" << std::endl;
    }
     

    // store QP
    CGAL_optimisation_precondition( m >  0);
    CGAL_optimisation_precondition( n >= m);
    qp_n = n;    qp_m = m;
    qp_A = A_it; qp_b = b_it; qp_c = c_it; qp_D = D_it;
    max_basis = max_b;

    
    // set up pricing strategy
    strategyP->set( *this, vout2);
                                  
}

// set-up of auxiliary problem
template < class Rep_ >
void
QP_solver<Rep_>::
set_up_auxiliary_problem( )
{
    int i;

    // delete artifical part of `A' and auxiliary `c', if necessary
    art_A.erase( art_A.begin(), art_A.end());
    aux_c.erase( aux_c.begin(), aux_c.end());

    // initialize artificial part of `A' and auxiliary `c'
    art_A.reserve( qp_m);
    aux_c.reserve( qp_n+qp_m);
    aux_c.insert( aux_c.end(), qp_n, nt_0);

    for ( i = 0; i < qp_m; ++i) {
        Artificial_column  art_col;
        art_col.push_back( std::make_pair( i,
                           qp_b[ i] >= nt_0 ? nt_1 : nt_minus_1));
        art_A.push_back( art_col);
        aux_c.push_back( nt_1);
    }

    // handling of zero `b_i's, if any
    if ( std::find( qp_b, qp_b+qp_m, nt_0) != qp_b+qp_m) {
        int  k = 0;
        while ( ( k < qp_m) && ( qp_b[ k] == nt_0)) { ++k; }
        CGAL_optimisation_precondition( k < qp_m);
        for ( i = 0; i < qp_m; ++i) {
            if ( qp_b[ i] == nt_0) {
                art_A[ k].push_back( std::make_pair( i, nt_minus_1));
            }
        }
    }
}

// set-up of basis and basis inverse
template < class Rep_ >
void
QP_solver<Rep_>::
set_up_basis( )
{
    int i;

    // initialize basis (with artificial variables)
    if ( B.size() > 0) B.erase( B.begin(), B.end());
    B.insert( B.end(), qp_m, 0);
    for ( i = 0; i < qp_m; ++i) B[ i] = qp_n+i;
    art_basic = qp_m;

    
    CGAL_optimisation_debug {
        vout1 << "  "; vout2 << "initial ";
        
        vout << "basis: ";
        if ( vout.verbose()) {
            std::copy( B.begin(), B.end(),
                       std::ostream_iterator<int>( vout.out(), " "));
            vout.out() << std::endl;
        }
         
        vout3 << "initial basis-inverse:" << std::endl;
    }
     

    // initialize positions in basis
    if ( in_B.size() > 0) in_B.erase( in_B.begin(), in_B.end());
    in_B.reserve( qp_n+qp_m);
    in_B.insert( in_B.end(), qp_n, -1);
    for ( i = 0; i < qp_m; ++i) in_B.push_back( i);

    // initialize basis inverse
    int  k = 0;
    while ( qp_b[ k] == nt_0) { ++k; }

    std::vector<NT>  u, w;
    u.reserve( qp_m);
    w.reserve( qp_m);
    NT  u_k = ( qp_b[ k] > nt_0 ? nt_1 : nt_minus_1);
    /***** NT  u_k = nt_0; *****/
    for ( i = 0; i < qp_m; ++i) {
        if ( qp_b[ i] < nt_0) {
            u.push_back( nt_minus_1);
            w.push_back( nt_0);
        } else {
            u.push_back( nt_1);
            w.push_back( qp_b[ i] > nt_0 ? nt_0 : u_k);
        }
    }

    inv_M_B.init( qp_m, k, u.begin(), w.begin());
}

// set-up of initial solution
template < class Rep_ >
void
QP_solver<Rep_>::
set_up_initial_solution( )
{
    // initialize exact version of `qp_b' (implicit conversion to ET)
    b = Values( qp_b, qp_b+qp_m);

    // initialize exact version of `-aux_c' (implicit conversion to ET)
    minus_c_B.erase( minus_c_B.begin(), minus_c_B.end());
    minus_c_B.reserve( qp_m);
    std::transform( aux_c.begin()+qp_n, aux_c.end(),
                    std::back_inserter( minus_c_B), std::negate<ET>());

    // allocate storage
    lambda = Values( qp_m);
    x_B    = Values( qp_m);

    // compute initial solution
    inv_M_B.multiply( b.begin(), minus_c_B.begin(),
                      lambda.begin(), x_B.begin());
}

// set-up of additional variables
template < class Rep_ >
void
QP_solver<Rep_>::
set_up_additional_variables( )
{
    
    A_j      = Values( qp_m);
    two_D_Bj = Values( qp_m);
    
    q_lambda = Values( qp_m);
    q_x      = Values( qp_m);
    
    
    strategyP->init();
}

// initialization of phase I
template < class Rep_ >
void
QP_solver<Rep_>::
init( )
{
    
    CGAL_optimisation_debug {
        vout2 << std::endl
              << "--------------" << std::endl
              << "Initialization" << std::endl
              << "--------------" << std::endl;
    }
     

    // set up auxiliary problem
    set_up_auxiliary_problem();

    // set up basis and basis inverse
    set_up_basis();

    // set up initial solution
    set_up_initial_solution();

    // set up additional variables
    set_up_additional_variables();

    // set up status
    m_phase      = 1;
    m_status     = UPDATE;
    m_iterations = 0;
    is_phase_I   = true;

    
    CGAL_optimisation_debug {
        vout2 << std::endl;
        
        vout2 << "     b: ";
        if ( vout2.verbose()) {
            std::copy( b.begin(), b.begin()+qp_m,
                       std::ostream_iterator<ET>( vout2.out(), " "));
            vout2.out() << std::endl;
        }
         
        
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
         
        vout1 << "  "; vout2 << std::endl << "initial ";
        
        vout << "solution: ";
        CGAL::Quotient<ET>  s = solution();
        vout  << s << "  ~= " << CGAL::to_double( s) << std::endl;
                                                                  
    }
     
}

// transition to phase II
template < class Rep_ >
void
QP_solver<Rep_>::
transition( )
{
    
    CGAL_optimisation_debug {
        vout1 << "  t"; vout2 << std::endl << "T";
        vout  <<  "ransition to phase II" << std::endl;
        vout2 << "----------------------" << std::endl;
        }
         

    // remove artificial variables
    in_B.erase( in_B.begin()+qp_n, in_B.end());

    // initialize exact version of `-qp_c' (implicit conversion to ET)
    Access_c_B  access_c_B( qp_c);
    std::transform( c_B_iterator( B.begin(), access_c_B),
                    c_B_iterator( B.end  (), access_c_B),
                    minus_c_B.begin(), std::negate<ET>());

    // compute initial solution of phase II
    inv_M_B.multiply( b.begin(), minus_c_B.begin(),
                      lambda.begin(), x_B.begin());

    // update status
    m_phase    = 2;
    is_phase_I = false;

    // notify pricing strategy
    strategyP->transition();

    
    CGAL_optimisation_debug {
        
        vout2 << "     b: ";
        if ( vout2.verbose()) {
            std::copy( b.begin(), b.begin()+qp_m,
                       std::ostream_iterator<ET>( vout2.out(), " "));
            vout2.out() << std::endl;
        }
         
        
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
         
        vout1 << "  "; vout2 << std::endl;
        
        vout << "solution: ";
        CGAL::Quotient<ET>  s = solution();
        vout  << s << "  ~= " << CGAL::to_double( s) << std::endl;
                                                                  
    }

    CGAL_optimisation_assertion( check_basis( Is_lp()));
}


// access functions
// ----------------

// numerator of current solution
template < class Rep_ >
QP_solver<Rep_>::ET
QP_solver<Rep_>::
solution_numerator( ) const
{
    Basic_variable_index_iterator        i_it,   j_it;
    Basic_variable_numerator_iterator  x_i_it, x_j_it;
    ET   s = et_0, sum;
    int  i;
    bool is_phase_I = (phase() == 1);

    // compute  c^T x + x^T D x  (D is symmetric)
    // ------------------------------------------

    // i: 0..|B|-1
      i_it = basic_variables_index_begin();
    x_i_it = basic_variables_numerator_begin();
    for ( ; i_it != basic_variables_index_end(); ++i_it, ++x_i_it) {
        sum = et_0;
        i   = *i_it;

        if ( ! ( CGAL::check_tag( Is_lp()) || is_phase_I)) {

            // j: 0..i-1
              j_it = basic_variables_index_begin();
            x_j_it = basic_variables_numerator_begin();
            for ( ; j_it != i_it; ++j_it, ++x_j_it) {

                // D_{i,j} x_j
                sum += ET( qp_D[ i][ *j_it]) * *x_j_it;
            }
            sum *= et_2;

            // D_{i,i} x_i
            sum += ET( qp_D[ i][ i]) * *x_i_it;
        }

        // d c_i
        sum += d * ( is_phase_I ? aux_c[ i] : qp_c[ i]);

        s += sum * *x_i_it;
    }
    return s;
}

// variables of dual LP
template < class Rep_ >
QP_solver<Rep_>::ET
QP_solver<Rep_>::
dual_variable( int i) const
{
    Assert_compile_time_tag( Tag_true(), Is_lp());
    Values  unity( qp_m), col;
    unity[ i] = et_1;
    col.reserve( qp_m);
    inv_M_B.multiply_x( unity.begin(), std::back_inserter( col));
    return std::inner_product( col.begin(), col.end(),
                               minus_c_B.begin(), et_0);
}


// pivot function
// --------------
template < class Rep_ >
void
QP_solver<Rep_>::
pivot_step( )
{
    ++m_iterations;

    
    CGAL_optimisation_debug {
        vout2 << std::endl
              << "----------" << std::endl
              << "Pivot Step" << std::endl
              << "----------" << std::endl;
        vout  << "[ phase " << ( is_phase_I ? "I" : "II")
              << ", iteration " << m_iterations << " ]" << std::endl;
    }
     

    // pricing
    // -------
    pricing();

    // check for optimality
    if ( 
         j < 0
              ) {

        // which phase?
        if ( is_phase_I) {          // phase I

            // check for infeasibility
            if ( 
                 art_basic > 0
                              ) {
                m_phase  = 3;
                m_status = INFEASIBLE;
                
                CGAL_optimisation_debug {
                   vout1 << "  "; vout << "problem is INFEASIBLE" << std::endl;
                }
                 
            } else {
                // QP feasible, transition to phase II
                transition();
            }
        } else {                    // phase II
            m_phase  = 3;
            m_status = OPTIMAL;
            
            CGAL_optimisation_debug {
                vout1 << "  "; vout << "solution is OPTIMAL" << std::endl;
            }
             
        }
        return;
    }

    /*
    // loop until new basis is found
    init_ratio_test_update_loop();
    do {

        // ratio test
        // ----------
        ratio_test();

        // check for unboundedness
        if ( 
             q_i == et_0
                        ) {
            m_phase  = 3;
            m_status = UNBOUNDED;
            
            CGAL_optimisation_debug {
                vout1 << "  "; vout << "problem is UNBOUNDED" << std::endl;
            }
             
            return;
        }

        // update
        // ------
        update();

    } while ( 
              j >= 0
                    );
    */

    // ratio test
    // ----------
    ratio_test();

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

    // update
    // ------
    update();

    // loop until new basis is found
    // -----------------------------
    iterated_ratio_test_update();

    
    CGAL_optimisation_debug {
        vout1 << std::endl << "  ";
        
        vout1 << "basis: ";
        if ( vout1.verbose()) {
            std::copy( B.begin(), B.end(),
                       std::ostream_iterator<int>( vout1.out(), " "));
            vout1.out() << std::endl;
        }
         
        vout1 << "  "; vout2 << std::endl << "new ";
        
        vout << "solution: ";
        CGAL::Quotient<ET>  s = solution();
        vout  << s << "  ~= " << CGAL::to_double( s) << std::endl;
                                                                  
    }
     
}
 
CGAL_END_NAMESPACE

// ===== EOF ==================================================================
