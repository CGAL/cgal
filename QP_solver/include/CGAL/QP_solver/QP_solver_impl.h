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

#include <CGAL/QP_solver/Initialization.h>

namespace CGAL {

// =============================
// class implementation (cont'd)
// =============================

// transition (to phase II)
// ------------------------
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
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
    //ensure_size(tmp_x_2, tmp_x.size());
    // update basis inverse
    CGAL_qpe_debug {
	vout4 << std::endl << "basis-inverse:" << std::endl;
    }
    transition( Is_linear());
    CGAL_qpe_debug {
        check_basis_inverse();
    }

    // initialize exact version of `-qp_c' (implicit conversion to ET)
    C_by_index_accessor  c_accessor( qp_c);
    std::transform( C_by_index_iterator( B_O.begin(), c_accessor),
                    C_by_index_iterator( B_O.end  (), c_accessor),
                    minus_c_B.begin(), std::negate<ET>());
    
    // compute initial solution of phase II
    compute_solution(Is_nonnegative());

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout.verbose()) print_solution();
    }

    // notify pricing strategy
    strategyP->transition();
}

// access
// ------
// numerator of current solution; the denominator is 2*d*d, so we should
// compute here d*d*(x^T2Dx + x^T2c + 2c0)
template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::
solution_numerator( ) const
{
    ET   s, z = et0;
    int  i, j;

    if (check_tag(Is_nonnegative()) || is_phaseI) {
      // standard form or phase I; it suffices to go
      // through the basic variables; all D- and c-entries
      // are obtained through the appropriate iterators 
      Index_const_iterator  i_it;
      Value_const_iterator  x_i_it, c_it;

      // foreach i
      x_i_it =       x_B_O.begin();
      c_it = minus_c_B  .begin();
      for ( i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++x_i_it, ++c_it){
        i = *i_it;

        // compute quadratic part: 2D_i x
        s = et0;
        if ( is_QP && is_phaseII) {
	  // half the off-diagonal contribution
	  s += std::inner_product(x_B_O.begin(), x_i_it,
				  D_pairwise_iterator(
					 B_O.begin(),
					 D_pairwise_accessor( qp_D, i)),
				  et0);
	  // the other half
	  s *= et2;
	  // diagonal contribution
	  s += ET( (*(qp_D + i))[ i]) * *x_i_it;
        }
        // add linear part: 2c_i
        s -= d * et2 * ET( *c_it);

        // add x_i(2D_i x + 2c_i)
        z += s * *x_i_it; // endowed with a factor of d*d now
      }
    } else {
      // nonstandard form and phase II, 
      // take all original variables into account; all D- and c-entries
      // are obtained from the input data directly
      // order in i_it and j_it matches original variable order
      if (is_QP) {
	// quadratic part
	i=0;
	for (Variable_numerator_iterator 
	       i_it = this->original_variables_numerator_begin(); 
	     i_it < this->original_variables_numerator_end(); ++i_it, ++i) {
	  // do something only if *i_it != 0
	  if (*i_it == et0) continue;
	  s = et0; // contribution of i-th row
	  Variable_numerator_iterator j_it = 
	    this->original_variables_numerator_begin();
	  // half the off-diagonal contribution
	  j=0;
	  for (; j<i; ++j_it, ++j)
	    s += ET(*((*(qp_D+i))+j)) * *j_it;
	  // the other half
	  s *= et2;
	  // the diagonal entry
	  s += ET(*((*(qp_D+i))+j)) * *j_it;
	  // accumulate
	  z += s * *i_it;
	}
      }
      // linear part
      j=0; s = et0;
      for (Variable_numerator_iterator 
	     j_it = this->original_variables_numerator_begin();
	   j_it < this->original_variables_numerator_end(); ++j_it, ++j)
	s +=  et2 * ET(*(qp_c+j)) * *j_it;
      z += d * s;
    }
    // finally, add the constant term (phase II only)
    if (is_phaseII) z += et2 * ET(qp_c0) * d * d;
    return z;
}

// pivot step
// ----------
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
pivot_step( )
{
    ++m_pivots;

    
  
    // diagnostic output
    CGAL_qpe_debug {
        vout2 << std::endl
              << "==========" << std::endl
              << "Pivot Step" << std::endl
              << "==========" << std::endl;
    }
  
  
    vout  << "[ phase " << ( is_phaseI ? "I" : "II")
	  << ", iteration " << m_pivots << " ]" << std::endl;
    
	    
    // pricing
    // -------
    pricing();

	    
    // check for optimality
    if ( j < 0) {

        if ( is_phaseI) {                               // phase I
	    // since we no longer assume full row rank and subsys assumption
	    // we have to strengthen the precondition for infeasibility
            if (this->solution_numerator() > et0) {    
	      // problem is infeasible
	        m_phase  = 3;
	        m_status = QP_INFEASIBLE;
	        
		vout << "  ";
		vout << "problem is INFEASIBLE" << std::endl;
	       
	    } else {  // Drive/remove artificials out of basis
	        expel_artificial_variables_from_basis();
	        transition();
	    }
        } else {                                        // phase II

	    // optimal solution found
	    m_phase  = 3;
            m_status = QP_OPTIMAL;
  
	    vout << "  ";
	    vout  << "solution is OPTIMAL" << std::endl;
            
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
            m_status = QP_UNBOUNDED;
            
	    vout << "  ";
	    vout << "problem is UNBOUNDED" << std::endl;
	    
	    CGAL_qpe_debug {
		//nu should be zero in this case
		// note: (-1)/hat{\nu} is stored instead of \hat{\nu}
		// todo kf: as this is just used for an assertion check,
		// all the following lines should only be executed if
		// assertions are enabled...
		nu = inv_M_B.inner_product(     A_Cj.begin(), two_D_Bj.begin(),
		    q_lambda.begin(),    q_x_O.begin());
	        if (is_QP) {
		    if (j < qp_n) {
		        nu -= d*ET(*((*(qp_D+j))+j));
		    }
		}
		CGAL_qpe_assertion(nu == et0);
            }
            return;
        }
	
        // update
        update_1();

    } while ( j >= 0);

    // ratio test & update (step 2)
    // ----------------------------
/*    
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
*/
    // instead of the above piece of code we now have
    // diagnostic output
    if (is_RTS_transition) {
        is_RTS_transition = false;
     
        CGAL_qpe_debug {
            vout2 << std::endl
		  << "----------------------------" << std::endl
		  << "Ratio Test & Update (Step 2)" << std::endl
		  << "----------------------------" << std::endl;
        }

        // compute index of entering variable
        j += static_cast<int>(in_B.size());

        ratio_test_2( Is_linear());
    
        while ((i >= 0) && basis_matrix_stays_regular()) {
        
	    update_2(Is_linear());
	
	    ratio_test_2(Is_linear());
	
        }
    } 



    // ratio test & update (step 3)
    // ----------------------------
    CGAL_qpe_assertion_msg( i < 0, "Step 3 should never be reached!");

    // diagnostic output
    
    if ( vout.verbose()) print_basis();
    if ( vout.verbose()) print_solution();

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
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
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
    j = strategyP->pricing(direction);

    // diagnostic output

    if ( vout.verbose()) {
      if ( j < 0) {
	CGAL_qpe_debug {
	  vout2 << "entering variable: none" << std::endl;
	}
      } else {
	vout  << "  ";
	vout  << "entering: ";
	vout  << j;
	CGAL_qpe_debug {
	  vout2 << " (" << variable_type( j) << ')' << std::endl;
	  vout2 << "direction: "
		<< ((direction == 1) ? "positive" : "negative") << std::endl;
	}
      }
    }
}

// initialization of ratio-test
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
ratio_test_init( )
{
    // store exact version of `A_Cj' (implicit conversion)
    ratio_test_init__A_Cj( A_Cj.begin(), j, no_ineq);

    // store exact version of `2 D_{B_O,j}'
    ratio_test_init__2_D_Bj( two_D_Bj.begin(), j, Is_linear());
}

template < typename Q, typename ET, typename Tags >                                         // no ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j_, Tag_true)
{
    // store exact version of `A_Cj' (implicit conversion)
    if ( j_ < qp_n) {                                   // original variable

	CGAL::cpp11::copy_n( *(qp_A + j_), qp_m, A_Cj_it);

    } else {                                            // artificial variable

	unsigned int  k = j_;
	k -= qp_n;
	std::fill_n( A_Cj_it, qp_m, et0);
	A_Cj_it[ k] = ( art_A[ k].second ? -et1 : et1);
    }
}

template < typename Q, typename ET, typename Tags >                                        // has ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j_, Tag_false)
{
    // store exact version of `A_Cj' (implicit conversion)
    if ( j_ < qp_n) {                                   // original variable
      A_by_index_accessor  a_accessor( *(qp_A + j_));
      std::copy( A_by_index_iterator( C.begin(), a_accessor),
                 A_by_index_iterator( C.end  (), a_accessor),
                 A_Cj_it);

    } else {
      unsigned int  k = j_;
      k -= qp_n;
      std::fill_n( A_Cj_it, C.size(), et0);

      if ( k < static_cast<unsigned int>(slack_A.size())) {                      // slack variable

        A_Cj_it[ in_C[ slack_A[ k].first]] = ( slack_A[ k].second ? -et1
						                      :  et1);

      } else {                                        // artificial variable
        k -= static_cast<unsigned int>(slack_A.size());

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
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
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
    ratio_test_1__q_x_S( no_ineq);

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

    // computation of t_{min}^{j}
    ratio_test_1__t_min_j(Is_nonnegative());
    CGAL_qpe_debug { // todo kf: at first sight, this debug message should
                     // only be output for problems in nonstandard form...
        if (vout2.verbose()) {
            vout2.out() << "t_min_j: " << x_i << '/' << q_i << std::endl;
            vout2.out() << std::endl;
        }
    }    

    // what happens, if all original variables are nonbasic?
/*
    ratio_test_1__t_i(   B_O.begin(),   B_O.end(),
		       x_B_O.begin(), q_x_O.begin(), Tag_false());
    ratio_test_1__t_i(   B_S.begin(),   B_S.end(),
		       x_B_S.begin(), q_x_S.begin(), no_ineq);
*/		       
    ratio_test_1__t_min_B(no_ineq);    

    // check `t_j'
    ratio_test_1__t_j( Is_linear());

    // diagnostic output
    CGAL_qpe_debug {
        if ( vout2.verbose()) {
            for ( unsigned int k = 0; k < static_cast<unsigned int>(B_O.size()); ++k) {
                print_ratio_1_original(k, x_B_O[k], q_x_O[k]);
            }     
            if ( has_ineq) {
                for ( unsigned int k = 0; k < static_cast<unsigned int>(B_S.size()); ++k) {
                    /*
                    vout2.out() << "t_S_" << k << ": "
				    << x_B_S[ k] << '/' << q_x_S[ k]
				    << ( ( q_i > et0) && ( i == B_S[ k]) ? " *":"")
				    << std::endl;
				    */
				    print_ratio_1_slack(k, x_B_S[k], q_x_S[k]);
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
    }
    if ( q_i > et0) {
      if ( i < 0) {
	vout2 << "leaving variable: none" << std::endl;
      } else {
	vout << ", ";
	vout  << "leaving: ";
	vout  << i;
	CGAL_qpe_debug {
	  if ( vout2.verbose()) {
	    if ( ( i < qp_n) || ( i >= static_cast<int>( qp_n+slack_A.size())) ) {
	      vout2.out() << " (= B_O[ " << in_B[ i] << "]: "
			  << variable_type( i) << ')';
	    } else {
	      vout2.out() << " (= B_S[ " << in_B[ i] << "]: slack)";
	    }
	  }
	  vout2 << std::endl;
	}
	
      }
    }
    
}


template < typename Q, typename ET, typename Tags >                         // Standard form
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_min_j(Tag_true /*is_nonnegative*/)
{
}

// By the pricing step we have the following precondition
// direction == +1 => x_O_v_i[j] == (LOWER v ZERO)
// direction == -1 => x_O_v_i[j] == (UPPER v ZERO) 
template < typename Q, typename ET, typename Tags >                         // Upper bounded
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_min_j(Tag_false /*is_nonnegative*/)
{
    if (j < qp_n) {                                 // original variable
        if (direction == 1) {
            if (x_O_v_i[j] == LOWER) {              // has lower bound value
                if (*(qp_fu+j)) {                   // has finite upper bound
                    x_i = (*(qp_u+j) - *(qp_l+j));
                    q_i = et1;
                    i = j;
                    ratio_test_bound_index = UPPER;
                } else {                            // has infinite upper bound
                    x_i = et1;
                    q_i = et0;
                }
            } else {                                // has value zero
                if (*(qp_fu+j)) {                   // has finite upper bound
                    x_i = *(qp_u+j);
                    q_i = et1;
                    i = j;
                    ratio_test_bound_index = UPPER;
                } else {                            // has infinite upper bound
                    x_i = et1;
                    q_i = et0;                    
                }
            }
        } else {                                    // direction == -1
            if (x_O_v_i[j] == UPPER) {              // has upper bound value
                if (*(qp_fl+j)) {                   // has finite lower bound
                    x_i = (*(qp_u+j) - *(qp_l+j));
                    q_i = et1;
                    i = j;
                    ratio_test_bound_index = LOWER;
                } else {                            // has infinite lower bound
                    x_i = et1;
                    q_i = et0;
                }
            } else {                                // has value zero
                if (*(qp_fl+j)) {                   // has finite lower bound
                    x_i = -(*(qp_l+j));
                    q_i = et1;
                    i = j;
                    ratio_test_bound_index = LOWER;
                } else {                            // has infinite lower bound
                    x_i = et1;
                    q_i = et0;
                }
            }
        }
    } else {                                        // slack or artificial var
        x_i = et1;
        q_i = et0;
    }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_min_B(Tag_true  /*has_equalities_only_and_full_rank*/)
{
    ratio_test_1_B_O__t_i(B_O.begin(), B_O.end(), x_B_O.begin(),
                        q_x_O.begin(), Is_nonnegative());
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_min_B(Tag_false /*has_equalities_only_and_full_rank*/)
{
    ratio_test_1_B_O__t_i(B_O.begin(), B_O.end(), x_B_O.begin(),
                        q_x_O.begin(), Is_nonnegative());
    ratio_test_1_B_S__t_i(B_S.begin(), B_S.end(), x_B_S.begin(),
                        q_x_S.begin(), Is_nonnegative());
}    

// ratio test for the basic original variables
template < typename Q, typename ET, typename Tags >                         // Standard form
void  QP_solver<Q, ET, Tags>::
ratio_test_1_B_O__t_i(Index_iterator i_it, Index_iterator end_it,
                    Value_iterator x_it, Value_iterator q_it,
                    Tag_true  /*is_nonnegative*/)
{    
    for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
        test_implicit_bounds_dir_pos(*i_it, *x_it, *q_it, i, x_i, q_i);
    }
}

// ratio test for the basic original variables                    
template < typename Q, typename ET, typename Tags >                         // Upper bounded
void  QP_solver<Q, ET, Tags>::
ratio_test_1_B_O__t_i(Index_iterator i_it, Index_iterator end_it,
                    Value_iterator x_it, Value_iterator q_it,
                    Tag_false /*is_nonnegative*/)
{
    if (is_phaseI) {
        if (direction == 1) {
            for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
                test_mixed_bounds_dir_pos(*i_it, *x_it, *q_it, i, x_i, q_i);
            }
        } else {
            for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
                test_mixed_bounds_dir_neg(*i_it, *x_it, *q_it, i, x_i, q_i);
            }
        }
    } else {
        if (direction == 1) {
            for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
                test_explicit_bounds_dir_pos(*i_it, *x_it, *q_it, i, x_i, q_i);
            }
        } else {
            for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
                test_explicit_bounds_dir_neg(*i_it, *x_it, *q_it, i, x_i, q_i);
            }
        }
    }
}

// ratio test for the basic slack variables
template < typename Q, typename ET, typename Tags >                         // Standard form
void  QP_solver<Q, ET, Tags>::
ratio_test_1_B_S__t_i(Index_iterator i_it, Index_iterator end_it,
                Value_iterator x_it, Value_iterator q_it,
                Tag_true  /*is_nonnegative*/)
{
    for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
        test_implicit_bounds_dir_pos(*i_it, *x_it, *q_it, i, x_i, q_i);
    }
}

// ratio test for the basic slack variables
template < typename Q, typename ET, typename Tags >                         // Upper bounded
void  QP_solver<Q, ET, Tags>::
ratio_test_1_B_S__t_i(Index_iterator i_it, Index_iterator end_it,
                Value_iterator x_it, Value_iterator q_it,
                Tag_false /*is_nonnegative*/)
{
    if (direction == 1) {
        for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
            test_implicit_bounds_dir_pos(*i_it, *x_it, *q_it, i, x_i, q_i);
        }
    } else {
        for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
            test_implicit_bounds_dir_neg(*i_it, *x_it, *q_it, i, x_i, q_i);
        }    
    }
}

// test for one basic variable with implicit bounds only,
// note that this function writes the member variables i, x_i, q_i
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_implicit_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{   
    if (q_k > et0) {
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((x_k * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (x_k * q_min == d_min * q_k)) ) {
        i_min = k;
        d_min = x_k;
        q_min = q_k;
      }
    }
}

// test for one basic variable with implicit bounds only,
// note that this function writes the member variables i, x_i, q_i
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_implicit_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
    if (q_k < et0) {
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((x_k * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (x_k * q_min == -(d_min * q_k))) ) {
        i_min = k;
        d_min = x_k;
        q_min = -q_k;
      }
    }
}

// test for one basic variable with explicit bounds only,
// note that this function writes the member variables i, x_i, q_i and
// ratio_test_bound_index, although the second and third variable name
// are in the context of upper bounding misnomers
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_explicit_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
    if (q_k > et0) {                                // check for lower bound
        if (*(qp_fl+k)) {
            ET  diff = x_k - (d * ET(*(qp_l+k)));
            // BLAND rule: In case the ratios are the same, only update if the new index
            // is smaller. The special artificial variable is always made to leave first.
            if ((diff * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == d_min * q_k)) ) {
                i_min = k;
                d_min = diff;
                q_min = q_k;
                ratio_test_bound_index = LOWER;
            }
        }
    } else {                                        // check for upper bound
        if ((q_k < et0) && (*(qp_fu+k))) {
            ET  diff = (d * ET(*(qp_u+k))) - x_k;
            // BLAND rule: In case the ratios are the same, only update if the new index
            // is smaller. The special artificial variable is always made to leave first.
            if ((diff * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == -(d_min * q_k))) ) {
                i_min = k;
                d_min = diff;
                q_min = -q_k;
                ratio_test_bound_index = UPPER;
            }    
        }
    }
}

// test for one basic variable with explicit bounds only,
// note that this function writes the member variables i, x_i, q_i and
// ratio_test_bound_index, although the second and third variable name
// are in the context of upper bounding misnomers
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_explicit_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
    if (q_k < et0) {                                // check for lower bound
        if (*(qp_fl+k)) {
            ET  diff = x_k - (d * ET(*(qp_l+k)));
            // BLAND rule: In case the ratios are the same, only update if the new index
            // is smaller. The special artificial variable is always made to leave first.
            if ((diff * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == -(d_min * q_k))) ) {
                i_min = k;
                d_min = diff;
                q_min = -q_k;
                ratio_test_bound_index = LOWER;
            }
        }
    } else {                                        // check for upper bound
        if ((q_k > et0) && (*(qp_fu+k))) {
            ET  diff = (d * ET(*(qp_u+k))) - x_k;
            // BLAND rule: In case the ratios are the same, only update if the new index
            // is smaller. The special artificial variable is always made to leave first.
            if ((diff * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == d_min * q_k)) ) {
                i_min = k;
                d_min = diff;
                q_min = q_k;
                ratio_test_bound_index = UPPER;
            }    
        }
    }
}

// test for one basic variable with mixed bounds,
// note that this function writes the member variables i, x_i, q_i and
// ratio_test_bound_index, although the second and third variable name
// are in the context of upper bounding misnomers
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_mixed_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
    if (q_k > et0) {                                // check for lower bound
        if (k < qp_n) {                             // original variable
            if (*(qp_fl+k)) {
                ET  diff = x_k - (d * ET(*(qp_l+k)));
                // BLAND rule: In case the ratios are the same, only update if the new index
                // is smaller. The special artificial variable is always made to leave first.
                if ((diff * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == d_min * q_k))) {
                  i_min = k;
                  d_min = diff;
                  q_min = q_k;
                  ratio_test_bound_index = LOWER;
              } // phase  I II switch
              
            }
        } else {                                    // artificial variable
            // BLAND rule: In case the ratios are the same, only update if the new index
            // is smaller. The special artificial variable is always made to leave first.
            if ((x_k * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (x_k * q_min == d_min * q_k))) {
                i_min = k;
                d_min = x_k;
                q_min = q_k;
            }
        }
    } else {                                        // check for upper bound
        if ((q_k < et0) && (k < qp_n) && *(qp_fu+k)) {
            ET  diff = (d * ET(*(qp_u+k))) - x_k;
            // BLAND rule: In case the ratios are the same, only update if the new index
            // is smaller. The special artificial variable is always made to leave first.
            if ((diff * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == -(d_min * q_k))) ) {
                i_min = k;
                d_min = diff;
                q_min = -q_k;
                ratio_test_bound_index = UPPER;
            }
        }
    }
}

// test for one basic variable with mixed bounds,
// note that this function writes the member variables i, x_i, q_i and
// ratio_test_bound_index, although the second and third variable name
// are in the context of upper bounding misnomers
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_mixed_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
    if (q_k < et0) {                                // check for lower bound
        if (k < qp_n) {                             // original variable
            if (*(qp_fl+k)) {
                ET  diff = x_k - (d * ET(*(qp_l+k)));
                // BLAND rule: In case the ratios are the same, only update if the new index
                // is smaller. The special artificial variable is always made to leave first.
                if ((diff * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == -(d_min * q_k))) ) {
                    i_min = k;
                    d_min = diff;
                    q_min = -q_k;
                    ratio_test_bound_index = LOWER;
                }
            }
        } else {                                    // artificial variable
            // BLAND rule: In case the ratios are the same, only update if the new index
            // is smaller. The special artificial variable is always made to leave first.
            if ((x_k * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (x_k * q_min == -(d_min * q_k))) ) {
                i_min = k;
                d_min = x_k;
                q_min = -q_k;
            }
        }
    } else {                                        // check for upper bound
        if ((q_k > et0) && (k < qp_n) && *(qp_fu+k)) {
            ET  diff = (d * ET(*(qp_u+k))) - x_k;
            // BLAND rule: In case the ratios are the same, only update if the new index
            // is smaller. The special artificial variable is always made to leave first.
            if ((diff * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == d_min * q_k)) ) {
                i_min = k;
                d_min = diff;
                q_min = q_k;
                ratio_test_bound_index = UPPER;
            }
        }
    }
}    


template < typename Q, typename ET, typename Tags >                                         // QP case
void
QP_solver<Q, ET, Tags>::
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
    ratio_test_2__p( no_ineq);
 
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
    
    // Idea here: At this point, the goal is to increase \mu_j until either we
    // become optimal (\mu_j=0), or one of the variables in x^*_\hat{B} drops
    // down to zero.
    //
    // Let us see first how this is done in the standard-form case (where
    // Sven's thesis applies).  Eq. (2.11) in Sven's thesis holds, and by
    // multlying it by $M_\hat{B}^{-1}$ we obtain an equation for \lambda and
    // x^*_\hat{B}.  The interesting equation (the one for x^*_\hat{B}) looks
    // more or less as follows:
    //
    //    x(mu_j)      = x(0) + mu_j      q_it                          (1)
    //
    // where q_it is the vector from (2.12).  In paritcular, for
    // mu_j=mu_j(t_1) (i.e., if we plug the value of mu_j at the beginning of
    // this ratio step 2 into (1)) we have
    //
    //    x(mu_j(t_1)) = x(0) + mu_j(t_1) q_it                          (2)
    //
    // where x(mu_j(t_1)) is the current solution of the solver at this point
    // (i.e., at the beginning of ratio step 2).
    //
    // By subtracting (2) from (1) we can thus eliminate the "unkown" x(0)
    // (which is cheaper than computing it):
    //
    //    x(mu_j) = x(mu_j(t_1)) + (mu_j-mu_j(t_1)) q_it
    //                             ----------------
    //                                  := delta
    //
    // In order to compute for each variable x_k in \hat{B} the value of mu_j
    // for which x_k(mu_j) = 0, we thus evaluate
    //
    //                x(mu_j(t_1))_k
    //    delta_k:= - --------------
    //                    q_it_k
    //
    // The first variable in \hat{B} that hits zero "in the future" is then
    // the one whose delta_k equals
    //
    //    delta_min:= min {delta_k | k in \hat{B} and (q_it)_k < 0 }
    //    
    // So in order to handle the standard-form case, we would compute this
    // minimum.  Once we have delta_min, we need to check whether we get
    // optimal BEFORE a variable drops to zero.  As delta = mu_j - mu_j(t_1),
    // the latter is precisely the case if delta_min >= -mu_j(t_1).
    //
    // (Note: please forget the crap identitiy between (2.11) and (2.12); the
    // notation is misleading.)
    //
    // Now to the nonstandard-form case.
    
    // fw: By definition delta_min >= 0, such that initializing
    // delta_min with -mu_j(t_1) has the desired effect that a basic variable
    // is leaving only if 0 <= delta_min < -mu_j(t_1).
    //
    // The only initialization of delta_min as fraction x_i/q_i that works is
    // x_i=mu_j(t_1); q_i=-1; (see below).
    //
    // Since mu_j(t_1) has been computed in ratio test step 1 we can
    // reuse it.
      
    x_i = mu;                                     // initialize minimum
    q_i = -et1;                                        // with -mu_j(t_1) 

    Value_iterator  x_it = x_B_O.begin();
    Value_iterator  q_it = q_x_O.begin();
    Index_iterator  i_it;
    for ( i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++x_it, ++q_it) {
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ( (*q_it < et0) && (
              (( *x_it * q_i) < ( x_i * *q_it)) ||
              ( (*i_it < i) && (i != art_s_i) && (( *x_it * q_i) == ( x_i * *q_it)) )
            )
          ) {
        i = *i_it; x_i = *x_it; q_i = *q_it;
      }
    }
    x_it = x_B_S.begin();
    q_it = q_x_S.begin();
    for ( i_it = B_S.begin(); i_it != B_S.end(); ++i_it, ++x_it, ++q_it) {
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ( ( *q_it < et0) && (
             (( *x_it * q_i) < ( x_i * *q_it)) ||
             ( (*i_it < i) && (i != art_s_i) && (( *x_it * q_i) == ( x_i * *q_it)) )
            )
          ){
          i = *i_it; x_i = *x_it; q_i = *q_it;
      }
    }

    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    for ( unsigned int k = 0; k < static_cast<unsigned int>(B_O.size()); ++k) {
		vout2.out() << "mu_j_O_" << k << ": - "
			    << x_B_O[ k] << '/' << q_x_O[ k]
			    << ( ( q_i < et0) && ( i == B_O[ k]) ? " *" : "")
			    << std::endl;
	    }
	    for ( unsigned int k = 0; k < static_cast<unsigned int>(B_S.size()); ++k) {
		vout2.out() << "mu_j_S_" << k << ": - "
			    << x_B_S[ k] << '/' << q_x_S[ k]
			    << ( ( q_i < et0) && ( i == B_S[ k]) ? " *" : "")
			    << std::endl;
	    }
	    vout2.out() << std::endl;
	}
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

// update (step 1)
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
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
    CGAL_qpe_assertion(check_basis_inverse());
    
    // check the updated vectors r_C and r_S_B
    CGAL_expensive_assertion(check_r_C(Is_nonnegative()));
    CGAL_expensive_assertion(check_r_S_B(Is_nonnegative()));
    
    // check the vectors r_B_O and w in phaseII for QPs
    CGAL_qpe_debug {
        if (is_phaseII && is_QP) {
            CGAL_expensive_assertion(check_r_B_O(Is_nonnegative()));
            CGAL_expensive_assertion(check_w(Is_nonnegative()));
        }
    }

    // compute current solution
    compute_solution(Is_nonnegative());
    
    // check feasibility 
    CGAL_qpe_debug {
      if (j < 0 && !is_RTS_transition) // todo kf: is this too conservative?
                 // Note: the above condition is necessary because of the
                 // following.  In theory, it is true that the current
                 // solution is at this point in the solver always
                 // feasible. However, the solution has its x_j-entry equal to
                 // the current t from the pricing, and is not zero (in the
                 // standard-form case) or the current lower/upper bound
                 // of the variable (in the non-standard-form case), resp., as
                 // the routines is_solution_feasible() and
                 // is_solution_feasible_for_auxiliary_problem() assume.
	if (is_phaseI) {
	  CGAL_expensive_assertion(
            is_solution_feasible_for_auxiliary_problem());
	} else {
	  CGAL_expensive_assertion(is_solution_feasible());
	}
      else
	vout2 << "(feasibility not checked in intermediate step)" << std::endl;
      CGAL_expensive_assertion(check_tag(Is_nonnegative()) ||
			       r_C.size() == C.size());
    }
	 
}

// update (step 2)
template < typename Q, typename ET, typename Tags >                                         // QP case
void
QP_solver<Q, ET, Tags>::
update_2( Tag_false)
{
    CGAL_qpe_debug {
	vout2 << std::endl
	      << "Update (Step 2)" << std::endl
	      << "---------------" << std::endl;
    }

    // leave variable from basis
    leave_variable();
    CGAL_qpe_debug {
        check_basis_inverse();
    }
    
    // check the updated vectors r_C, r_S_B, r_B_O and w
    CGAL_expensive_assertion(check_r_C(Is_nonnegative()));
    CGAL_expensive_assertion(check_r_S_B(Is_nonnegative()));
    CGAL_expensive_assertion(check_r_B_O(Is_nonnegative()));
    CGAL_expensive_assertion(check_w(Is_nonnegative()));

    // compute current solution
    compute_solution(Is_nonnegative());
}
 
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
expel_artificial_variables_from_basis( )
{
    int row_ind;
    ET r_A_Cj;
    
    CGAL_qpe_debug {
        vout2 << std::endl
	      << "---------------------------------------------" << std::endl
	      << "Expelling artificial variables from the basis" << std::endl
	      << "---------------------------------------------" << std::endl;
    }
  

  
    for (unsigned int i_ = 0; i_ < static_cast<unsigned int>(qp_n + slack_A.size()); ++i_) {
      if (!is_basic(i_)) { 
      ratio_test_init__A_Cj( A_Cj.begin(), i_, no_ineq);
      }
    }

    
    // try to pivot the artificials out of the basis
    // Note that we do not notify the pricing strategy about variables
    // leaving the basis, furthermore the pricing strategy does not
    // know about variables entering the basis.
    // The partial pricing strategies that keep the set of nonbasic vars
    // explicitly are synchronized during transition from phaseI to phaseII 
    for (unsigned int i_ = static_cast<unsigned int>(qp_n + slack_A.size()); i_ < static_cast<unsigned int>(in_B.size()); ++i_) {
      if (is_basic(i_)) { 					// is basic
        row_ind = in_B[i_];  
        
        //CGAL_qpe_assertion(row_ind >= 0);
        
	    
        // determine first possible entering variable,
        // if there is any
        for (unsigned int j_ = 0; j_ < static_cast<unsigned int>(qp_n + slack_A.size()); ++j_) {
	        if (!is_basic(j_)) {  				// is nonbasic
            ratio_test_init__A_Cj( A_Cj.begin(), j_, no_ineq);
                        
            r_A_Cj = inv_M_B.inv_M_B_row_dot_col(row_ind, A_Cj.begin());
                       
            if (r_A_Cj != et0) {
              ratio_test_1__q_x_O(Is_linear());
              i = i_;
              j = j_;
              update_1(Is_linear());
              break;
            } 
          }
        }
      }
    }
  
    if ((art_basic != 0) && no_ineq) {
      // the vector in_C was not used in phase I, but now we remove redundant
      // constraints and switch to has_ineq treatment, hence we need it to
      // be correct at this stage
      for (int i=0; i<qp_m; ++i)
        in_C.push_back(i);
      }
      diagnostics.redundant_equations = (art_basic != 0);

      // now reset the no_ineq and has_ineq flags to match the situation
      no_ineq = no_ineq && !diagnostics.redundant_equations;
      has_ineq = !no_ineq;
    
      // remove the remaining ones with their corresponding equality constraints
      // Note: the special artificial variable can always be driven out of the
      // basis
      for (unsigned int i_ = static_cast<unsigned int>(qp_n + slack_A.size()); i_ < static_cast<unsigned int>(in_B.size()); ++i_) {
        if (in_B[i_] >= 0) {
          i = i_;
          CGAL_qpe_debug {
            vout2 << std::endl
            << "~~> removing artificial variable " << i
            << " and its equality constraint" << std::endl
            << std::endl;
          }
        remove_artificial_variable_and_constraint();
      }
    }
}


// replace variable in basis
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
replace_variable( )
{
    CGAL_qpe_debug {
	vout2 <<   "<--> nonbasic (" << variable_type( j) << ") variable " << j
	      << " replaces basic (" << variable_type( i) << ") variable " << i
	      << std::endl << std::endl;
    }

    // replace variable
    replace_variable( no_ineq);

    // pivot step done
    i = j = -1;
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
replace_variable_original_original( )
{
    // updates for the upper bounded case
    replace_variable_original_original_upd_r(Is_nonnegative());
    
    int  k = in_B[ i];

    // replace original variable [ in: j | out: i ]
    in_B  [ i] = -1;
    in_B  [ j] = k;
       B_O[ k] = j;

    minus_c_B[ k] = 
      ( is_phaseI ? 
	( j < qp_n ? et0 : -aux_c[j-qp_n-slack_A.size()]) : -ET( *(qp_c+ j)));

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

// update of the vector r for U_5 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
replace_variable_original_original_upd_r(Tag_true )
{
}

// update of the vector r for U_5 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
replace_variable_original_original_upd_r(Tag_false )
{
    ET      x_j, x_i;
    
    if (is_artificial(j)) {
        if (!is_artificial(i)) {
            x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+i) : *(qp_u+i);
            update_r_C_r_S_B__i(x_i);
            // update x_O_v_i
            x_O_v_i[i] = ratio_test_bound_index;
        }
    } else {
        x_j = nonbasic_original_variable_value(j);
        if (is_artificial(i)) {
            update_r_C_r_S_B__j(x_j);
        } else {
            x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+i) : *(qp_u+i);
            update_r_C_r_S_B__j_i(x_j, x_i);
            // update x_O_v_i
            x_O_v_i[i] = ratio_test_bound_index;
        }
        // update x_O_v_i
        x_O_v_i[j] = BASIC;
    }
}


template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_slack( )
{
    
    // updates for the upper bounded case
    replace_variable_slack_slack_upd_r(Is_nonnegative()); 
    
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

     b_C[ k] = ET( *(qp_b+ new_row));

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }

    // update basis inverse
    A_row_by_index_accessor  a_accessor =
      boost::bind( A_accessor( qp_A, 0, qp_n), _1, new_row);
    std::copy( A_row_by_index_iterator( B_O.begin(), a_accessor),
	       A_row_by_index_iterator( B_O.end  (), a_accessor),
	       tmp_x.begin());
    if ( art_s_i > 0) {                                 // special artificial
	tmp_x[ in_B[ art_s_i]] = ET( art_s[ new_row]);
    }
    inv_M_B.enter_slack_leave_slack( tmp_x.begin(), k);
}

// update of the vector r for U_6 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_slack_upd_r(Tag_true )
{
}

// update of the vector r for U_6 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_slack_upd_r(Tag_false )
{
    int     sigma_j = slack_A[ j-qp_n].first;
    
    // swap r_gamma_C(sigma_j) in r_C with r_gamma_S_B(sigma_i) in r_S_B
    std::swap(r_C[in_C[sigma_j]], r_S_B[in_B[i]]); 
}


template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_original( )
{
    // updates for the upper bounded case
    replace_variable_slack_original_upd_r(Is_nonnegative()); 
     
    int  k = in_B[ i];

    // leave original variable [ out: i ]
    in_B  [ B_O.back()] = k;
       B_O[ k] = B_O.back();
       in_B  [ i         ] = -1;
       B_O.pop_back();

    minus_c_B[ k] = minus_c_B[ B_O.size()];

    if ( is_phaseI && ( i >= qp_n)) --art_basic;

    // enter slack variable [ in: j ]
    int  old_row = slack_A[ j-qp_n].first;
    in_B  [ j] = static_cast<int>(B_S.size());
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

// update of the vector r for U_8 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_original_upd_r(Tag_true )
{
}

// update of the vector r for U_8 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_original_upd_r(Tag_false )
{
    if (!is_artificial(i)) {
        ET  x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+i) : *(qp_u+i);
        update_r_C_r_S_B__i(x_i);
    }
    
    int     sigma_j = slack_A[ j-qp_n].first;
    
    // append r_gamma_C(sigma_j) from r_C to r_S_B:
    r_S_B.push_back(r_C[in_C[sigma_j]]);
    
    // remove r_gamma_C(sigma_j) from r_C:
    r_C[in_C[sigma_j]] = r_C.back();
    r_C.pop_back();
    
    // update x_O_v_i
    if (!is_artificial(i)) // original and not artificial?
      x_O_v_i[i] = ratio_test_bound_index;
}


template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
replace_variable_original_slack( )
{
    // updates for the upper bounded case
    replace_variable_original_slack_upd_r(Is_nonnegative());
    
    int  k = in_B[ i];

    // enter original variable [ in: j ]

    minus_c_B[ B_O.size()]
      = ( is_phaseI ? 
	  ( j < qp_n ? et0 : -aux_c[j-qp_n-slack_A.size()]) 
	  : -ET( *(qp_c+ j)));
    

    in_B  [ j] = static_cast<int>(B_O.size());
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

     b_C[ C.size()] = ET( *(qp_b+ new_row));
    in_C[ new_row ] = static_cast<int>(C.size());
       C.push_back( new_row);
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }

    // update basis inverse
    A_row_by_index_accessor  a_accessor =
      boost::bind (A_accessor( qp_A, 0, qp_n), _1, new_row);
    std::copy( A_row_by_index_iterator( B_O.begin(), a_accessor),
	       A_row_by_index_iterator( B_O.end  (), a_accessor),
	       tmp_x.begin());
    if ( art_s_i > 0) {                                 // special art.
	tmp_x[ in_B[ art_s_i]] = ET( art_s[ new_row]);
    }
    inv_M_B.enter_original_leave_slack( q_x_O.begin(), tmp_x.begin());
    
}

// update of the vector r for U_7 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
replace_variable_original_slack_upd_r(Tag_true )
{
}

// update of the vector r for U_7 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
replace_variable_original_slack_upd_r(Tag_false )
{
    if (!is_artificial(j)) {
        ET x_j = nonbasic_original_variable_value(j);
        update_r_C_r_S_B__j(x_j);
    }
    
    // append r_gamma_S_B(sigma_i) from r_S_B to r_C
    r_C.push_back(r_S_B[in_B[i]]);
    
    // remove r_gamma_S_B(sigma_i) from r_S_B
    r_S_B[in_B[i]] = r_S_B.back();
    r_S_B.pop_back();
    
    // update x_O_v_i
    if (!is_artificial(j)) {
        x_O_v_i[j] = BASIC;
    }
}


template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
remove_artificial_variable_and_constraint( )
{
    // updates for the upper bounded case
    remove_artificial_variable_and_constraint_upd_r(Is_nonnegative());
    
    int  k = in_B[ i];

    // leave artificial (original) variable [ out: i ]
    in_B  [ B_O.back()] = k;
       B_O[ k] = B_O.back();
       in_B  [ i         ] = -1;
       B_O.pop_back();

    minus_c_B[ k] = minus_c_B[ B_O.size()];

    if ( is_phaseI && ( i >= qp_n)) --art_basic;

    int old_row = art_A[i - qp_n - slack_A.size()].first;

    // leave its equality constraint 
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

// update of the vector r with upper bounding for the removal of an
// artificial variable with its equality constraint, note that we 
// need the headings C before it is updated
template < typename Q, typename ET, typename Tags >                                 // Standard form
void  QP_solver<Q, ET, Tags>::
remove_artificial_variable_and_constraint_upd_r(Tag_true )
{
}

// update of the vector r with upper bounding for the removal of an
// artificial variable with its equality constraint, note that we 
// need the headings C before it is updated
template < typename Q, typename ET, typename Tags >                                 // Upper bounded
void  QP_solver<Q, ET, Tags>::
remove_artificial_variable_and_constraint_upd_r(Tag_false )
{
    int sigma_i = art_A[i - qp_n - slack_A.size()].first;
    
    // remove r_gamma_C(sigma_i) from r_C
    r_C[in_C[sigma_i]] = r_C.back();
    r_C.pop_back();
}

// update that occurs only with upper bounding in ratio test step 1
template < typename Q, typename ET, typename Tags >            
void  QP_solver<Q, ET, Tags>::
enter_and_leave_variable( )
{
    
    CGAL_qpe_assertion((i == j) && (i >= 0));
    
    CGAL_qpe_debug {
	vout2 <<   "<--> nonbasic (" << variable_type( j) << ") variable " << j
	      << " enters and leaves basis" << std::endl << std::endl;
    }

    
    ET diff;
    ET x_j = nonbasic_original_variable_value(j);
    
    if (ratio_test_bound_index == LOWER) {
        diff = x_j - ET(*(qp_l+j));
    } else {
        diff = x_j - ET(*(qp_u+j));
    }
    
    if (is_phaseI) {
        update_r_C_r_S_B__j(diff);
    } else {
        update_w_r_B_O__j(diff);
        update_r_C_r_S_B__j(diff);
    }
    
    x_O_v_i[j] = ratio_test_bound_index;
    
    // notify pricing strategy (it has called enter_basis on i before)
    strategyP->leaving_basis (i);

    // variable entered and left basis
    i = -1; j = -1;
}



// enter variable into basis
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
enter_variable( )
{
  CGAL_qpe_assertion (is_phaseII);
    CGAL_qpe_debug {
	vout2 << "--> nonbasic (" << variable_type( j) << ") variable "
	      << j << " enters basis" << std::endl << std::endl;
    }

    // update basis & basis inverse:
    if (no_ineq || (j < qp_n)) {              // original variable
    
        // updates for the upper bounded case:
        enter_variable_original_upd_w_r(Is_nonnegative());

	// enter original variable [ in: j ]:
	if (minus_c_B.size() <= B_O.size()) { // Note: minus_c_B and the
					      // containers resized in this
					      // if-block are only enlarged
					      // and never made smaller
					      // (while B_O always has the
					      // correct size). We check here
					      // whether we need to enlarge
					      // them.
	  CGAL_qpe_assertion(minus_c_B.size() == B_O.size());
	    minus_c_B.push_back(et0);
	        q_x_O.push_back(et0);
	      tmp_x  .push_back(et0);
	      tmp_x_2.push_back(et0);
	     two_D_Bj.push_back(et0);
	        x_B_O.push_back(et0);
	}
	minus_c_B[B_O.size()] = -ET(*(qp_c+ j)); // Note: B_O has always the
					       // correct size.
	
	in_B[j] = static_cast<int>(B_O.size());
	B_O.push_back(j);

	// diagnostic output
	CGAL_qpe_debug {
	    if (vout2.verbose())
	      print_basis();
	}
	    
	// update basis inverse
	// note: (-1)\hat{\nu} is stored instead of \hat{\nu}
	inv_M_B.enter_original(q_lambda.begin(), q_x_O.begin(), -nu);
	
    } else {                                  // slack variable

        // updates for the upper bounded case:
        enter_variable_slack_upd_w_r(Is_nonnegative());

	// enter slack variable [ in: j ]:
	in_B  [ j] = static_cast<int>(B_S.size());
	   B_S.push_back( j);
	   S_B.push_back( slack_A[ j-qp_n].first);

	// leave inequality constraint [ out: j ]:
	int old_row = slack_A[ j-qp_n].first;
	int k = in_C[old_row];
	
	// reflect change of active constraints heading C in b_C:
	b_C[ k] = b_C[C.size()-1];
		
	   C[ k] = C.back();
	in_C[ C.back()      ] = k;
	in_C[ old_row       ] = -1;
	   C.pop_back();
	
	// diagnostic output:
	CGAL_qpe_debug {
	    if (vout2.verbose())
	      print_basis();
	}

	// update basis inverse:
	inv_M_B.swap_constraint(k);  // swap to back
	inv_M_B.enter_slack();       // drop drop
    }

    // variable entered:
    j -= static_cast<int>(in_B.size());
}

// update of the vectors w and r for U_1 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
enter_variable_original_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_1 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded     
void  QP_solver<Q, ET, Tags>::
enter_variable_original_upd_w_r(Tag_false )
{

    ET      x_j = nonbasic_original_variable_value(j);

    // Note: w needs to be updated before r_C, r_S_B
    update_w_r_B_O__j(x_j);
    update_r_C_r_S_B__j(x_j);
    
    // append w_j to r_B_O
    if (!check_tag(Is_linear())) // (kf.)
      r_B_O.push_back(w[j]);
    
    // update x_O_v_i
    x_O_v_i[j] = BASIC;
}

// update of the vectors w and r for U_3 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
enter_variable_slack_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_3 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded     
void  QP_solver<Q, ET, Tags>::
enter_variable_slack_upd_w_r(Tag_false )
{
    
    int     sigma_j = slack_A[ j-qp_n].first;       
    
    // append r_gamma_C(sigma_j) to r_S_B
    r_S_B.push_back(r_C[in_C[sigma_j]]);
    
    // remove r_gamma_C(sigma_j) from r_C   
    r_C[in_C[sigma_j]] = r_C.back();
    r_C.pop_back();
}

// leave variable from basis
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
leave_variable( )
{
    CGAL_qpe_debug {
	vout2 << "<-- basic (" << variable_type( i) << ") variable "
	      << i << " leaves basis" << std::endl << std::endl;
    }

    // update basis & basis inverse
    int  k = in_B[ i];
    if ( no_ineq || ( i < qp_n)) {                      // original variable
        
        // updates for the upper bounded case
        leave_variable_original_upd_w_r(Is_nonnegative());

	// leave original variable [ out: i ]
	in_B  [ B_O.back()] = k;
	in_B  [ i         ] = -1; 
	//in_B  [ B_O.back()] = k;
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
        
        // updates for the upper bounded case
        leave_variable_slack_upd_w_r(Is_nonnegative());

	// leave slack variable [ out: i ]
	in_B  [ B_S.back()] = k;      // former last var moves to position k
	in_B  [ i         ] = -1;     // i gets deleted
	   B_S[ k] = B_S.back(); B_S.pop_back();
	   S_B[ k] = S_B.back(); S_B.pop_back();

	// enter inequality constraint [ in: i ]
	int new_row = slack_A[ i-qp_n].first;

	A_Cj[ C.size()] = ( j < qp_n ? ET( *((*(qp_A + j))+ new_row)) : et0);

	 b_C[ C.size()] = ET( *(qp_b+ new_row));
	in_C[ new_row ] = static_cast<int>(C.size());
	   C.push_back( new_row);

	// diagnostic output
	CGAL_qpe_debug {
	    if ( vout2.verbose()) print_basis();
	}

	// update basis inverse
	A_row_by_index_accessor  a_accessor =
	  boost::bind (A_accessor( qp_A, 0, qp_n), _1, new_row);
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

// update of the vectors w and r for U_2 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
leave_variable_original_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_2 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
leave_variable_original_upd_w_r(Tag_false )
{

    ET      x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+i) : *(qp_u+i);
    
    // Note: w needs to be updated before r_C, r_S_B
    update_w_r_B_O__i(x_i);
    update_r_C_r_S_B__i(x_i);    
    
    // remove r_beta_O(i) from r_B_O
    if (!check_tag(Is_linear())) { // (kf.)
      r_B_O[in_B[i]] = r_B_O.back();
      r_B_O.pop_back();
    }
    
    // update x_O_v_i
    x_O_v_i[i] = ratio_test_bound_index;
}

// update of the vectors w and r for U_4 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
leave_variable_slack_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_4 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
leave_variable_slack_upd_w_r(Tag_false )
{
    
    // append r_gamma_S_B(sigma_i) to r_C
    r_C.push_back(r_S_B[in_B[i]]);
    
    // remove r_gamma_S_B(sigma_i) from r_S_B
    r_S_B[in_B[i]] = r_S_B.back();
    r_S_B.pop_back();
}


// replace variable in basis QP-case, transition to Ratio Test Step 2
template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
z_replace_variable( )
{
    CGAL_qpe_debug {
	vout2 <<   "<--> nonbasic (" << variable_type( j) << ") variable " << j
	      << " z_replaces basic (" << variable_type( i) << ") variable " << i
	      << std::endl << std::endl;
    }

    // replace variable
    z_replace_variable( no_ineq);

    // pivot step not yet completely done
    i = -1;
    j -= static_cast<int>(in_B.size());
    is_RTS_transition = true;
}


template < typename Q, typename ET, typename Tags >  inline                           // no inequalities
void QP_solver<Q, ET, Tags>::
z_replace_variable( Tag_true)
{
    z_replace_variable_original_by_original();
    strategyP->leaving_basis(i);

}


template < typename Q, typename ET, typename Tags >  inline                          // has inequalities
void QP_solver<Q, ET, Tags>::
z_replace_variable( Tag_false)
{
    // determine type of variables
    bool  enter_original = ( (j < qp_n) || (j >= static_cast<int>( qp_n+slack_A.size())));
    bool  leave_original = ( (i < qp_n) || (i >= static_cast<int>( qp_n+slack_A.size())));

    // update basis and basis inverse
    if ( leave_original) {
        if ( enter_original) {               
	    z_replace_variable_original_by_original();
	} else {                             
	    z_replace_variable_original_by_slack();
	}
    } else {
        if ( enter_original) {
	    z_replace_variable_slack_by_original();
	} else {
	    z_replace_variable_slack_by_slack();
	}
    }
    strategyP->leaving_basis( i);
}


// replacement with precond det(M_{B \setminus \{i\}})=0
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_original( )
{
    // updates for the upper bounded case
    z_replace_variable_original_by_original_upd_w_r(Is_nonnegative());
    
    int  k = in_B[ i];

    // replace original variable [ in: j | out: i ]
    in_B  [ i] = -1;
    in_B  [ j] = k;
       B_O[ k] = j;

    minus_c_B[ k] = -ET( *(qp_c+ j));

    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }
	
    // compute s_delta
    D_pairwise_accessor  d_accessor( qp_D, j);
    ET                   s_delta =d_accessor(j)-d_accessor(i); 
    	    
    // update basis inverse
    // note: (-1)\hat{\nu} is stored instead of \hat{\nu}
    inv_M_B.z_replace_original_by_original( q_lambda.begin(), q_x_O.begin(), 
        s_delta, -nu, k);

}

// update of the vectors w and r for U_Z_1 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_original_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_Z_1 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                           
// Upper bounded      
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_original_upd_w_r(Tag_false )
{

    ET      x_j = nonbasic_original_variable_value(j);
    ET      x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+i) : *(qp_u+i);
    
    // Note: w needs to be updated before r_C, r_S_B
    update_w_r_B_O__j_i(x_j, x_i);
    update_r_C_r_S_B__j_i(x_j, x_i);
    
    // replace r_beta_O(i) with w_j    
    r_B_O[in_B[i]] = w[j];
    
    // update x_O_v_i
    x_O_v_i[j] = BASIC;
    x_O_v_i[i] = ratio_test_bound_index;    
}


// replacement with precond det(M_{B \setminus \{i\}})=0
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_slack( )
{
    // updates for the upper bounded case
    z_replace_variable_original_by_slack_upd_w_r(Is_nonnegative());
    
    int  k = in_B[ i];

    // leave original variable [ out: i ]
    in_B  [ B_O.back()] = k;
       B_O[ k] = B_O.back();
       in_B  [ i         ] = -1;
       B_O.pop_back();

    minus_c_B[ k] = minus_c_B[ B_O.size()];

    // enter slack variable [ in: j ]
    int  old_row = slack_A[ j-qp_n].first;
    in_B  [ j] = static_cast<int>(B_S.size());
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
    inv_M_B.z_replace_original_by_slack( );

}

// update of the vectors w and r for U_Z_2 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                         // Standard form
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_slack_upd_w_r(Tag_true )
{
}


// update of the vectors w and r for U_Z_2 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                         // Upper bounded
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_slack_upd_w_r(Tag_false )
{

    ET      x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+i) : *(qp_u+i);
    
    // Note: w needs to be updated before r_C, r_S_B
    update_w_r_B_O__i(x_i);
    update_r_C_r_S_B__i(x_i);
    
    int     sigma_j = slack_A[ j-qp_n].first;
    
    // append r_gamma_C(sigma_j) to r_S_B
    r_S_B.push_back(r_C[in_C[sigma_j]]);
    
    // remove r_gamma_C(sigma_j) from r_C
    r_C[in_C[sigma_j]] = r_C.back();
    r_C.pop_back();
    
    // remove r_beta_O(i) from r_B_O    
    if (!check_tag(Is_linear())) { // (kf.)
      r_B_O[in_B[i]] = r_B_O.back();
      r_B_O.pop_back();
    }
    
    // update x_O_v_i
    x_O_v_i[i] = ratio_test_bound_index;
}


// replacement with precond det(M_{B \setminus \{i\}})=0
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_original( )
{
    // updates for the upper bounded case
    z_replace_variable_slack_by_original_upd_w_r(Is_nonnegative());
    
    int  k = in_B[ i];
    if (minus_c_B.size() <= B_O.size()) {  // Note: minus_c_B and the
					   // containers resized in this
					   // if-block are only enlarged
					   // and never made smaller
					   // (while B_O always has the
					   // correct size). We check here
					   // whether we need to enlarge
					   // them.
      CGAL_qpe_assertion(minus_c_B.size() == B_O.size());
	 minus_c_B.push_back(et0);
	     q_x_O.push_back(et0);
	   tmp_x  .push_back(et0);
	   tmp_x_2.push_back(et0);
	  two_D_Bj.push_back(et0);
	     x_B_O.push_back(et0);
    }

    // enter original variable [ in: j ]

    minus_c_B[ B_O.size()] = -ET( *(qp_c+ j));
    

    in_B  [ j] = static_cast<int>(B_O.size());
       B_O.push_back( j);


    // leave slack variable [ out: i ]
       B_S[ k         ] = B_S.back();
       S_B[ k         ] = S_B.back();
    in_B  [ B_S.back()] = k;
    in_B  [ i         ] = -1; 
       B_S.pop_back();
       S_B.pop_back();

    // enter inequality constraint [ in: i ]
    int new_row = slack_A[ i-qp_n].first;

     b_C[ C.size()] = ET( *(qp_b+ new_row));
    in_C[ new_row ] = static_cast<int>(C.size());
       C.push_back( new_row);
    
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }

    // update basis inverse
    // --------------------

    // prepare u
    A_row_by_index_accessor  a_accessor =
      boost::bind (A_accessor( qp_A, 0, qp_n), _1, new_row);
    std::copy( A_row_by_index_iterator( B_O.begin(), a_accessor),
	       A_row_by_index_iterator( B_O.end  (), a_accessor),
	       tmp_x.begin());
    
    
    // prepare kappa
    ET  kappa = d * ET(*((*(qp_A + j))+new_row)) 
                - inv_M_B.inner_product_x( tmp_x.begin(), q_x_O.begin());
		
	// note: (-1)/hat{\nu} is stored instead of \hat{\nu}	 
    inv_M_B.z_replace_slack_by_original( q_lambda.begin(), q_x_O.begin(),
                                          tmp_x.begin(), kappa, -nu);    
}

// update of the vectors w and r for U_Z_3 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_original_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_Z_3 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                             // Upper bounded
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_original_upd_w_r(Tag_false )
{

    ET      x_j = nonbasic_original_variable_value(j);

    // Note: w needs to be updated before r_C, r_S_B    
    update_w_r_B_O__j(x_j);
    update_r_C_r_S_B__j(x_j);
        
    // append r_gamma_S_B(sigma_i) to r_C
    r_C.push_back(r_S_B[in_B[i]]);
    
    // remove r_gamma_S_B(sigma_i) from r_S_B
    r_S_B[in_B[i]] = r_S_B.back();
    r_S_B.pop_back();
    
    // append w_j to r_B_O    
    if (!check_tag(Is_linear())) // (kf.)
      r_B_O.push_back(w[j]);
    
    // update x_O_v_i
    x_O_v_i[j] = BASIC;
}


// replacement with precond det(M_{B \setminus \{i\}})=0
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_slack( )
{
    // updates for the upper bounded case
    z_replace_variable_slack_by_slack_upd_w_r(Is_nonnegative());
    
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

     b_C[ k] = ET( *(qp_b+ new_row));

    // diagnostic output
    CGAL_qpe_debug {
	    if ( vout2.verbose()) print_basis();
    }

    // update basis inverse
    // --------------------
    A_row_by_index_accessor  a_accessor =
      boost::bind ( A_accessor( qp_A, 0, qp_n), _1, new_row);
    std::copy( A_row_by_index_iterator( B_O.begin(), a_accessor),
	       A_row_by_index_iterator( B_O.end  (), a_accessor),
	       tmp_x.begin());


    inv_M_B.z_replace_slack_by_slack( tmp_x.begin(), k);
}

// update of the vectors w and r for U_Z_4 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                             // Standard form
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_slack_upd_w_r(Tag_true )
{
}


// update of the vectors w and r for U_Z_4 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                             // Upper bounded
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_slack_upd_w_r(Tag_false )
{
    
    int     sigma_j = slack_A[ j-qp_n].first;
    
    // swap r_sigma_j in r_C with r_sigma_i in r_S_B
    std::swap(r_C[in_C[sigma_j]], r_S_B[in_B[i]]);           
}

// update of the vectors r_C and r_S_B with "x_j" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_r_C_r_S_B__j(ET& x_j)
{
    // update of vector r_{C}
    A_by_index_accessor     a_accessor_j(*(qp_A + j));
    
    A_by_index_iterator     A_C_j_it(C.begin(), a_accessor_j);
    
    for (Value_iterator r_C_it = r_C.begin(); r_C_it != r_C.end();
                                                    ++r_C_it, ++A_C_j_it) {
        *r_C_it -= x_j * ET(*A_C_j_it); 
    }
    
    // update of r_{S_{B}}
    A_by_index_iterator     A_S_B_j_it(S_B.begin(), a_accessor_j);
        
    for (Value_iterator r_S_B_it = r_S_B.begin(); r_S_B_it != r_S_B.end();
                                                ++r_S_B_it, ++A_S_B_j_it) {
        *r_S_B_it -= x_j * ET(*A_S_B_j_it);
    }   
}

// update of the vectors r_C and r_S_B with "x_j" and "x_i" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_r_C_r_S_B__j_i(ET& x_j, ET& x_i)
{
    // update of vector r_{C}
    A_by_index_accessor     a_accessor_j(*(qp_A+j));
    A_by_index_accessor     a_accessor_i(*(qp_A+i));
    
    A_by_index_iterator     A_C_j_it(C.begin(), a_accessor_j);
    A_by_index_iterator     A_C_i_it(C.begin(), a_accessor_i);
    
    for (Value_iterator r_C_it = r_C.begin(); r_C_it != r_C.end();
                                        ++r_C_it, ++A_C_j_it, ++A_C_i_it) {
        *r_C_it += (x_i * ET(*A_C_i_it)) - (x_j * ET(*A_C_j_it)); 
    }
    
    // update of vector r_{S_{B}}
    A_by_index_iterator     A_S_B_j_it(S_B.begin(), a_accessor_j);
    A_by_index_iterator     A_S_B_i_it(S_B.begin(), a_accessor_i);

    for (Value_iterator r_S_B_it = r_S_B.begin(); r_S_B_it != r_S_B.end();
                                    ++r_S_B_it, ++A_S_B_j_it, ++A_S_B_i_it) {
        *r_S_B_it += (x_i * ET(*A_S_B_i_it)) - (x_j * ET(*A_S_B_j_it)); 
    }
}

// update of the vectors r_C and r_S_B with "x_i'" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_r_C_r_S_B__i(ET& x_i)
{
    // update of vector r_{C}
    A_by_index_accessor     a_accessor_i(*(qp_A+i));
    A_by_index_iterator     A_C_i_it(C.begin(), a_accessor_i);
    
    for (Value_iterator r_C_it = r_C.begin(); r_C_it != r_C.end(); 
                                                    ++r_C_it, ++A_C_i_it) {
        *r_C_it += x_i * ET(*A_C_i_it); 
    }
    
    // update of r_{S_{B}}
    A_by_index_iterator     A_S_B_i_it(S_B.begin(), a_accessor_i);
    
    for (Value_iterator r_S_B_it = r_S_B.begin(); r_S_B_it != r_S_B.end();
                                                ++r_S_B_it, ++A_S_B_i_it) {
        *r_S_B_it += x_i * ET(*A_S_B_i_it);
    }   
}


// Update of w and r_B_O with "x_j" column.
//
// todo: could be optimized slightly by factoring out the factor 2
// (which is implicitly contained in the pairwise accessor for D).
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_w_r_B_O__j(ET& x_j)
{
  // assertion checking:
  CGAL_expensive_assertion(!check_tag(Is_nonnegative()));

  // Note: we only do anything it we are dealing with a QP.
  if (!check_tag(Is_linear())) {
    D_pairwise_accessor     d_accessor_j(qp_D, j);
    
    // update of vector w:
    for (int it = 0; it < qp_n; ++it)
      w[it] -= d_accessor_j(it) * x_j;
    
    // update of r_B_O:
    D_pairwise_iterator D_B_O_j_it(B_O.begin(), d_accessor_j);
    for (Value_iterator r_B_O_it = r_B_O.begin();
	 r_B_O_it != r_B_O.end();
	 ++r_B_O_it, ++D_B_O_j_it)
        *r_B_O_it -= *D_B_O_j_it * x_j;
  }
}

// update of w and r_B_O with "x_j" and "x_i" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_w_r_B_O__j_i(ET& x_j, ET& x_i)
{
  // assertion checking:
  CGAL_expensive_assertion(!check_tag(Is_nonnegative()));

  // Note: we only do anything it we are dealing with a QP.
  if (!check_tag(Is_linear())) {
    D_pairwise_accessor     d_accessor_i(qp_D, i);
    D_pairwise_accessor     d_accessor_j(qp_D, j);
    
    // update of vector w
    for (int it = 0; it < qp_n; ++it)
        w[it] += (d_accessor_i(it) * x_i) - (d_accessor_j(it) * x_j);
    
    // update of r_B_O
    D_pairwise_iterator D_B_O_j_it(B_O.begin(), d_accessor_j);
    D_pairwise_iterator D_B_O_i_it(B_O.begin(), d_accessor_i);
    for (Value_iterator r_B_O_it = r_B_O.begin();
	 r_B_O_it != r_B_O.end();
	 ++r_B_O_it, ++D_B_O_j_it, ++D_B_O_i_it)
        *r_B_O_it += (*D_B_O_i_it * x_i) - (*D_B_O_j_it * x_j);
  }
}

// update of w and r_B_O with "x_i" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_w_r_B_O__i(ET& x_i)
{
  CGAL_expensive_assertion(!check_tag(Is_nonnegative()));

  // Note: we only do anything it we are dealing with a QP.
  if (!check_tag(Is_linear())) {

    // update of vector w:
    D_pairwise_accessor     d_accessor_i(qp_D, i);
    
    for (int it = 0; it < qp_n; ++it)
        w[it] += d_accessor_i(it) * x_i;
    
    // update of r_B_O:
    D_pairwise_iterator D_B_O_i_it(B_O.begin(), d_accessor_i);
    for (Value_iterator r_B_O_it = r_B_O.begin();
	 r_B_O_it != r_B_O.end();
	 ++r_B_O_it, ++D_B_O_i_it)
        *r_B_O_it += *D_B_O_i_it * x_i;
  }
}

// Compute solution, meaning compute the solution vector x and the KKT
// coefficients lambda.
template < typename Q, typename ET, typename Tags >                                      // Standard form
void  QP_solver<Q, ET, Tags>::
compute_solution(Tag_true)
{
  // compute current solution, original variables and lambdas
  inv_M_B.multiply( b_C.begin(), minus_c_B.begin(),
		    lambda.begin(), x_B_O.begin());
  
  // compute current solution, slack variables
  compute__x_B_S(no_ineq, Is_nonnegative());
}

// Compute solution, meaning compute the solution vector x and the KKT
// coefficients lambda.
template < typename Q, typename ET, typename Tags >                                      // Upper bounded
void  QP_solver<Q, ET, Tags>::
compute_solution(Tag_false)
{ 
  // compute the difference b_C - r_C
  std::transform(b_C.begin(), b_C.begin()+C.size(), // Note: r_C.size() ==
						    // C.size() always holds,
						    // whereas b_C.size() >=
						    // C.size() in general.
		 r_C.begin(),  tmp_l.begin(), std::minus<ET>());
        
  // compute the difference minus_c_B - r_B_O:
  if (is_phaseII && is_QP) {
    std::transform(minus_c_B.begin(), minus_c_B.begin() + B_O.size(), 
		   r_B_O.begin(), tmp_x.begin(), std::minus<ET>());
    
    // compute current solution, original variables and lambdas:
    inv_M_B.multiply( tmp_l.begin(), tmp_x.begin(),
		      lambda.begin(), x_B_O.begin());
  } else {                                          // r_B_O == 0
    
    // compute current solution, original variables and lambdas        
    inv_M_B.multiply( tmp_l.begin(), minus_c_B.begin(),
		      lambda.begin(), x_B_O.begin());
  }
  
  // compute current solution, slack variables
  compute__x_B_S( no_ineq, Is_nonnegative());
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__A_S_BxB_O(Value_iterator in, Value_iterator out) const
{
  // initialize with zero vector
  std::fill_n( out, B_S.size(), et0);
  
  // foreach original column of A in B_O (artificial columns are zero in S_B
  A_column              a_col;                             // except special)
  Index_const_iterator  row_it, col_it;
  Value_iterator        out_it;
  //ET                    in_value;
  for ( col_it = B_O.begin(); col_it != B_O.end(); ++col_it, ++in) {
    const ET in_value = *in;
    out_it   = out;
    
    if ( *col_it < qp_n) {	                        // original variable
      a_col = *(qp_A+ *col_it);
      
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

// compare the updated vector r_{C} with t_r_C=A_{C, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                         // Standard form
bool  QP_solver<Q, ET, Tags>::
check_r_C(Tag_true) const
{
    return true;
}

// compare the updated vector r_{C} with t_r_C=A_{C, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                         // Upper bounded
bool  QP_solver<Q, ET, Tags>::
check_r_C(Tag_false) const
{
    Values                  t_r_C;
    // compute t_r_C from scratch
    t_r_C.resize(C.size(), et0);
    multiply__A_CxN_O(t_r_C.begin());
    
    // compare r_C and the from scratch computed t_r_C
    bool failed = false;
    int count = 0;
    Value_const_iterator t_r_C_it = r_C.begin();
    for (Value_const_iterator r_C_it = r_C.begin(); r_C_it != r_C.end();
                                    ++r_C_it, ++t_r_C_it, ++count) {
        if (*r_C_it != *t_r_C_it) {
            failed = true;
        }
    }
    return (!failed);
}

// compare the updated vector r_{S_B} with t_r_S_B=A_{S_B, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                             // Standard form
bool  QP_solver<Q, ET, Tags>::
check_r_S_B(Tag_true) const
{
    return true;
}

// compare the updated vector r_{S_B} with t_r_S_B=A_{S_B, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                             // Upper bounded
bool  QP_solver<Q, ET, Tags>::
check_r_S_B(Tag_false) const
{
    Values                  t_r_S_B;
    // compute t_r_S_B from scratch
    t_r_S_B.resize(S_B.size(), et0);
    multiply__A_S_BxN_O(t_r_S_B.begin());
    
    // compare r_S_B and the from scratch computed t_r_S_B
    bool failed = false;
    int count = 0;
    Value_const_iterator    t_r_S_B_it = t_r_S_B.begin();
    for (Value_const_iterator r_S_B_it = r_S_B.begin(); r_S_B_it != r_S_B.end();
                                    ++r_S_B_it, ++t_r_S_B_it, ++count) {
        if (*r_S_B_it != *t_r_S_B_it) {
            failed = true;
        }
    }
    return (!failed);
}


// computes r_{B_{O}}:=2D_{B_O, N_O}x_{N_O} with upper bounding
// OPTIMIZATION: If D is symmetric we can multiply by two at the end of the
// computation of entry of r_B_O instead of each access to D
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__2D_B_OxN_O(Value_iterator out) const
{
    //initialize
    std::fill_n( out, B_O.size(), et0);

    Value_iterator          out_it;
    ET                      value;
    
    // foreach entry in r_B_O
    out_it = out;
    for (Index_const_iterator row_it = B_O.begin(); row_it != B_O.end();
                                                        ++row_it, ++out_it) {
        D_pairwise_accessor d_accessor_i(qp_D, *row_it);
        for (int i = 0; i < qp_n; ++i) {
            if (!is_basic(i)) {
                value = nonbasic_original_variable_value(i);
                *out_it += d_accessor_i(i) * value;
            }
        }    
    }
}

// compares the updated vector r_{B_O} with t_r_B_O=2D_{B_O, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                                // Standard form
bool  QP_solver<Q, ET, Tags>::
check_r_B_O(Tag_true) const
{
    return true;
}

// compares the updated vector r_{B_O} with t_r_B_O=2D_{B_O, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                                 // Upper bounded
bool  QP_solver<Q, ET, Tags>::
check_r_B_O(Tag_false) const
{
    Values                  t_r_B_O;
    // compute t_r_B_O from scratch
    t_r_B_O.resize(B_O.size(), et0);
    multiply__2D_B_OxN_O(t_r_B_O.begin());
    
    // compare r_B_O and the from scratch computed t_r_B_O
    bool failed = false;
    int count = 0;
    Value_const_iterator    t_r_B_O_it = t_r_B_O.begin();
    for (Value_const_iterator r_B_O_it = r_B_O.begin(); r_B_O_it != r_B_O.end();
                                    ++r_B_O_it, ++t_r_B_O_it, ++count) {
        if (*r_B_O_it != *t_r_B_O_it) {
            failed = true;
        }
    }
    return (!failed);   
}

// compares the updated vector w with t_w=2D_{O,N_O}*x_N_O
template < typename Q, typename ET, typename Tags >                             // Standard form
bool  QP_solver<Q, ET, Tags>::
check_w(Tag_true) const
{
    return true;
}

// compares the updated vector w with t_w=2D_O_N_O*x_N_O
template < typename Q, typename ET, typename Tags >                             // Upper bounded
bool  QP_solver<Q, ET, Tags>::
check_w(Tag_false) const
{
    Values              t_w;
    // compute t_w from scratch
    t_w.resize( qp_n, et0);
    multiply__2D_OxN_O(t_w.begin());
    
    // compare w and the from scratch computed t_w
    bool  failed = false;
    Value_const_iterator    t_w_it = t_w.begin();
    for (int i = 0; i < qp_n; ++i, ++t_w_it) {
        if (w[i] != *t_w_it) {
            failed = true;
        }
    }
    return (!failed);
}

// check basis inverse
template < typename Q, typename ET, typename Tags >
bool
QP_solver<Q, ET, Tags>::
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

template < typename Q, typename ET, typename Tags >                                         // LP case
bool  QP_solver<Q, ET, Tags>::
check_basis_inverse( Tag_true)
{
    CGAL_qpe_debug {
	vout4 << std::endl;
    }
    bool res = true;
    unsigned int    row, rows =   static_cast<unsigned int>(C.size());
    unsigned int    col, cols = static_cast<unsigned int>(B_O.size());
    Index_iterator  i_it = B_O.begin();
    Value_iterator  q_it;

    
    // BG: is this a real check?? How does the special artifical
    // variable come in, e.g.? OK: it comes in through
    // ratio_test_init__A_Cj
    for ( col = 0; col < cols; ++col, ++i_it) {
	ratio_test_init__A_Cj( tmp_l.begin(), *i_it,
			       no_ineq);
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

template < typename Q, typename ET, typename Tags >                                         // QP case
bool  QP_solver<Q, ET, Tags>::
check_basis_inverse( Tag_false)
{
    bool res = true;
    unsigned int    row, rows =   static_cast<unsigned int>(C.size());
    unsigned int    col, cols = static_cast<unsigned int>(B_O.size());
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
	    *v_it = ( *i_it < qp_n ? *((*(qp_A+ *i_it))+ row) :  // original
		      art_A[ *i_it - qp_n].first != static_cast<int>(row) ? et0 :// artific.
		      ( art_A[ *i_it - qp_n].second ? -et1 : et1));
	}
// 	if ( art_s_i >= 0) {              // special artificial variable?
// 	  CGAL_qpe_assertion (static_cast<int>(in_B.size()) == art_s_i+1);  
// 	  // the special artificial variable has never been
// 	  // removed from the basis, consider it
// 	  tmp_x[ in_B[ art_s_i]] = art_s[ row];
// 	} 
// BG: currently, this check only runs in phase II, where we have no
// articficials
	CGAL_qpe_assertion (art_s_i < 0);
	
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
			     boost::bind2nd( std::not_equal_to<ET>(), et0));
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
	ratio_test_init__A_Cj  ( tmp_l.begin(), *i_it, 
				 no_ineq);
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
			     boost::bind2nd( std::not_equal_to<ET>(), et0));
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

// filtered strategies are only allowed if the input type as
// indicated by C_entry is double; this may still fail if the
// types of some other program entries are not double, but 
// since the filtered stuff is internal anyway, we can probably
// live with this simple check
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
set_pricing_strategy
( Quadratic_program_pricing_strategy strategy)
{
    CGAL_qpe_assertion( phase() != 1);
    CGAL_qpe_assertion( phase() != 2);

    if (strategy == QP_DANTZIG)
      strategyP = new QP_full_exact_pricing<Q, ET, Tags>;
    else if (strategy == QP_FILTERED_DANTZIG)    
      // choose between FF (double) and FE (anything else)
      strategyP = 
	new typename QP_solver_impl::Filtered_pricing_strategy_selector
	<Q, ET, Tags, C_entry>::FF;
    else if (strategy == QP_PARTIAL_DANTZIG)
      strategyP = new QP_partial_exact_pricing<Q, ET, Tags>;
    else if (strategy == QP_PARTIAL_FILTERED_DANTZIG
	     || strategy == QP_CHOOSE_DEFAULT)   
      // choose between PF (double) and PE (anything else)
      strategyP = 
	new typename QP_solver_impl::Filtered_pricing_strategy_selector
	<Q, ET, Tags, C_entry>::PF;
    else if (strategy == QP_BLAND) 
      strategyP = new QP_exact_bland_pricing<Q, ET, Tags>;
    CGAL_qpe_assertion(strategyP != static_cast<Pricing_strategy*>(0));
   
    if ( phase() != -1) strategyP->set( *this, vout2);
}

// diagnostic output
// -----------------
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
set_verbosity( int verbose, std::ostream& stream)
{
    vout  = Verbose_ostream( verbose >  0, stream);
    vout1 = Verbose_ostream( verbose == 1, stream);
    vout2 = Verbose_ostream( verbose >= 2, stream);
    vout3 = Verbose_ostream( verbose >= 3, stream);
    vout4 = Verbose_ostream( verbose == 4, stream);
    vout5 = Verbose_ostream( verbose == 5, stream);
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_program( ) const
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
    vout4.out() << "(+ " << qp_c0 << ") ";
    vout4.out() << std::endl;
    vout4.out() << std::endl;

    // constraints
    vout4.out() << "constraints:" << std::endl;
    for ( row = 0; row < qp_m; ++row) {

	// original variables
	for (i = 0; i < qp_n; ++i) 
	  vout4.out() << *((*(qp_A+ i))+row) << ' ';

	// slack variables
	if ( ! slack_A.empty()) {
	    vout4.out() << " |  ";
	    for ( i = 0; i < static_cast<int>(slack_A.size()); ++i) {
		vout4.out() << ( slack_A[ i].first != row ? " 0" :
		               ( slack_A[ i].second ? "-1" : "+1")) << ' ';
	    }
	}

	// artificial variables
	if ( ! art_A.empty()) {
	    vout4.out() << " |  ";
	    for ( i = 0; i < static_cast<int>(art_A.size()); ++i) {
	      if (art_s_i == i+qp_n+static_cast<int>(slack_A.size()))
		vout4.out() << " * ";          // for special artificial column
	      vout4.out() << ( art_A[ i].first != row ? " 0" :
		             ( art_A[ i].second ? "-1" : "+1")) << ' ';
	    }
	}
	if ( ! art_s.empty()) vout4.out() << " |  " << art_s[ row] << ' ';

	// rhs
	vout4.out() << " |  "
		    << ( *(qp_r+ row) == CGAL::EQUAL      ? ' ' :
		       ( *(qp_r+ row) == CGAL::SMALLER ? '<' : '>')) << "=  "
		    << *(qp_b+ row);
		    if (!is_nonnegative) {
		        vout4.out() << " - " << multiply__A_ixO(row);
		    }
		    vout4.out() << std::endl;
    }
    vout4.out() << std::endl;
    
    // explicit bounds
    if (!is_nonnegative) {
        vout4.out() << "explicit bounds:" << std::endl; 
        for (int i = 0; i < qp_n; ++i) {
            if (*(qp_fl+i)) {                   // finite lower bound
                vout4.out() << *(qp_l+i);
            } else {                            // infinite lower bound
                vout4.out() << "-inf";
            }
            vout4.out() << " <= x_" << i << " <= ";
            if (*(qp_fu+i)) {
                vout4.out() << *(qp_u+i);
            } else {
                vout4.out() << "inf";
            }
            vout4.out() << std::endl;
        }
    }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_basis( ) const
{
    char label;
    vout << "  basis: ";
    CGAL_qpe_debug {
      vout2 << "basic variables" << ( has_ineq ? "  " : "") << ":  ";
    }
    std::copy( B_O.begin(), B_O.end  (),
	       std::ostream_iterator<int>( vout.out(), " "));
    CGAL_qpe_debug {
      if ( vout2.verbose()) {
	if ( has_ineq && ( ! slack_A.empty())) {
	  vout2.out() << " |  ";
	  std::copy( B_S.begin(), B_S.end(),
		     std::ostream_iterator<int>( vout2.out(), " "));
	}
	if ( is_phaseI) {
	  vout2.out() << " (# of artificials: " << art_basic << ')';
	}
	if ( has_ineq) {
	  vout2.out() << std::endl
		      << "basic constraints:  ";
	  for (Index_const_iterator i_it = 
		 C.begin(); i_it != C.end(); ++i_it) {
	    label = (*(qp_r+ *i_it) == CGAL::EQUAL) ? 'e' : 'i';
	    vout2.out() << *i_it << ":" << label << " ";
	  }
	  /*
	    std::copy( C.begin(), C.begin()+(C.size()-slack_A.size()),
	    std::ostream_iterator<int>( vout2.out(), " "));
	    if ( ! slack_A.empty()) {
	    vout2.out() << " |  ";
	    std::copy( C.end() - slack_A.size(), C.end(),
	    std::ostream_iterator<int>( vout2.out(), " "));
	    }
	  */
	}
	if ( vout3.verbose()) {
	  vout3.out() << std::endl
		      << std::endl
		      << "    in_B: ";
	  std::copy( in_B.begin(), in_B.end(),
		     std::ostream_iterator<int>( vout3.out(), " "));
	  if ( has_ineq) {
	    vout3.out() << std::endl
			<< "    in_C: ";
	    std::copy( in_C.begin(), in_C.end(),
		       std::ostream_iterator<int>( vout3.out(), " "));
	  }
	}
      }
    }
    vout.out() << std::endl;
    CGAL_qpe_debug {
      vout4 << std::endl << "basis-inverse:" << std::endl;
    }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_solution( ) const
{
  CGAL_qpe_debug {
    if ( vout3.verbose()) {
	   vout3.out() << std::endl
            << "     b_C: ";
	   std::copy( b_C.begin(), b_C.begin()+C.size(),
		   std::ostream_iterator<ET>( vout3.out()," "));
	   vout3.out() << std::endl
            << "  -c_B_O: ";
	   std::copy( minus_c_B.begin(), minus_c_B.begin()+B_O.size(),
		   std::ostream_iterator<ET>( vout3.out()," "));
	   if (!is_nonnegative) {
	       vout3.out() << std::endl
                << "     r_C: ";
	       std::copy( r_C.begin(), r_C.begin()+r_C.size(),
	           std::ostream_iterator<ET>( vout3.out(), " "));
	       vout3.out() << std::endl
	           << "   r_B_O: ";
	       std::copy( r_B_O.begin(), r_B_O.begin()+r_B_O.size(),
	           std::ostream_iterator<ET>( vout3.out(), " "));
	       if (r_B_O.size() == 0)
		 vout3.out() << "< will only be allocated in phase II>";

	   }
	   vout3.out() << std::endl;
    }
    if ( vout2.verbose()) {
        vout2.out() << std::endl << "  lambda: ";
        std::copy( lambda.begin(), lambda.begin()+C.size(),
            std::ostream_iterator<ET>( vout2.out(), " "));
        vout2.out() << std::endl << "   x_B_O: ";
        std::copy( x_B_O.begin(), x_B_O.begin()+B_O.size(),
            std::ostream_iterator<ET>( vout2.out(), " "));
        vout2.out() << std::endl;
        if (!is_nonnegative) {
            vout2.out() << "   x_N_O: ";
            for (int i = 0; i < qp_n; ++i) {
                if (!is_basic(i)) {
                    vout2.out() << d * 
		      nonbasic_original_variable_value(i);
                    vout2.out() << " ";    
                }
            }
            vout2.out() << std::endl;
        }
        if ( has_ineq) {
            vout2.out() << "   x_B_S: ";
            std::copy( x_B_S.begin(), x_B_S.begin()+B_S.size(),
                std::ostream_iterator<ET>( vout2.out()," "));
            vout2.out() << std::endl;
        }
        const ET denom = inv_M_B.denominator();
        vout2.out() << "   denominator: " << denom << std::endl;
        vout2.out() << std::endl;
    }
  }
  vout << "  ";
  vout.out() 
    << "solution: " 
    << solution_numerator() << " / " << solution_denominator() 
    << "  ~= " 
    << to_double
    (CGAL::Quotient<ET>(solution_numerator(), solution_denominator()))
    << std::endl;
  CGAL_qpe_debug {
      vout2 << std::endl;
  }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_ratio_1_original(int k, const ET& x_k, const ET& q_k)
{
    if (is_nonnegative) {                      // => direction== 1
        if (q_k > et0) {                            // check for lower bound
            vout2.out() << "t_O_" << k << ": "
            << x_k << '/' << q_k
            << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
            << std::endl;
        } else if (q_k < et0) {                     // check for upper bound
            vout2.out() << "t_O_" << k << ": "
            << "inf" << '/' << q_k
            << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
            << std::endl;
        } else {                                    // q_k == 0
            vout2.out() << "t_O_" << k << ": "
            << "inf"
            << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
            << std::endl;
        }
    } else {                                        // upper bounded
        if (q_k * ET(direction) > et0) {            // check for lower bound
            if (B_O[k] < qp_n) {                         // original variable
                if (*(qp_fl+B_O[k])) {                   // finite lower bound
                    vout2.out() << "t_O_" << k << ": "
                    << x_k - (d * ET(*(qp_l+B_O[k]))) << '/' << q_k
                    << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
                    << std::endl;
                } else {                            // lower bound -infinity
                    vout2.out() << "t_O_" << k << ": "
                    << "-inf" << '/' << q_k
                    << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
                    << std::endl;                
                }
            } else {                                // artificial variable
                vout2.out() << "t_O_" << k << ": "
                << x_k << '/' << q_k
                << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
                << std::endl;
            }
        } else if (q_k * ET(direction) < et0) {     // check for upper bound
            if (B_O[k] < qp_n) {                         // original variable
                if (*(qp_fu+B_O[k])) {                   // finite upper bound
                    vout2.out() << "t_O_" << k << ": "
                    << (d * ET(*(qp_l+B_O[k]))) - x_k << '/' << q_k
                    << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
                    << std::endl;                    
                } else {                            // upper bound infinity
                    vout2.out() << "t_O_" << k << ": "
                    << "inf" << '/' << q_k
                    << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
                    << std::endl;
                }
            } else {                                // artificial variable
                vout2.out() << "t_O_" << k << ": "
                << "inf" << '/' << q_k
                << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
                << std::endl;
            }
        } else {                                    // q_k == 0
            vout2.out() << "t_O_" << k << ": "
            <<  "inf"
            << ( ( q_i != et0) && ( i == B_O[ k]) ? " *" : "")
            << std::endl;
        }
    }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_ratio_1_slack(int k, const ET& x_k, const ET& q_k)
{
    if (is_nonnegative) {                      // => direction == 1
        if (q_k > et0) {                            // check for lower bound
            vout2.out() << "t_S_" << k << ": "
            << x_k << '/' << q_k
            << ( ( q_i != et0) && ( i == B_S[ k]) ? " *" : "")
            << std::endl;
        } else {                                    // check for upper bound
            vout2.out() << "t_S_" << k << ": "
            << "inf" << '/' << q_k
            << ( ( q_i != et0) && ( i == B_S[ k]) ? " *" : "")
            << std::endl;        
        }
    } else {                                        // upper bounded
        if (q_k * ET(direction) > et0) {            // check for lower bound
            vout2.out() << "t_S_" << k << ": "
            << x_k << '/' << q_k
            << ( ( q_i != et0) && ( i == B_S[ k]) ? " *" : "")
            << std::endl;            
        } else if (q_k * ET(direction) < et0) {     // check for upper bound
            vout2.out() << "t_S_" << k << ": "
            << "inf" << '/' << q_k
            << ( ( q_i != et0) && ( i == B_S[ k]) ? " *" : "")
            << std::endl;
        } else {                                    // q_k == 0
            vout2.out() << "t_S_" << k << ": "
            << "inf"
            << ( ( q_i != et0) && ( i == B_S[ k]) ? " *" : "")
            << std::endl;
        }
    }
}


template < typename Q, typename ET, typename Tags >
const char*  QP_solver<Q, ET, Tags>::
variable_type( int k) const
{
    return ( k <        qp_n                 ? "original"  :
	   ( k < static_cast<int>( qp_n+slack_A.size()) ? "slack"     :
	                                       "artificial"));
}

template < typename Q, typename ET, typename Tags > 
bool QP_solver<Q, ET, Tags>::
is_artificial(int k) const
{
    return (k >= static_cast<int>(qp_n+slack_A.size())); 
}

template < typename Q, typename ET, typename Tags > 
int QP_solver<Q, ET, Tags>::
get_l() const
{
    return l;
}

} //namespace CGAL

// ===== EOF ==================================================================
