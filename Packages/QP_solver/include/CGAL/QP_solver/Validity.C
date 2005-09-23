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
// file          : include/CGAL/QP_solver/Validity.C
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Validity checks
// ============================================================================

CGAL_BEGIN_NAMESPACE

template < class Rep_ >
bool QP_solver<Rep_>::is_valid()
{
    switch(this->m_status) {
    case UPDATE:
      return false;
    case OPTIMAL:
      {
	const bool f = this->is_solution_feasible();
	const bool o = this->is_solution_optimal(Is_linear());
	CGAL_qpe_debug {
	  vout << "feasible: " << f << std::endl;
	  vout << "optimal: " << o << std::endl;
	}
	return f && o;
      }
    case INFEASIBLE:
      {
	const bool f = this->is_solution_feasible_for_auxiliary_problem();
	const bool o = this->is_solution_feasible_for_auxiliary_problem();
	const bool aux_positive = this->solution() > et0;
	CGAL_qpe_debug {
	  vout << "feasible_aux: " << f << std::endl;
	  vout << "optimal_aux: " << o << std::endl;
	  vout << "objective value positive: " << aux_positive << std::endl;
	}
	return f && o && aux_positive;
      }
    case UNBOUNDED:    
      {
	const bool f = this->is_solution_feasible();
	const bool u = this->is_solution_unbounded(Is_linear());
	CGAL_qpe_debug {
	  vout << "feasible: " << f << std::endl;
	  vout << "unbounded: " << u << std::endl;
	}
	return f && u;
      }
    default: 	      
      return false;
    }
}

template < class Rep_ >
bool QP_solver<Rep_>::is_solution_feasible_for_auxiliary_problem()
{
  // some simple consistency checks:
  CGAL_qpe_assertion(is_phaseI);

  // check nonnegativity of original- and artificial variables:
  for (Value_const_iterator v_it = x_B_O.begin(); v_it != x_B_O.end(); ++v_it)
    if (*v_it < et0)
      return false;
  
  // check nonegativity of slack variables:
  for (Value_const_iterator v_it = x_B_S.begin(); v_it != x_B_S.end(); ++v_it)
    if (*v_it < et0)
      return false;
  
  // Check whether the current solution is feasible for the auxiliary
  // problem.  For this, we use the fact that the auxiliary problem
  // looks as follows:
  //
  //                          [  x_original  ]
  //    [ A A_art a_spec A_s] [ x_artificial ] = b               (C1)
  //                          [   x_special  ]
  //                          [    x_slack   ]
  //
  // Here, every column i of A_art corresponds to an equality
  // constraint of the original problem and contains either e_i or
  // -e_i, a_spec is the special artificial column (which only
  // contains zeros, ones, or minus ones), and A_s contains a column
  // (e_i or -e_i) for every slack variable.  Observe in the code
  // below that the right-hand size is multiplied by d because we
  // maintain the denominator of the solution (which is d) separately.

  // compute left-hand side of (C1) (up to slackies):
  Values lhs_col(qp_m, et0);
  Value_const_iterator x_it = x_B_O.begin();
  for (Index_const_iterator i_it = B_O.begin();
       i_it != B_O.end();
       ++i_it, ++x_it)                    // iterate over all nonzero vars
    if (*i_it < qp_n)                     // ordinary original variable?
      for (int i=0; i<qp_m; ++i)
	lhs_col[i] += (*x_it) * qp_A[*i_it][i];
    else                                  // artificial variable?
      if (*i_it != art_s_i)   {           // normal artificial variable?
	const int k = *i_it - qp_n - slack_A.size();
	lhs_col[art_A[k].first] += (art_A[ k].second ? -et1 :  et1) * (*x_it);
      } else                              // special artificial variable
	for (int i=0; i<qp_m; ++i)
	  lhs_col[i] += (*x_it) * art_s[i];

  // compute left-hand side of (C1) (part for slackies):
  x_it = x_B_S.begin();
  for (Index_const_iterator i_it = B_S.begin();
       i_it != B_S.end();
       ++i_it, ++x_it) {
    const int k = *i_it - qp_n;
    lhs_col[slack_A[k].first] += (slack_A[ k].second ? -et1 : et1) * (*x_it);
  }
  
  // check equality (C1);
  for (int i=0; i<qp_m; ++i)
    if (lhs_col[i] != ET(qp_b[i]) * d)
      return false;
  
  return true;
}

template < class Rep_ >
bool QP_solver<Rep_>::is_solution_optimal_for_auxiliary_problem()
{
  CGAL_expensive_precondition(is_solution_feasible_for_auxiliary_problem());

  // As described in equation (C1) in routine is_solution_feasible_-
  // for_auxiliary_problem(), the auxiliary problem looks as follows:
  //
  //                   minimize     aux_c^T x
  //                   subject to   aux_A x =  b,
  //                                      x >= 0,
  //
  // where aux_A = [ A A_art a_spec A_s].
  //
  // Using the KKT conditions, we obtain that a feasible point x is
  // optimal iff there exists a m-vector \lambda and a |x|-vector \tau
  // such that
  //
  //              \tau^T = c^T + \lambda^T aux_A,                (C2)
  //              \tau^T x =  0,
  //              \tau     >= 0.
  //
  // These are the conditions we are going to check below. (Notice
  // here that x is the vector containing the original variables, the
  // artificial variables, the special artificial variable, and the
  // slack variables.)

  // get number of working variables:
  // todo: do not understand this (kf): why is art_s_i not a working variable?
  const int no_of_wo_vars = this->number_of_working_variables() +
    (art_s_i == -2)? 1 : 0;    // Note: if there ever was a special
			       // artifial variable, it has to be
			       // considered for this optimality test,
			       // even if it has already left the
			       // basis.
  
  // collect solution vector of auxiliary problem:
  Values x_aux(no_of_wo_vars, et0);
  Value_const_iterator v_it = x_B_O.begin();
  for (Index_const_iterator i_it = B_O.begin();
       i_it != B_O.end(); ++i_it, ++v_it)
    x_aux[*i_it] = *v_it;
  v_it = x_B_S.begin();
  for (Index_const_iterator i_it = B_S.begin();
       i_it != B_S.end(); ++i_it, ++v_it)
    x_aux[*i_it] = *v_it;
  
  // Note: lambda[i] <= 0 for qp_r[i] == "GREATER_EQUAL"
  // todo: (ask frans) what does the above note mean?
  Values lambda_aux(qp_m, et0);
  v_it = lambda.begin();
  for (Index_const_iterator i_it = C.begin();
       i_it != C.end(); ++i_it, ++v_it)
    lambda_aux[*i_it] = *v_it;

  // output for debugging:
  CGAL_qpe_debug {
    for (int col = 0; col < qp_m; ++col)
      vout5 << "lambda_aux[" << col << "]= " << lambda_aux[col] << std::endl;
  }
  
  // compute \tau^T = c^T + \lambda^T * aux_A (see conditions (C2) above):
  //
  // Note: as the \lambda we have access to is actual d times the lambda
  // from (C2), we will not compute \tau but \tau * d.
  //
  // (a) compute \tau^T = c^T:
  Values tau_aux(no_of_wo_vars, et0);
  C_auxiliary_iterator c_it = aux_c.begin();
  for (int col = qp_n+slack_A.size(); col < no_of_wo_vars; ++c_it, ++col)
    tau_aux[col] = ET(*c_it) * d;
  
  // (b) compute \tau^T = c^T + \lambda^T * aux_A:
  for (int col = 0; col < no_of_wo_vars; ++col) {
    if (col < qp_n)                  // ordinary original variable
      for (int i=0; i<qp_m; ++i)
	tau_aux[col] += lambda_aux[i] * qp_A[col][i];
    else {
      int k = col - qp_n;
      if (k < slack_A.size())        // slack variable
	tau_aux[col] = (slack_A[k].second? -et1 :  et1)
	  * lambda_aux[slack_A[k].first];
      else {                         // artificial variable
	k -= slack_A.size();
	if (col != art_s_i)          // normal artificial variable
	  tau_aux[col] += (art_A[ k].second? -et1 : et1)
	    * lambda_aux[art_A[k].first];
	else                         // special artificial variable
	  for (int i=0; i<qp_m; ++i)
	    tau_aux[col] += lambda_aux[i] * art_s[i];
      }
    }
    CGAL_qpe_debug {
      vout5 << "tau_aux[" << col << "]= " << tau_aux[col] << std::endl;
    }
  }

  // check last two lines of (C2):
  for (int col = 0; col < no_of_wo_vars; ++col) {
    if (tau_aux[col] >= et0)
      if (tau_aux[col] * x_aux[col] != et0)
	return false;
    return false;
  }

  return true;     
}

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
