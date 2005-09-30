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
bool QP_solver<Rep_>::has_finite_lower_bound(int i) const
  // Given an index of an original or slack variable, returns whether
  // or not the variable has a finite lower bound.
{
  CGAL_qpe_assertion(i < qp_n + static_cast<int>(slack_A.size()));
  return i>=qp_n || check_tag(Is_in_standard_form()) || *(qp_fl+i);
}

template < class Rep_ >
bool QP_solver<Rep_>::has_finite_upper_bound(int i) const
  // Given an index of an original or slack variable, returns whether
  // or not the variable has a finite upper bound.
{
  CGAL_qpe_assertion(i < qp_n + static_cast<int>(slack_A.size()));
  return i<qp_n && !check_tag(Is_in_standard_form()) && *(qp_fu+i);
}

template < class Rep_ >
typename QP_solver<Rep_>::ET QP_solver<Rep_>::lower_bound(int i) const
  // Given an index of an original or slack variable, returns its
  // lower bound.
{
  CGAL_qpe_assertion(i < qp_n + static_cast<int>(slack_A.size()));
  if (i < qp_n)                     // original variable?
    if (check_tag(Is_in_standard_form()))
      return et0;
    else {
      CGAL_qpe_assertion(has_finite_lower_bound(i));
      return *(qp_l+i);
    }
  else                              // slack variable?
    return et0;
}

template < class Rep_ >
typename QP_solver<Rep_>::ET QP_solver<Rep_>::upper_bound(int i) const
  // Given an index of an original or slack variable, returns its
  // upper bound.
{
  CGAL_qpe_assertion(i < qp_n); // Note: slack variables cannot have
				// finite upper bounds.
  CGAL_qpe_assertion(has_finite_upper_bound(i));
  return *(qp_u+i);
}

template < class Rep_ >
bool QP_solver<Rep_>::is_valid()
{
    switch(this->m_status) {
    case UPDATE:
      return false;
    case OPTIMAL:
      {
	const bool f = this->is_solution_feasible();
	const bool o = this->is_solution_optimal();
	CGAL_qpe_debug {
	  vout << "is in phase II: " << is_phaseII << std::endl;
	  vout << "feasible: " << f << std::endl;
	  vout << "optimal: " << o << std::endl;
	}
	return is_phaseII && f && o;
      }
    case INFEASIBLE:
      {
	const bool f = this->is_solution_feasible_for_auxiliary_problem();
	const bool o = this->is_solution_optimal_for_auxiliary_problem();
	const bool aux_positive = this->solution() > et0;
	CGAL_qpe_debug {
	  vout << "is in phase I: " << is_phaseI << std::endl;
	  vout << "feasible_aux: " << f << std::endl;
	  vout << "optimal_aux: " << o << std::endl;
	  vout << "objective value positive: " << aux_positive << std::endl;
	}
	return is_phaseI && f && o && aux_positive;
      }
    case UNBOUNDED:    
      {
	const bool f = this->is_solution_feasible();
	const bool u = this->is_solution_unbounded();
	CGAL_qpe_debug {
	  vout << "is in phase II: " << is_phaseII << std::endl;
	  vout << "feasible: " << f << std::endl;
	  vout << "unbounded: " << u << std::endl;
	}
	return is_phaseII && f && u;
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

  // check bounds on original- and artificial variables:
  Value_const_iterator v_it = x_B_O.begin();
  for (Index_const_iterator i_it = B_O.begin();
       i_it != B_O.end();
       ++i_it, ++v_it)
    if (*i_it < qp_n) {                 // original variable?
      if (has_finite_lower_bound(*i_it) && (*v_it < lower_bound(*i_it) * d) ||
	    has_finite_upper_bound(*i_it) && (*v_it > upper_bound(*i_it) * d))
        return false;
    } else                              // artificial variable?
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

  // compute left-hand side of (C1) (up to slackies and nonbasics):
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

  // compute left-hand side of (C1) (part for nonbasics):
  for (int i=0; i<qp_n; ++i)
    if (!is_basic(i)) {
      const ET var = nonbasic_original_variable_value(i) * d;
      for (int j=0; j<qp_m; ++j)
	lhs_col[j] += var * qp_A[i][j];
    }

  // compute left-hand side of (C1) (part for slackies):
  x_it = x_B_S.begin();
  for (Index_const_iterator i_it = B_S.begin();
       i_it != B_S.end();
       ++i_it, ++x_it) {
    const int k = *i_it - qp_n;
    lhs_col[slack_A[k].first] += (slack_A[ k].second ? -et1 : et1) * (*x_it);
  }

  // check equality (C1);
  for (int i=0; i<qp_m; ++i) {
    if ((qp_r[i] == Rep::EQUAL)         && (lhs_col[i] != ET(qp_b[i]) * d) ||
        (qp_r[i] == Rep::LESS_EQUAL)    && (lhs_col[i]  > ET(qp_b[i]) * d) ||
        (qp_r[i] == Rep::GREATER_EQUAL) && (lhs_col[i]  < ET(qp_b[i]) * d))
      return false;
  }
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
  const int no_of_wo_vars = this->number_of_working_variables() +
    (art_s_i == -2)? 1 : 0; // Note: if there ever was a special
                            // artifical variable, it has to be
                            // considered for this optimality test,
                            // even if it has already left the basis.
                            // (If there was initially a special
                            // artificial variable and it never left
                            // the basis, then number_of_-
                            // working_variables() counts it, so we
                            // add 0; it it left the basis (and thus
                            // art_s_i == -2) then number_of_-
                            // working_variables() does not count it,
                            // add we need to add 1.)
  
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
  if (!check_tag(Is_in_standard_form()))
    for (int i=0; i<qp_n; ++i)
      if (!is_basic(i))
	x_aux[i] = nonbasic_original_variable_value(i) * d;
  
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
      if (k < static_cast<int>(slack_A.size()))
                                     // slack variable
	tau_aux[col] = (slack_A[k].second? -et1 :  et1)
	  * lambda_aux[slack_A[k].first];
      else {                         // artificial variable
	k -= slack_A.size();
	if ((art_s_i == -1) || (col < no_of_wo_vars - 1))
	                             // normal artificial variable
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
  for (int col = 0; col < no_of_wo_vars; ++col)
    if (tau_aux[col] < et0 || 
	tau_aux[col] * x_aux[col] != et0)
	return false;

  return true;     
}

template < class Rep_ >
bool QP_solver<Rep_>::is_solution_feasible()
{
  Values lhs_col(qp_m, et0);
  for (int i=0; i<qp_n; ++i) {

    // check bounds on original variables:
    const ET var = is_basic(i)? x_B_O[in_B[i]] :  // should be ET &
      nonbasic_original_variable_value(i) * d;
    if (has_finite_lower_bound(i) && var < lower_bound(i) * d ||
	has_finite_upper_bound(i) && var > upper_bound(i) * d)
      return false;


    // compute A x times d:
    for (int j=0; j<qp_m; ++j)
      lhs_col[j] += var * qp_A[i][j];
  }
  
  // check A x = b (where in the code both sides are multiplied by d):
  for (int row = 0; row < qp_m; ++row) {
    const ET rhs = ET(qp_b[row])*d;
    if (qp_r[row] == Rep::EQUAL         && lhs_col[row] != rhs ||
	qp_r[row] == Rep::LESS_EQUAL    && lhs_col[row] >  rhs ||
	qp_r[row] == Rep::GREATER_EQUAL && lhs_col[row] <  rhs)
      return false;
  }
  
  return true;
}

template < class Rep_ >
bool QP_solver<Rep_>::
is_solution_optimal()
{
  CGAL_expensive_precondition(is_solution_feasible());

  // As described in documentation/UpperBounding.tex, the optimality
  // conditions for a QP that is not (necessarily) in standard form
  // read as follows. A feasible solution x is optimal if and only if
  // there exists \lambda and \tau such that
  //
  //    \tau^T = c^T + \lambda^T A + 2 D x                       (C3)
  //         0 = \lambda^T (A x - b)                             (C4)
  //
  // and
  //
  //      \lambda_i >= 0                for <=-constraints i,    (C5)
  //      \lambda_i <= 0                for >=-constraints i,    (C6)
  //
  // and
  //
  //             |  >= 0                if x_j = l_j < u_j,
  //     \tau_j <   =  0                if l_j < x_j < u_j,      (C7)
  //             |  <= 0                if x_j = u_j > l_j.
  //
  // We are going to check these conditions below.

  // get solution vector of original problem (multiplied by d):
  Values x(qp_n, et0);
  for (int i=0; i<qp_n; ++i)
    x[i] = is_basic(i)? x_B_O[in_B[i]] : nonbasic_original_variable_value(i)*d;

  // get \lambda:
  // Note: as the \lambda we have access to is actual d times the \lambda
  // from (C3), we will not compute \tau but \tau * d.
  Values lambda_prime(qp_m, et0);
  Value_const_iterator v_it = lambda.begin();
  for (Index_const_iterator i_it = C.begin();
       i_it != C.end(); ++i_it, ++v_it)
    lambda_prime[*i_it] = *v_it;
  
  // debug output:
  CGAL_qpe_debug {
    for (int col = 0; col < qp_m; ++col)
      vout5 << "lambda'[" << col << "]= " << lambda_prime[col] << std::endl;
  }
    
  // compute \tau^T = c^T + \lambda^T A + 2 D x (see conditions (C3) above):
  Values tau(qp_n,et0);
  for (int col = 0; col < qp_n; ++col) {
    tau[col] = ET(qp_c[col]) * d;
    for (int i=0; i<qp_m; ++i)
      tau[col] += lambda_prime[i] * qp_A[col][i];
  }
  if (!check_tag(Is_linear())) {
    for (int col = 0; col < qp_n; ++col) {
      D_pairwise_accessor twoD(qp_D,col);
      for (int i=0; i<qp_n; ++i)
	tau[col] += x[i] * twoD(i);
    }
  }

  for (int col = 0; col < qp_n; ++col)
    CGAL_qpe_debug {
      vout5 << "tau[" << col << "]= " << tau[col] << std::endl;
    }

  // compute A x:
  Values lhs_col(qp_m, et0);
  for (int i=0; i<qp_n; ++i)
    for (int j=0; j<qp_m; ++j)
      lhs_col[j] += x[i] * qp_A[i][j];

  // check (C5) and (C6), more precisely:
  // - check \lambda[i] >= 0 for i in LE:={j|qp_r[j]==LESS_EQUAL}
  // - check \lambda[i] <= 0 for i in GE:={j|qp_r[j]==GREATER_EQUAL}, and
  // - check \lambda'[i] * (Ax-b) == 0 for i in (LE union GE)
  for (int row = 0; row < qp_m; ++row)
    if (qp_r[row] != Rep::EQUAL) {
      const bool is_active = (lhs_col[row] == (ET(qp_b[row]) * d));
      if (qp_r[row] == Rep::LESS_EQUAL) {
	if (lambda_prime[row] < et0 ||
	    (lambda_prime[row] > et0) && !is_active)
	  return false;
      } else {           // GREATER_EQUAL?
	if (lambda_prime[row] > et0 ||
	    (lambda_prime[row] < et0) && !is_active)
	  return false;
      }
    }
  
  // check condition (C7) on \tau:
  for (int col = 0; col < qp_n; ++col) {
    const bool lower_tight = has_finite_lower_bound(col) &&
      ET(lower_bound(col) * d) == x[col];
    const bool upper_tight = has_finite_upper_bound(col) &&
      ET(upper_bound(col) * d) == x[col];
    if ( lower_tight && !upper_tight && tau[col] <et0 ||
	!lower_tight && !upper_tight && tau[col]!=et0 ||
        !lower_tight &&  upper_tight && tau[col] >et0) 
	return false;
  }

  return true;
}

template < class Rep_ >
bool QP_solver<Rep_>::is_solution_unbounded()
{
  // (This is documented in documentation/Test_suite.tex for the case
  // when the program is in standard form.)
  //
  // An "unbounded direction" w is defined to be a vector w such that
  // x_t(t):= x'-tw yields for any t>=0 feasible solutions with unbounded
  // objective value (assuming x' denotes the current solution).
  //
  // In order to check that the vector w returned by
  // unbounded_direction_begin() and unbounded_direction_end() is
  // indeed un anbounded direction, we need to check (i) feasibility
  // for all t>0 and (ii) that the objective value decreases while t
  // increases.  As to (i), we need
  //
  //     A x_t(t) <=> b   for all t
  //
  // which is equivalent to 
  //
  //     row j is GREATER_EQUAL then (Aw)_j <= 0,
  //     row j is         EQUAL then (Aw)_j  = 0,                (C8)
  //     row j is   LOWER_EQUAL then (Aw)_j >= 0.
  //
  // Feasibility of the constraints for the variables x is similar: we
  // need l <= x_t(t) <= u for all t, which is equivalent to
  //
  //     x_j has finite lower bound then w_j <= 0,               (C9)
  //     x_j has finite upper bound then w_j >= 0.               (C10)
  //
  // As to unboundedness (ii), we have (see equation (10) in
  // documentation/Test_suite.tex)
  //
  //     f(x_t(t)) = f(x') + t^2 w^TDw - t(c^T+2x^TD)w,
  //
  // that is, the objective function f behaves like a parabola.  So if
  // it should be unbounded then it has to be a lower parabola which
  // means w^TDw<=0, but D is positive semidefinite, so w^TDw = 0.
  // Thus, f must be linear in t, and unboundedness then implies
  // (c^T+2x^TD)w > 0.  Subsuming, (ii) requires
  //
  //      w^TDw = 0,                                             (C11)
  //      (c^T+2x^TD)w > 0.                                      (C12)

  CGAL_expensive_precondition(is_solution_feasible());
  
  // get solution vector of original problem (multiplied by d):
  Values x(qp_n, et0);
  for (int i=0; i<qp_n; ++i)
    x[i] = is_basic(i)? x_B_O[in_B[i]] : nonbasic_original_variable_value(i)*d;

  // check that the direction is not the zero vector:
  Unbounded_direction_iterator w = unbounded_direction_begin();
  bool all_zero = true;
  for (int i=0; i<qp_n; ++i) {
    if (w[i] != et0)
      all_zero = false;
    CGAL_qpe_debug {
      vout5 << "w[" << i << "]= " << w[i] << std::endl;
    }
  }
  if (all_zero)
    return false;

  // compute A w into aw:
  Values aw(qp_m, et0);
  for (int i=0; i<qp_n; ++i)
    for (int j=0; j<qp_m; ++j)
      aw[j] += w[i] * qp_A[i][j];

  // check feasibility (C8):
  for (int row=0; row<qp_m; ++row)
    if (qp_r[row] == Rep::GREATER_EQUAL && aw[row]  > et0 ||
	qp_r[row] == Rep::EQUAL         && aw[row] != et0 ||
	qp_r[row] == Rep::LESS_EQUAL    && aw[row]  < et0) 
      return false;

  // check feasibility (C9) and (C10):
  for (int i=0; i<qp_n; ++i)
    if (has_finite_lower_bound(i) && w[i] > et0 ||
	has_finite_upper_bound(i) && w[i] < et0)
      return false;

  // check unboundedness w^TDw=0 (C11):
  Values Dw(qp_n, et0);     // Note: will be reused for (C12) below.
  if (!check_tag(Is_linear())) {
    for (int i=0; i<qp_n; ++i)
      for (int j=0; j<qp_n; ++j)
	Dw[j] += w[i] * qp_D[j][i];
    ET sum = et0;           // will hold w^T * Dw...
    for (int i=0; i<qp_n; ++i)
      sum += w[i] * Dw[i];
    if (sum != et0)
      return false;
  }

  // check unboundedness (c^T+2x^TD)w > 0 (C12):
  ET m1 = et0, m2 = et0;
  Value_const_iterator x_it = x_B_O.begin();
  for (int i=0; i<qp_n; ++i)
    m1 += x[i] * Dw[i];
  m1 *= et2;                  // Note: m1 contains 2x^TDw*d (and not 2x^TDw).
  for (int i=0; i<qp_n; ++i)
    m2 += w[i] * qp_c[i];
  if (m2*d + m1 <= et0)
    return false;

  return true;
}

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
