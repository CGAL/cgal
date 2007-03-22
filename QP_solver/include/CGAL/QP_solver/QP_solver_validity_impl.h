// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp <fransw@inf.ethz.ch>
//                 Kaspar Fischer <fischerk@inf.ethz.ch>

CGAL_BEGIN_NAMESPACE
template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::optimality_certificate_numerator(int i) const
{
  // we use the vector lambda which conforms to C (basic constraints)
  CGAL_qpe_precondition (i >= 0);
  CGAL_qpe_precondition (i <= qp_m);
  if (no_ineq)
    return lambda[i];
  else {
    int k = in_C[i];     // position of i in C
    if (k != -1) 
      return lambda[k];
    else 
      return et0;
  }   
}


template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::is_valid() const
{
  CGAL_qpe_debug {
    vout << std::endl;
    vout << "========" << std::endl;
    vout << "Validity" << std::endl;
    vout << "========" << std::endl;
  }
    switch(this->m_status) {
    case QP_UPDATE:
	CGAL_qpe_debug {
	  vout << " still in Update state!" << std::endl;
	}
      return false;
    case QP_OPTIMAL:
      {
	const bool f = this->is_solution_feasible();
	const bool o = this->is_solution_optimal();
	const bool v = this->is_value_correct();
	CGAL_qpe_debug {
	  vout << std::endl
	       << "----------" << std::endl
	       << "Validation" << std::endl
	       << "----------" << std::endl;
	  vout << " is in phase II: " << is_phaseII << std::endl;
	  vout << "       feasible: " << f << std::endl;
	  vout << "        optimal: " << o << std::endl;
	  vout << "  correct value: " << v << std::endl;
	}
	return is_phaseII && f && o && v;
      }
    case QP_INFEASIBLE:
      {
	const bool f = this->is_solution_feasible_for_auxiliary_problem();
	const bool o = this->is_solution_optimal_for_auxiliary_problem();
	const bool aux_positive = 
	  ( (this->solution_numerator() > et0) && (d > et0) );
	CGAL_qpe_debug {
	  vout << std::endl
	       << "----------" << std::endl
	       << "Validation" << std::endl
	       << "----------" << std::endl;
	  vout << "      is in phase I: " << is_phaseI << std::endl;
	  vout << "       feasible_aux: " << f << std::endl;
	  vout << "        optimal_aux: " << o << std::endl;
	  vout << " obj. val. positive: " << aux_positive << std::endl;
	}
	return is_phaseI && f && o && aux_positive;
      }
    case QP_UNBOUNDED:    
      {
	const bool f = this->is_solution_feasible();
	const bool u = this->is_solution_unbounded();
	CGAL_qpe_debug {
	  vout << std::endl
	       << "----------" << std::endl
	       << "Validation" << std::endl
	       << "----------" << std::endl;
	  vout << " is in phase II: " << is_phaseII << std::endl;
	  vout << "       feasible: " << f << std::endl;
	  vout << "      unbounded: " << u << std::endl;
	}
	return is_phaseII && f && u;
      }
    default: 	      
	CGAL_qpe_debug {
	  vout << " unknown state!" << std::endl;
	}
      return false;
    }
}

template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::is_solution_feasible_for_auxiliary_problem() const
{
  // some simple consistency checks:
  CGAL_qpe_assertion(is_phaseI);

  // check bounds on original- and artificial variables:
  Value_const_iterator v_it = x_B_O.begin();
  for (Index_const_iterator i_it = B_O.begin();
       i_it != B_O.end();
       ++i_it, ++v_it) {
    CGAL_qpe_assertion(B_O[in_B[*i_it]] == *i_it);
    if (*i_it < qp_n) {                 // original variable?
      if ((has_finite_lower_bound(*i_it) && (*v_it < lower_bound(*i_it) * d))||
	  (has_finite_upper_bound(*i_it) && (*v_it > upper_bound(*i_it) * d)))
        return false;
    } else                              // artificial variable?
      if (*v_it < et0)
        return false;
  }
  
  // check nonegativity of slack variables (the basic ones suffice):
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
  // below that the right-hand side is multiplied by d because we
  // maintain the denominator of the solution (which is d) separately.

  // compute left-hand side of (C1) (up to slackies and nonbasics):
  Values lhs_col(qp_m, et0);
  Value_const_iterator x_it = x_B_O.begin();
  for (Index_const_iterator i_it = B_O.begin();
       i_it != B_O.end();
       ++i_it, ++x_it)                    // iterate over all basic vars
    if (*i_it < qp_n)                     // ordinary original variable?
      for (int i=0; i<qp_m; ++i)
	lhs_col[i] += (*x_it) * ET(A_column(qp_A[*i_it])[i]);
    else                                  // artificial variable?
      if (*i_it != art_s_i)   {           // normal artificial variable?
	const int k = *i_it - qp_n - slack_A.size();
	lhs_col[art_A[k].first] += ET(art_A[ k].second ? -et1 : et1) * (*x_it);
      } else                              // special artificial variable
	for (int i=0; i<qp_m; ++i)
	  lhs_col[i] += (*x_it) * ET(art_s[i]);

  // compute left-hand side of (C1) (part for nonbasics):
  for (int j=0; j<qp_n; ++j)
    if (!is_basic(j)) {
      const ET var = nonbasic_original_variable_value(j) * d;
      for (int i=0; i<qp_m; ++i)
	lhs_col[i] += var * ET(A_column(qp_A[j])[i]);
    }

  // compute left-hand side of (C1) (part for slackies):
  x_it = x_B_S.begin();
  for (Index_const_iterator i_it = B_S.begin();
       i_it != B_S.end();
       ++i_it, ++x_it) {
    CGAL_qpe_assertion(B_S[in_B[*i_it]] == *i_it);
    const int k = *i_it - qp_n;
    lhs_col[slack_A[k].first] += ET(slack_A[ k].second ? -et1 : et1) * (*x_it);
  }

  // check equality (C1);
  for (int i=0; i<qp_m; ++i)
    if (lhs_col[i] != ET(qp_b[i]) * d)
      return false;
  return true;
}

template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::is_value_correct() const
{
  // checks whether solution_numerator() returns the right value
  // by computing 2 f(x) = x^T 2D x + 2x^T c + 2c0 from scratch; 
  // we have a common numerator d for the variables, so in our computations, 
  // - x^T 2D x carries a factor of d^2, 
  // - x^T 2c carries a factor of d,
  // - 2c0 carries no factor 
  CGAL_qpe_assertion(is_phaseII);
  int i = 0;
  ET z = et0; // represents x_i (2D_i x + 2c_i)
  for (Variable_numerator_iterator 
	 i_it = original_variables_numerator_begin(); 
       i_it < original_variables_numerator_end(); ++i_it, ++i) {
    if (*i_it == et0) continue; // no contribution from this variable
    ET s = et0; // represents 2D_i x + 2c_i
    int j = 0;
    if (is_QP) {
      Variable_numerator_iterator j_it = 
	original_variables_numerator_begin();
      // half the offdiagonal contribution
      for (;j<i; ++j_it, ++j)
         s += ET((*(qp_D+i))[j]) * *j_it; // 2D_ij x_j * d
      // the other half
      s *= et2;
      // the diagonal
      s += ET((*(qp_D+i))[j]) * *j_it;
    } 
    s += et2 * d * ET(qp_c[i]);       // add 2c_i * d
    z += s * *i_it;                   // add d*(2D_i x + 2 c_i) * d*x_i
  }
  // add constant term
  z += d * d * et2 * ET(qp_c0);       // now z = 2*d^2 * f(x) 
  return (z * solution_denominator() == d * d * et2 * solution_numerator());
}

template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::is_solution_optimal_for_auxiliary_problem() const
{
  // First, a note about artificials and how they need to be handled in this
  // optimality check. Observe that the (normal) artificial variables are
  // merely introduced to have an initial feasible solution: for each
  // infeasible equality, we introduce an artificial so that the equality gets
  // feasible. (Also, we need at least one artificial to have an initial
  // basis, see comment in set_up_auxiliary_problem().)  So once an artificial
  // drops to zero during phase I, we can argue as follows to prove that we
  // can FIX THIS ARTIFICIAL TO ZERO FOREVER (and, as a consequence, do not
  // need to price it in the sequel): since the artificial has dropped to
  // zero, we now have a point that fulfills the artificial's equality with
  // equality in the original problem, and so we can conceptually set up a NEW
  // auxiliary problem, in which no artificial is needed for the equality in
  // question.  This shows that there is no need to ever change an artificial
  // again once it has come down to zero (and that is why we do not price
  // artificials, see the pricing strategies, and why we do not need to
  // consider it for the optimality check).
  //
  // Note that the same observation also applies to the special artificial
  // variable.  Once the latter drops to zero, we have found a point that
  // fulfills all inequality constraints of the original problem.  So we can
  // again (conceptually) set up a new auxiliary problem --- this time without
  // a special artificial.
  //
  // With these considerations in mind, the check we perform below is the
  // following.  We test whether the auxiliary problem
  //
  //                   minimize     aux_c^T x                    (C13)
  //                   subject to   aux_A x =  b,
  //                                l <= x <= u,
  //
  // where aux_A = [ A A_art a_spec A_s], is optimal.  Here, the column a_spec
  // is only present in aux_A if the special artificial is existent and
  // nonzero; similarly, A_art only contains the columns of the (ordinary)
  // artificials that are nonzero.  Three cases may occur:
  //
  // - If the current solution does not optimally solve (C13), the test fails
  //   (and we have found a bug).
  //
  // - If the current solution optimally solves (C13) and has an objective
  //   value that is nonzero (i.e., larger than zero), we can conclude that
  //   the original problem input by the user is infeasible: if it were
  //   feasible there would be a solution x' to the initial auxiliary problem
  //   where all artificials (including the special one) are zero. But
  //   (C13) IS the initial auxiliary problem with the additional constraints
  //   that some auxiliaries (possibly including the special one) are fixed to
  //   zero; so x' would be a feasible point of (C13) with objective value 0,
  //   contradiction.
  //
  // - If the current solution optimally solves (C13) and has objective value
  //   zero, we have found a feasible point of the initial problem input by
  //   user. So we can continue with phase II.
  //
  // Note: if we checked optimality of (C13) with ALL artificial included (not
  // only the nonzero ones), this does not always work. The example
  // Ub_LP_Chvatal_8_1_d_QPE_solver.mps in the testsuite will then cause
  // problems: during some pricing step in phase I, the solver sees that no
  // non-artificial variable can be taken into the basis, and so it concludes
  // that the current solution optimally solves auxiliary problem. As the
  // objective function is nonzero, the original problem is recognized as
  // infeasible (which is correct), but the validiaty check fails because in
  // this instance, the entry in \tau from (C2) below for some artificial is
  // negative.
  //
  // Implementation of this check: using the KKT conditions, we obtain that a
  // feasible point x is optimal iff there exists a m-vector \lambda and a
  // |x|-vector \tau such that
  //
  //              \tau^T = aux_c^T + \lambda^T aux_A,            (C2)
  //
  //                      / >= 0   if x_j = l_j and l_j < u_j,
  //              \tau_j  |  = 0   if l_j < x_j < u_j,           (C2')
  //                      \ <= 0   if x_j = u_j and l_j < u_j,
  //
  // for all j. (This is Lemma 2 from documentation/UpperBounding.tex for the
  // special case where D=0.)
  //
  // These are the conditions we are going to check below. (Notice here that x
  // is the vector containing the original variables, the nonzero artificial
  // variables, the special artificial variable if it exists and is nonzero,
  // and the slack variables. In the code below, x will contain all ordinary
  // artifical variables, whether zero or not, and we will check \tau >= 0
  // only for those components of \tau that correspond to nonzero
  // aritificials.)

  // get number of original and (zero or nonzero, ordinary or special)
  // artificial variables:
  //
  // Note: if there was initially a special artificial variable and it never
  // left the basis, then number_of_working_variables() counts it; if it left
  // the basis (and thus art_s_i == -2) then number_of_working_variables()
  // does not count it.
  const int no_of_wo_vars = this->number_of_working_variables();
  CGAL_qpe_debug {
    vout5 << "number_of_working_variables: " << no_of_wo_vars << std::endl
	  << "art_s_i: " << art_s_i << std::endl << std::endl;
  }
  
  // collect solution vector of auxiliary problem:
  // todo: this calls for a nicer method to query the solution ...
  Values x_aux(no_of_wo_vars, et0);
  Value_const_iterator v_it = x_B_O.begin();
  for (Index_const_iterator i_it = B_O.begin();
       i_it != B_O.end(); ++i_it, ++v_it)
    x_aux[*i_it] = *v_it;
  v_it = x_B_S.begin();
  for (Index_const_iterator i_it = B_S.begin();
       i_it != B_S.end(); ++i_it, ++v_it)
    x_aux[*i_it] = *v_it;
  if (!check_tag(Is_nonnegative()))
    for (int j=0; j<qp_n; ++j)
      if (!is_basic(j))
	x_aux[j] = nonbasic_original_variable_value(j) * d;
  
  // Note: lambda[i] <= 0 for qp_r[i] == "GREATER_EQUAL"
  // todo: (ask frans) what does the above note mean here?
  Values lambda_aux(qp_m, et0);
  v_it = lambda.begin();
  for (Index_const_iterator i_it = C.begin();
       i_it != C.end(); ++i_it, ++v_it)
    {
      lambda_aux[*i_it] = *v_it;
    }

  // output for debugging:
  CGAL_qpe_debug {
    for (int i = 0; i < no_of_wo_vars; ++i)
      vout5 << "x_aux[" << i << "]= " << x_aux[i] << std::endl;
    for (int col = 0; col < qp_m; ++col)
      vout5 << "lambda_aux[" << col << "]= " << lambda_aux[col] << std::endl;
  }
  
  // compute \tau^T = aux_c^T + \lambda^T * aux_A (see conditions (C2) above):
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
	tau_aux[col] += lambda_aux[i] * ET(A_column(qp_A[col])[i]);
    else {
      int k = col - qp_n;
      if (k < static_cast<int>(slack_A.size()))
                                     // slack variable
	tau_aux[col] += ET(slack_A[k].second? -et1 :  et1)
	  * lambda_aux[slack_A[k].first];
      else {                         // artificial variable
	k -= slack_A.size();
	if (art_s_i == -1 ||         // no spec. art. ever => all art. normal
	    art_s_i == -2 ||         // spec. art. out now => all art. normal
	    col < no_of_wo_vars-1)   // spec. art still here => check
	                             // case of normal artificial variable
	  tau_aux[col] += ET(art_A[ k].second? -et1 : et1)
	    * lambda_aux[art_A[k].first];
	else                         // case of special artificial variable
	  for (int i=0; i<qp_m; ++i)
	    tau_aux[col] += lambda_aux[i] * ET(art_s[i]);
      }
    }
    CGAL_qpe_debug {
      vout5 << "tau_aux[" << col << "]= " << tau_aux[col] << std::endl;
    }
  }

  // check (C2'):
  for (int col = 0; col < no_of_wo_vars; ++col) {
    // actually, basic variables should have tau == 0
    CGAL_qpe_assertion(!is_basic(col) ||  tau_aux[col] == 0);
    if (!is_artificial(col) || x_aux[col] != et0) { // is it a slack or
						    // original variable, or a
						    // nonzero aritificial?
      const Bnd l_bnd = lower_bnd(col) * d;
      const Bnd u_bnd = upper_bnd(col) * d;
      if ((l_bnd == x_aux[col] && l_bnd < u_bnd && tau_aux[col] < et0) ||
	  (l_bnd  < x_aux[col] && u_bnd > x_aux[col] && tau_aux[col] != et0) ||
	  (u_bnd == x_aux[col] && l_bnd < u_bnd && tau_aux[col] > et0))
	return false;
    }
  }

  return true;
}

template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::is_solution_feasible() const
{  
  // some simple consistency checks:
  CGAL_qpe_assertion(is_phaseII);

  Values lhs_col(qp_m, et0);
  Variable_numerator_iterator it = original_variables_numerator_begin();
  for (int i=0; i<qp_n; ++i, ++it) {

    // check bounds on original variables:
    const ET var = is_basic(i)? x_B_O[in_B[i]] :  // should be ET &
      nonbasic_original_variable_value(i) * d;
    if (var != *it)
      return false; // original_variables_numerator_begin() inconsistent
    if ((has_finite_lower_bound(i) && var < lower_bound(i) * d) ||
	(has_finite_upper_bound(i) && var > upper_bound(i) * d))
      return false;

    // compute A x times d:
    for (int j=0; j<qp_m; ++j)
      lhs_col[j] += var * ET(A_column(qp_A[i])[j]);
  }
  
  // check A x = b (where in the code both sides are multiplied by d):
  for (int row = 0; row < qp_m; ++row) {
    const ET rhs = ET(qp_b[row])*d;
    if ((qp_r[row] == CGAL::EQUAL         && lhs_col[row] != rhs) ||
	(qp_r[row] == CGAL::SMALLER    && lhs_col[row] >  rhs) ||
	(qp_r[row] == CGAL::LARGER && lhs_col[row] <  rhs))
      return false;
  }
  
  return true;
}

template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::
is_solution_optimal() const
{
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
  Optimality_certificate_numerator_iterator 
    itb = optimality_certificate_numerator_begin();
  Optimality_certificate_numerator_iterator 
    ite = optimality_certificate_numerator_end();
  Optimality_certificate_iterator 
    itq = optimality_certificate_begin();
  Values lambda_prime;
  for (Optimality_certificate_numerator_iterator 
	 it = itb; it != ite; ++it, ++itq) {
    CGAL_qpe_assertion (*it * (*itq).denominator() == (*itq).numerator() * d);
    lambda_prime.push_back(*it);
  }
  
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
      tau[col] += lambda_prime[i] * ET(A_column(qp_A[col])[i]);
  }
  if (!check_tag(Is_linear())) {
    for (int col = 0; col < qp_n; ++col) {
      D_pairwise_accessor twoD(qp_D,col);
      for (int i=0; i<qp_n; ++i)
	tau[col] += x[i] * ET(twoD(i));
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
      lhs_col[j] += x[i] * ET(A_column(qp_A[i])[j]);

  // check (C5) and (C6), more precisely:
  // - check \lambda[i] >= 0 for i in LE:={j|qp_r[j]==LESS_EQUAL}
  // - check \lambda[i] <= 0 for i in GE:={j|qp_r[j]==GREATER_EQUAL}, and
  // - check \lambda'[i] * (Ax-b) == 0 for i in (LE union GE)
  for (int row = 0; row < qp_m; ++row)
    if (qp_r[row] != CGAL::EQUAL) {
      const bool is_active = (lhs_col[row] == (ET(qp_b[row]) * d));
      if (qp_r[row] == CGAL::SMALLER) {
	if (lambda_prime[row] < et0 ||
	    ((lambda_prime[row] > et0) && !is_active))
	  return false;
      } else {           // GREATER_EQUAL?
	if (lambda_prime[row] > et0 ||
	    ((lambda_prime[row] < et0) && !is_active))
	  return false;
      }
    }
  
  // check condition (C7) on \tau:
  for (int col = 0; col < qp_n; ++col) {
    const bool lower_tight = has_finite_lower_bound(col) &&
      ET(lower_bound(col) * d) == x[col];
    const bool upper_tight = has_finite_upper_bound(col) &&
      ET(upper_bound(col) * d) == x[col];
    if ( (lower_tight && !upper_tight && tau[col] <et0) ||
	(!lower_tight && !upper_tight && tau[col]!=et0) ||
	 (!lower_tight &&  upper_tight && tau[col] >et0)) 
	return false;
  }

  return true;
}

template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::is_solution_unbounded() const
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
      aw[j] += w[i] * ET(A_column(qp_A[i])[j]);

  // check feasibility (C8):
  for (int row=0; row<qp_m; ++row)
    if ((qp_r[row] == CGAL::LARGER && aw[row]  > et0) ||
	(qp_r[row] == CGAL::EQUAL         && aw[row] != et0) ||
	(qp_r[row] == CGAL::SMALLER    && aw[row]  < et0)) 
      return false;

  // check feasibility (C9) and (C10):
  for (int i=0; i<qp_n; ++i)
    if ((has_finite_lower_bound(i) && w[i] > et0) ||
	(has_finite_upper_bound(i) && w[i] < et0))
      return false;

  // check unboundedness 2 Dw=0 (C11):
  if (!check_tag(Is_linear()))
    for (int i=0; i<qp_n; ++i) {
      ET sum = et0;
      D_pairwise_accessor twoD(qp_D,i);
      for (int j=0; j<qp_n; ++j)
	sum += w[j] * ET(twoD(j));
      if (sum != et0)
	return false;
    }
  
  // check unboundedness c^Tw > 0 (C12):
  ET m = et0;
  for (int i=0; i<qp_n; ++i)
    m += w[i] * ET(qp_c[i]);
  if (m <= et0)
    return false;

  return true;
}

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
