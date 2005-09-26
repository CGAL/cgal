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
	  vout << "feasible: " << f << std::endl;
	  vout << "optimal: " << o << std::endl;
	}
	return f && o;
      }
    case INFEASIBLE:
      {
	const bool f = this->is_solution_feasible_for_auxiliary_problem();
	const bool o = this->is_solution_optimal_for_auxiliary_problem();
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
    // check bounds on original variables:
    Value_const_iterator v_it = x_B_O.begin();
    for (Index_const_iterator i_it = B_O.begin();
	 i_it != B_O.end();
	 ++i_it, ++v_it)
      if (has_finite_lower_bound(*i_it) && *v_it < lower_bound(*i_it) ||
	  has_finite_upper_bound(*i_it) && *v_it > upper_bound(*i_it))
	return false;

    // compute A x times d:
    Values lhs_col(qp_m, et0);
    v_it = x_B_O.begin();
    for (Index_const_iterator i_it = B_O.begin();
	 i_it != B_O.end();
	 ++i_it, ++v_it)
      for (int i=0; i<qp_m; ++i)
	lhs_col[i] += *v_it * qp_A[*i_it][i];
    
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
void QP_solver<Rep_>::
is_solution_optimal__add_2_D_x(Values& tau,const Values& x,Tag_true) // LP
{
}

template < class Rep_ >
void QP_solver<Rep_>::
is_solution_optimal__add_2_D_x(Values& tau,const Values& x,Tag_false) // QP
{
  for (int col = 0; col < qp_n; ++col) {
    D_pairwise_accessor twoD(qp_D,col);
    for (int i=0; i<qp_n; ++i)
      tau[col] += x[i] * twoD(i);
  }
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
  Value_const_iterator v_it = x_B_O.begin();
  for (Index_const_iterator i_it = B_O.begin();
       i_it != B_O.end();
       ++i_it, ++v_it)
    x[*i_it] = *v_it;

  // get \lambda:
  // Note: as the \lambda we have access to is actual d times the \lambda
  // from (C3), we will not compute \tau but \tau * d.
  Values lambda_prime(qp_m, et0);
  v_it = lambda.begin();
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
  is_solution_optimal__add_2_D_x(tau,x,Is_linear()); // Note: we only
						     // need to do
						     // this in the
						     // QP-case, so we
						     // need a
						     // helper...
  for (int col = 0; col < qp_n; ++col)
    CGAL_qpe_debug {
      vout5 << "tau[" << col << "]= " << tau[col] << std::endl;
    }

  // compute A x:
  Values lhs_col(qp_m, et0);
  v_it = x_B_O.begin();
  for (Index_const_iterator i_it = B_O.begin();
       i_it != B_O.end();
       ++i_it, ++v_it) 
    for (int i=0; i<qp_m; ++i)
      lhs_col[i] += *v_it * qp_A[*i_it][i];

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
      ET(lower_bound(col)) == x[col];
    const bool upper_tight = has_finite_upper_bound(col) &&
      ET(upper_bound(col)) == x[col];
    if ( lower_tight && !upper_tight && tau[col] <et0 ||
	!lower_tight && !upper_tight && tau[col]!=et0 ||
        !lower_tight &&  upper_tight && tau[col] >et0)
	return false;
  }

  return true;
}

template < class Rep_ >
bool QP_solver<Rep_>::
is_solution_unbounded(Tag_false)		//QP case
{

    // doc:Test_suite.tex
    // Note: the basic feasible direction w is defined such that 
    // x'-tw yields feasible solutions where x' denotes the current
    // solution.
    // The following conditions for the basic feasible direction
    // w defined as w^T=[q_B_O^T|q_B_S^T|-e_j^T], j in N, must be met
    // 1.  w <= 0
    // 2. A * w = 0
    // 3. w^T * D * w > 0
    // 4. (c^T + 2 * x'^T * D) * w > 0
    //  
    Index_iterator i_it, M_i_it;
    Value_iterator v_it;
    Indices        M;
    bool           unbounded, basic_vars_nonpos;
    int sl_ind;
    
    CGAL_expensive_precondition(is_solution_feasible());
    
    //check nonpositivity of w_B_O
    basic_vars_nonpos = true;
    v_it = q_x_O.begin();
    for (i_it = B_O.begin(); i_it != B_O.end(); ++v_it, ++i_it) {
        basic_vars_nonpos = basic_vars_nonpos && ((*v_it) <= et0);
    }
    
    //check nonpositivity of w_B_S
    v_it = q_x_S.begin();
    for (i_it = B_S.begin(); i_it != B_S.end(); ++v_it, ++i_it) {
        basic_vars_nonpos = basic_vars_nonpos && ((*v_it) <= et0);
    }

    //check basic feasible direction w for A * w=0
    Values lhs_col(qp_m, et0);
    //initialize M
    for (int i = 0; i < qp_m; ++i) {
        M.push_back(i);
    }
    
    // A_CuS_B,B_O * w_B_O
    v_it = q_x_O.begin();
    for (i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) {
        A_by_index_accessor  a_accessor( qp_A[*i_it]);
	for (M_i_it = M.begin(); M_i_it != M.end(); ++M_i_it) {
	    lhs_col[*M_i_it] += (*v_it)
	        * (*A_by_index_iterator( M_i_it, a_accessor));
	}
    }
    
    // A_CuS_B,B_S * w_B_S
    v_it = q_x_S.begin();
    for (i_it = B_S.begin(); i_it != B_S.end(); ++i_it, ++v_it) {
        sl_ind = *i_it - qp_n;
	if (slack_A[sl_ind].second) {
	    lhs_col[slack_A[sl_ind].first] -= *v_it;
	} else {
	    lhs_col[slack_A[sl_ind].first] += *v_it;
	}
    }
    
    //A_CuS_B,N * w_N
    //
    // nonbasic nonzero component of w_N corresponds to original variable            
    if (j < qp_n) {
        A_by_index_accessor  a_accessor( qp_A[j]);
        for (M_i_it = M.begin(); M_i_it != M.end(); ++M_i_it) {
	    lhs_col[*M_i_it] -= d
	        * (*A_by_index_iterator( M_i_it, a_accessor));
        }
    // nonbasic nonzero component of w_N corresponds to slack variable	
    } else {
        sl_ind = j - qp_n;
	if (slack_A[sl_ind].second) {
	    lhs_col[slack_A[sl_ind].first] += d;
	} else {
	    lhs_col[slack_A[sl_ind].first] -= d;
	}
    }  
    unbounded = basic_vars_nonpos;    
    for (int row = 0; row < qp_m; ++row) {
        unbounded = unbounded && (lhs_col[row] == et0); 
    }
    
    // check w^T * D * w = 0
    // D_w contains D * w, will be reused later for the computation
    // of (c + 2 * x'^T * D) * w
    // Note: only original variables contribute
    Values D_w(qp_n, et0);
    for (int row = 0; row < qp_n; ++row) {
       v_it = q_x_O.begin();
       // nonbasic nonzero component of w_N corresponds to original variable
       if (j < qp_n) D_w[row] = -qp_D[row][j];
       for (i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) {
           D_w[row] += (*v_it) * qp_D[row][*i_it];
       }    
    }
    // w^T * D_w
    ET sum = et0;
    v_it =q_x_O.begin();
    for (i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) {
	sum += (*v_it) * D_w[*i_it];
    }
    // nonbasic nonzero component of w_N corresponds to original variable 
    if (j < qp_n) sum -= D_w[j];
    unbounded = unbounded && (sum == et0);
    
    // check c * w + 2 * x'^T * D * w > 0
    // reuse of D_w
    // Note: only original variables contribute
    sum = et0;
    v_it = x_B_O.begin();
    for (i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) { 
        sum += (*v_it) * D_w[*i_it];
    }
    sum *= et2;
    // c * w
    v_it = q_x_O.begin();
    for (i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) {
        sum += (*v_it) * qp_c[*i_it];
    }
    // nonbasic nonzero component of w_N corresponds to original variable
    if (j < qp_n) sum -= d * qp_c[j];
    unbounded = unbounded && (sum > et0);
    return unbounded;
}


template < class Rep_ >
bool QP_solver<Rep_>::
is_solution_unbounded(Tag_true)		//LP case
{

    // doc:Test_suite.tex
    // Note: the basic feasible direction w is defined such that 
    // x'-tw yields feasible solutions where x' denotes the current
    // solution.
    // The following conditions for the basic feasible direction
    // w defined as w^T=[q_B_O^T|q_B_S^T|-e_j^T], j in N, must be met
    // 1.  w<=0
    // 2. A * w=0
    // 3. c^T * w>0
    //
    Index_iterator i_it, M_i_it;
    Value_iterator v_it;
    Indices        M;
    bool           unbounded, basic_vars_nonpos;
    int sl_ind;
    
    CGAL_expensive_precondition(is_solution_feasible());
    
    //check nonpositivity of w_B_O
    basic_vars_nonpos = true;
    v_it = q_x_O.begin();
    for (i_it = B_O.begin(); i_it != B_O.end(); ++v_it, ++i_it) {
        basic_vars_nonpos = basic_vars_nonpos && ((*v_it) <= et0);
    }
    
    //check nonpositivity of w_B_S
    v_it = q_x_S.begin();
    for (i_it = B_S.begin(); i_it != B_S.end(); ++v_it, ++i_it) {
        basic_vars_nonpos = basic_vars_nonpos && ((*v_it) <= et0);
    }

    //check basic feasible direction w for A * w=0
    Values lhs_col(qp_m, et0);
    //initialize M
    for (int i = 0; i < qp_m; ++i) {
        M.push_back(i);
    }
    
    // A_CuS_B,B_O * w_B_O
    v_it = q_x_O.begin();
    for (i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) {
        A_by_index_accessor  a_accessor( qp_A[*i_it]);
	for (M_i_it = M.begin(); M_i_it != M.end(); ++M_i_it) {
	    lhs_col[*M_i_it] += (*v_it)
	        * (*A_by_index_iterator( M_i_it, a_accessor));
	}
    }
    
    // A_CuS_B,B_S * w_B_S
    v_it = q_x_S.begin();
    for (i_it = B_S.begin(); i_it != B_S.end(); ++i_it, ++v_it) {
        sl_ind = *i_it - qp_n;
	if (slack_A[sl_ind].second) {
	    lhs_col[slack_A[sl_ind].first] -= *v_it;
	} else {
	    lhs_col[slack_A[sl_ind].first] += *v_it;
	}
    }
    
    //A_CuS_B,N * w_N
    //
    // nonbasic nonzero component of w_N corresponds to original variable    
    if (j < qp_n) {
        A_by_index_accessor  a_accessor( qp_A[j]);
        for (M_i_it = M.begin(); M_i_it != M.end(); ++M_i_it) {
	    lhs_col[*M_i_it] -= d
	        * (*A_by_index_iterator( M_i_it, a_accessor));
        }
    // nonbasic nonzero component of w_N corresponds to slack variable
    } else {
        sl_ind = j - qp_n;
	if (slack_A[sl_ind].second) {
	    lhs_col[slack_A[sl_ind].first] += d;
	} else {
	    lhs_col[slack_A[sl_ind].first] -= d;
	}
    }  
    unbounded = basic_vars_nonpos;    
    for (int row = 0; row < qp_m; ++row) {
        unbounded = unbounded && (lhs_col[row] == et0); 
    }
        
    //check c^T*w > 0
    // Note: only original variables contribute 
    ET sum = et0;
    v_it = q_x_O.begin();
    for (i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) {
        sum += (*v_it) * qp_c[*i_it];
    }
    // nonbasic nonzero component of w_N corresponds to original variable
    if (j < qp_n) sum -= d * qp_c[j];
    unbounded = unbounded && (sum > et0);
    return unbounded;
}

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
