// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp
//                 Kaspar Fischer

#include<CGAL/QP_functions.h>
#include<CGAL/NT_converter.h>

namespace CGAL {

// creation & initialization
// -------------------------
template < typename Q, typename ET, typename Tags >
QP_solver<Q, ET, Tags>::
QP_solver(const Q& qp, const Quadratic_program_options& options)
  : et0(0), et1(1), et2(2),
    strategyP(0),
    inv_M_B(vout4),
    d(inv_M_B.denominator()),
    m_phase(-1), is_phaseI(false), is_phaseII(false),
    is_RTS_transition(false),
    is_LP(check_tag(Is_linear())), is_QP(!is_LP),
    //no_ineq(check_tag(Has_equalities_only_and_full_rank())),
    no_ineq(QP_functions_detail::is_in_equational_form(qp)),
    // may change after phase I
    has_ineq(!no_ineq),
    is_nonnegative(check_tag(Is_nonnegative()))
{
  // init diagnostics
  diagnostics.redundant_equations = false;

  // initialization as in the standard-form case:
  set_verbosity(options.get_verbosity());
  // only if C_entry is double, we actually get filtered strategies,
  // otherwise we fall back to the respective non-filtered ones
  set_pricing_strategy(options.get_pricing_strategy());

  // Note: we first set the bounds and then call set() because set()
  // accesses qp_fl, qp_l, etc.
  set_explicit_bounds(qp);
  set(qp);

  // initialize and solve immediately:
  init();
  solve();
}

// set-up of QP

template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
set_D(const Q& /*qp*/, Tag_true /*is_linear*/)
{
  // dummy value, never used
  qp_D = 0;
}

template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
set_D(const Q& qp, Tag_false /*is_linear*/)
{
  qp_D = qp.get_d();
}

template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
set(const Q& qp)
{
  // assertions:
  CGAL_qpe_assertion(qp.get_n() >= 0);
  CGAL_qpe_assertion(qp.get_m() >= 0);

  // store QP
  qp_n = qp.get_n(); qp_m = qp.get_m();
  qp_A = qp.get_a(); qp_b = qp.get_b(); qp_c = qp.get_c(); qp_c0 = qp.get_c0();

  set_D(qp, Is_linear());
  qp_r = qp.get_r();

  // set up slack variables and auxiliary problem
  // --------------------------------------------

  // reserve memory for slack and artificial part of `A':
  if (has_ineq) {
    const unsigned int eq = static_cast<unsigned int>(std::count(qp_r, qp_r+qp_m, CGAL::EQUAL));
    slack_A.reserve(qp_m - eq);
    art_A.reserve  (       eq);
    art_s.insert(art_s.end(), qp_m, A_entry(0));
  } else
    art_A.reserve( qp_m);

  // decide on which bound the variables sit initially:
  if (!check_tag(Is_nonnegative()))
    init_x_O_v_i();

  set_up_auxiliary_problem();

  e = static_cast<int>(qp_m-slack_A.size()); // number of equalities
  l = (std::min)(qp_n+e+1, qp_m);  // maximal size of basis in phase I

  // diagnostic output:
  CGAL_qpe_debug {
    if (vout.verbose()) {
      if (vout2.verbose()) {
        vout2.out() << "======" << std::endl
                    << "Set-Up" << std::endl
                    << "======" << std::endl;
      }
    }
  }
  vout    << "[ " << (is_LP ? "LP" : "QP")
          << ", " << qp_n << " variables, " << qp_m << " constraints"
          << " ]" << std::endl;
  CGAL_qpe_debug {
      if (vout2.verbose() && (!slack_A.empty())) {
        vout2.out() << " (" << slack_A.size() << " inequalities)";
      }
      if (vout2.verbose()) {
        if (has_ineq)
          vout2.out() << "flag: has inequalities or rank not full"
                      << std::endl;
        if (vout4.verbose()) print_program();
      }
  }

  // set up pricing strategy:
  if (strategyP != static_cast< Pricing_strategy*>(0))
    strategyP->set(*this, vout2);

  // set up basis inverse:
  inv_M_B.set(qp_n, qp_m, e);

  // set phase:
  m_phase    = 0;
  is_phaseI  = false;
  is_phaseII = false;
}

template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
set_explicit_bounds(const Q& qp)
{
  set_explicit_bounds (qp, Is_nonnegative());
}

template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
set_explicit_bounds(const Q& /*qp*/, Tag_true) {
  // dummy values, never used
  qp_fl = 0;
  qp_l = 0;
  qp_fu = 0;
  qp_u = 0;
}

template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
set_explicit_bounds(const Q& qp, Tag_false) {
  qp_fl = qp.get_fl();
  qp_l = qp.get_l();
  qp_fu = qp.get_fu();
  qp_u = qp.get_u();
}



template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
init_x_O_v_i()
{
  // allocate storage:
  x_O_v_i.reserve(qp_n);
  x_O_v_i.resize (qp_n);

  // constants for comparisions:
  const L_entry l0(0);
  const U_entry u0(0);

  // our initial solution will have all original variables nonbasic,
  // and so we initialize them to zero (if the bound on the variable
  // allows it), or to the variable's lower or upper bound:
  for (int i = 0; i < qp_n; ++i) {
    CGAL_qpe_assertion( !*(qp_fl+i) || !*(qp_fu+i) || *(qp_l+i)<=*(qp_u+i));

    if (*(qp_fl+i))                    // finite lower bound?
      if (*(qp_fu+i))                  // finite lower and finite upper bound?
        if (*(qp_l+i) == *(qp_u+i))    // fixed variable?
          x_O_v_i[i] = FIXED;
        else                           // finite lower and finite upper?
          if (*(qp_l+i) <= l0 && u0 <= *(qp_u+i))
            x_O_v_i[i] = ZERO;
          else
            x_O_v_i[i] = LOWER;
      else                             // finite lower and infinite upper?
        if (*(qp_l+i) <= l0)
          x_O_v_i[i] = ZERO;
        else
          x_O_v_i[i] = LOWER;
    else                               // infinite lower bound?
      if (*(qp_fu+i))                  // infinite lower and finite upper?
        if (u0 <= *(qp_u+i))
          x_O_v_i[i] = ZERO;
        else
          x_O_v_i[i] = UPPER;
      else                             // infinite lower and infinite upper?
        x_O_v_i[i] = ZERO;
  }
}

template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
set_up_auxiliary_problem()
{
  ET            b_max(et0);
  const C_entry c1(1);
  int           i_max = -1; // i_max-th inequality is the most infeasible one
  int           i_max_absolute = -1; // absolute index of most infeasible ineq

  // TAG: TODO using variable i here, which is also the index of the entering
  // variable.
  for (int i = 0; i < qp_m; ++i) {
    // Note: For nonstandard form problems, our initial solution is not the
    // zero vector (but the vector with values original_variable_value(i),
    // 0<=i<qp_n), and therefore, rhs=b-Ax is not simply b as in the standard
    // form case, but Ax_init-b:
    const ET rhs = check_tag(Is_nonnegative())?
    ET(*(qp_b+i)) : ET(*(qp_b+i)) - multiply__A_ixO(i);

    if (has_ineq && (*(qp_r+i) != CGAL::EQUAL)) { // inequality constraint, so we
      // add a slack variable, and (if
      // needed) a special artificial
      if (*(qp_r+i) == CGAL::SMALLER) {        // '<='

        // add special artificial ('< -0') in case the inequality is
        // infeasible for our starting point (which is the origin):
        if (rhs < et0) {
          art_s[i] = -c1;
          if (-rhs > b_max) {
            i_max = static_cast<int>(slack_A.size());
            i_max_absolute = i;
            b_max = -rhs;
          }
        }


        // slack variable:
        slack_A.push_back(std::make_pair(i, false));
      } else {                                 // '>='

        // add special artificial ('> +0') in case the inequality is
        // infeasible for our starting point (which is the origin):
        if (rhs > et0) {
          art_s[i] = c1;
          if (rhs > b_max) {
            i_max = static_cast<int>(slack_A.size());
            i_max_absolute = i;
            b_max = rhs;
          }
        }

        // store slack column
        slack_A.push_back(std::make_pair(i, true));
      }

    } else {                                     // equality constraint, so we
      // add an artificial variable
      // (Note: if rhs==et0 here then the artificial variable is (at the
      // moment!) not needed. However, we nonetheless add it, for the following
      // reason. If we did and were given an equality problem with the zero
      // vector as the right-hand side then NO artificials would be added at
      // all; so our initial basis would be empty, something we do not want.)
      art_A.push_back(std::make_pair(i, rhs < et0));
    }
  } // end for

  // Note: in order to make our initial starting point (which is the origin) a
  // feasible point of the auxiliary problem, we need to initialize the
  // special artificial value correctly, namely to
  //
  //   max { |b_i| | i is index of an infeasible inequality constraint }. (C1)
  //
  // The index of this "most infeasible" constraint is, at this point of the
  // code, i_max (or i_max is -1 in which case all inequality constraints are
  // feasible and hence no special artififial column is needed at all).

  // prepare initialization of special artificial column:
  // Note: the real work is done in init_basis() below.
  if (i_max >= 0) {
    art_s_i = i_max;                           // Note: the actual
    art_basic = i_max_absolute;                // initialization of art_s_i
                                               // will be done in init_basis()
                                               // below. We misuse art_s_i to
                                               // remember i_max and art_basic
                                               // to remember i_max_absolute
  } else {                                     // no special art col needed
    art_s_i = -1;
    art_s.clear();
  }
}

// initialization (phase I)
template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
init()
{
  CGAL_qpe_debug {
    vout2 << std::endl
          << "==============" << std::endl
          << "Initialization" << std::endl
          << "==============" << std::endl;

  }

  // set status:
  m_phase    = 1;
  m_status   = QP_UPDATE;
  m_pivots   = 0;
  is_phaseI  = true;
  is_phaseII = false;

  // initial basis and basis inverse
  init_basis();

  // initialize additional data members
  init_additional_data_members();

  // initial solution
  init_solution();

  // initialize pricing strategy
  CGAL_qpe_assertion(strategyP != static_cast< Pricing_strategy*>(0));
  strategyP->init(0);

  // basic feasible solution already available?
  if (art_basic == 0) {

    // transition to phase II
    CGAL_qpe_debug {
      if (vout2.verbose()) {
        vout2.out() << std::endl
                    << "no artificial variables at all "
                    << "--> skip phase I"
                    << std::endl;
      }
    }
    transition();
  }
}

// Set up the initial basis and basis inverse.
template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
init_basis()
{
  int s_i = -1;
  int s_i_absolute = -1;
  const int s = static_cast<int>(slack_A.size());

  // has special artificial column?
  if (!art_s.empty()) {

    // Note: we maintain the information about the special artificial column in
    // the variable art_s_i and the vector s_art; in addition, however, we also
    // add a special "fake" column to art_A. This "fake" column has (in
    // constrast to the special artificial column) only one nonzero entry,
    // namely a +-1 for the most infeasible row (see (C1) above).

    // add "fake" column to art_A:
    s_i = art_s_i;               // s_i-th ineq. is most infeasible, see (C1)
    s_i_absolute = art_basic;    // absolute index of most infeasible ineq
    art_s_i = static_cast<int>(qp_n+s+art_A.size());    // number of special artificial var
    // BG: By construction of art_s_i (= i_max) in set_up_auxiliary_problem(),
    // s_i conforms with the indexing of slack_A, and the sign of the +-1
    // entry is just the negative of the corresponding slackie; this explains
    // the second parameter of make_pair. But the index passed as the
    // first parameter must refer to the ABSOLUTE index of the most
    // infeasible row. Putting s_i here is therefore a mistake unless
    // we only have equality constraints

    // art_A.push_back(std::make_pair(s_i, !slack_A[s_i].second));
    CGAL_qpe_assertion(s_i_absolute >= 0);
    CGAL_qpe_assertion(s_i_absolute == slack_A[s_i].first);
    art_A.push_back(std::make_pair(s_i_absolute, !slack_A[s_i].second));
  }

  // initialize indices of basic variables:
  if (!in_B.empty()) in_B.clear();
  in_B.reserve(qp_n+s+art_A.size());
  in_B.insert(in_B.end(), qp_n, -1);  // no original variable is basic

  init_basis__slack_variables(s_i, no_ineq);

  if (!B_O.empty()) B_O.clear();
  B_O.reserve(qp_n);                  // all artificial variables are basic
  for (int i = 0; i < static_cast<int>(art_A.size()); ++i) {
    B_O .push_back(qp_n+s+i);
    in_B.push_back(i);
  }
  art_basic = static_cast<int>(art_A.size());

  // initialize indices of 'basic' and 'nonbasic' constraints:
  if (!C.empty()) C.clear();
  init_basis__constraints(s_i, no_ineq);

  // diagnostic output:
  CGAL_qpe_debug {
    if (vout.verbose()) print_basis();
  }

  // initialize basis inverse (explain: 'art_s' not needed here (todo kf: don't
  // understand this note)):
  // BG: as we only look at the basic constraints, the fake column in art_A
  // will do as nicely as the actual column arts_s
  inv_M_B.init(static_cast<unsigned int>(art_A.size()), art_A.begin());
}

template < typename Q, typename ET, typename Tags >  inline                                 // no ineq.
void QP_solver<Q, ET, Tags>::
init_basis__slack_variables( int, Tag_true)
{
    // nop
}

template < typename Q, typename ET, typename Tags >                                        // has ineq.
void QP_solver<Q, ET, Tags>::
init_basis__slack_variables(int s_i, Tag_false)  // Note: s_i-th inequality is
                                                 // the most infeasible one,
                                                 // see (C1).
{
  const int s = static_cast<int>(slack_A.size());

  // reserve memory:
  if (!B_S.empty()) B_S.clear();
  B_S.reserve(s);

  // all slack variables are basic, except the slack variable corresponding to
  // special artificial variable (which is nonbasic): (todo kf: I do not
  // understand this)
  // BG: the s_i-th inequality is the most infeasible one, and the i-th
  // inequality corresponds to the slackie of index qp_n + i
  for (int i = 0; i < s; ++i) // go through all inequalities
    if (i != s_i) {
      in_B.push_back(static_cast<typename Indices::value_type>(B_S.size()));
      B_S .push_back(i+qp_n);
    } else
      in_B.push_back(-1);
}

template < typename Q, typename ET, typename Tags >  inline                                 // no ineq.
void QP_solver<Q, ET, Tags>::
init_basis__constraints( int, Tag_true)
{
  // reserve memory:
  C.reserve(qp_m);
  in_C.reserve(qp_m);

  // As there are no inequalities, C consists of all inequality constraints
  // only, so we add them all:
  for (int i = 0; i < qp_m; ++i) {
    C.push_back(i);
  }
}

template < typename Q, typename ET, typename Tags >                                        // has ineq.
void QP_solver<Q, ET, Tags>::
init_basis__constraints(int s_i, Tag_false)  // Note: s_i-th inequality is the
                                             // most infeasible one, see (C1).
{
  int i, j;

  // reserve memory:
  if (!in_C.empty()) in_C.clear();
  if (! S_B.empty())  S_B.clear();
  C.reserve(l);
  S_B.reserve(slack_A.size());

  // store constraints' indices:
  in_C.insert(in_C.end(), qp_m, -1);
  if (s_i >= 0) s_i = slack_A[s_i].first;    // now s_i is absolute index
                                             // of most infeasible row
  for (i = 0, j = 0; i < qp_m; ++i)
    if (*(qp_r+i) == CGAL::EQUAL) {             // equal. constraint basic
      C.push_back(i);
      in_C[i] = j;
      ++j;
    } else {                                  // ineq. constraint nonbasic
      if (i != s_i)                           // unless it's most infeasible
        S_B.push_back(i);
    }
    // now handle most infeasible inequality if any
  if (s_i >= 0) {
      C.push_back(s_i);
      in_C[s_i] = j;
  }
}

// Initialize r_C.
template < typename Q, typename ET, typename Tags >                 // Standard form
void  QP_solver<Q, ET, Tags>::
init_r_C(Tag_true)
{
}

// Initialize r_C.
template < typename Q, typename ET, typename Tags >                 // Upper bounded
void  QP_solver<Q, ET, Tags>::
init_r_C(Tag_false)
{
  r_C.resize(C.size());
  multiply__A_CxN_O(r_C.begin());
}

// Initialize r_S_B.
template < typename Q, typename ET, typename Tags >                 // Standard form
void  QP_solver<Q, ET, Tags>::
init_r_S_B(Tag_true)
{
}

// Initialize r_S_B.
template < typename Q, typename ET, typename Tags >                 // Upper bounded
void  QP_solver<Q, ET, Tags>::
init_r_S_B(Tag_false)
{
  r_S_B.resize(S_B.size());
  multiply__A_S_BxN_O(r_S_B.begin());
}

template < typename Q, typename ET, typename Tags >  inline                                 // no ineq.
void  QP_solver<Q, ET, Tags>::
init_solution__b_C(Tag_true)
{
  b_C.reserve(qp_m);
  std::transform(qp_b, qp_b+qp_m, std::back_inserter(b_C),
      NT_converter<B_entry,ET>());
}

template < typename Q, typename ET, typename Tags >  inline                                 // has ineq.
void  QP_solver<Q, ET, Tags>::
init_solution__b_C(Tag_false)
{
  b_C.insert(b_C.end(), l, et0);
  B_by_index_accessor  b_accessor(qp_b); // todo kf: is there some boost
                                         // replacement for this accessor?
  typedef typename std::iterator_traits<B_by_index_iterator>::value_type RT;
  std::transform(B_by_index_iterator(C.begin(), b_accessor),
                 B_by_index_iterator(C.end  (), b_accessor),
                 b_C.begin(),
                 NT_converter<RT,ET>());
}

// initial solution
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
init_solution()
{
  // initialize exact version of `qp_b' restricted to basic constraints C
  // (implicit conversion to ET):
  if (!b_C.empty()) b_C.clear();
  init_solution__b_C(no_ineq);

  // initialize exact version of `aux_c' and 'minus_c_B', the
  // latter restricted to basic variables B_O:
  if (!minus_c_B.empty()) minus_c_B.clear();
  minus_c_B.insert(minus_c_B.end(), l, -et1);   // todo: what is minus_c_B?
  CGAL_qpe_assertion(l >= static_cast<int>(art_A.size()));
  if (art_s_i > 0)
    minus_c_B[art_A.size()-1] *= ET(qp_n+qp_m); // Note: the idea here is to
                                                // give more weight to the
                                                // special artifical variable
                                                // so that it gets removed very
                                                // early, - todo kf: why?

  // ...and now aux_c: as we want to make all artificial variables (including
  // the special one) zero, we weigh these variables with >= 1 in the objective
  // function (and leave the other entries in the objective function at zero):
  aux_c.reserve(art_A.size());
  aux_c.insert(aux_c.end(), art_A.size(), 0);
  for (int col=static_cast<int>(qp_n+slack_A.size()); col<number_of_working_variables(); ++col)
    if (col==art_s_i)                           // special artificial?
      aux_c[col-qp_n-slack_A.size()]=  qp_n+qp_m;
    else                                        // normal artificial
      aux_c[col-qp_n-slack_A.size()]= 1;

  // allocate memory for current solution:
  if (!lambda.empty()) lambda.clear();
  if (!x_B_O .empty()) x_B_O .clear();
  if (!x_B_S .empty()) x_B_S .clear();
  lambda.insert(lambda.end(), l, et0);
  x_B_O .insert(x_B_O .end(), l, et0);
  x_B_S .insert(x_B_S .end(), slack_A.size(), et0);

  #if 0 // todo kf: I guess the following can be removed...
  //TESTING the updates of r_C, r_S_B, r_B_O, w
  //    ratio_test_bound_index = LOWER;
  //direction = 1;
  #endif

  // The following sets the pricing direction to "up" (meaning that
  // the priced variable will be increased and not decreased); the
  // statement is completely useless except that it causes debugging
  // output to be consistent in case we are running in standard form.
  // (If we are in standard form, the variable 'direction' is never
  // touched; otherwise, it will be set to the correct value during
  // each pricing step.)
  direction = 1;

  // initialization of vectors r_C, r_S_B:
  init_r_C(Is_nonnegative());
  init_r_S_B(Is_nonnegative());

  // compute initial solution:
  compute_solution(Is_nonnegative());

  // diagnostic output:
  CGAL_qpe_debug {
    if (vout.verbose()) print_solution();
  }
}

// Initialize additional data members.
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
init_additional_data_members()
{
  // todo kf: do we really have to insert et0, or would it suffice to just
  // resize() in the following statements?
  // BG: no clue, but it's at least safe that way

  if (!A_Cj.empty()) A_Cj.clear();
  A_Cj.insert(A_Cj.end(), l, et0);
  if (!two_D_Bj.empty()) two_D_Bj.clear();
  two_D_Bj.insert(two_D_Bj.end(), l, et0);

  if (!q_lambda.empty()) q_lambda.clear();
  q_lambda.insert(q_lambda.end(), l, et0);
  if (!q_x_O.empty()) q_x_O.clear();
  q_x_O.insert(q_x_O.end(), l, et0);
  if (!q_x_S.empty()) q_x_S.clear();
  q_x_S.insert(q_x_S.end(), slack_A.size(), et0);

  if (!tmp_l.empty()) tmp_l.clear();
  tmp_l.insert(tmp_l.end(), l, et0);
  if (!tmp_l_2.empty()) tmp_l_2.clear();
  tmp_l_2.insert(tmp_l_2.end(), l, et0);
  if (!tmp_x.empty()) tmp_x.clear();
  tmp_x.insert(tmp_x.end(), l, et0);
  if (!tmp_x_2.empty()) tmp_x_2.clear();
  tmp_x_2.insert(tmp_x_2.end(), l, et0);
}

} //namespace CGAL

// ===== EOF ==================================================================
