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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>
#ifndef CGAL_QP_SOLUTION_IMPL_H
#define CGAL_QP_SOLUTION_IMPL_H

#include <CGAL/license/QP_solver.h>


#include<iterator>

namespace CGAL {
// checks whether the solution actually solves program p  
// this performs exactly the checks described in the doc
// of the class Quadratic_program_solution
// -----------------------------------------------------
template <typename ET>
template <class Program, typename Is_linear, typename Is_nonnegative>
bool Quadratic_program_solution<ET>::solves_program 
(const Program& p, 
 Is_linear is_linear, Is_nonnegative is_nonnegative)
{    
  CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");

  // first check whether the dimensions agree
  int n = static_cast<int>(variable_numerators_end() - variable_numerators_begin());
  if (n != p.get_n()) 
    return error ("wrong number of variables");
  int m = static_cast<int>(solver()->lambda_end() - solver()->lambda_begin());
  if (m != p.get_m()) 
    return error("wrong number of constraints");

  // now distinguish between the three possible cases
  switch (status()) {
  case QP_OPTIMAL: {
    std::vector<ET> ax_minus_b (m, et0); // d(Ax-b)
    std::vector<ET> two_Dx (n, et0);     // d(2Dx)
    return 
      // the following fills ax_minus_b...
      is_feasible (p, ax_minus_b, is_nonnegative) &&          // feasible?
      // the following fills two_Dx...
      is_value_correct (p, two_Dx, is_linear) &&              // obj. value?
      is_optimal_1 (p) &&                                     // condition 1?
      // ...and the following uses ax_minus_b
      is_optimal_2 (p, ax_minus_b)         &&                 // condition 2?
      // ...and the following uses two_Dx
      is_optimal_3 (p, two_Dx, is_linear, is_nonnegative);    // condition 3?
  }
  case QP_INFEASIBLE: {
    std::vector<ET> lambda_a (n, et0); // lambda^TA
    return 
      is_infeasible_1 (p) &&                                  // condition 1?
      // the following fills lambda_a...
      is_infeasible_2 (p, lambda_a, is_nonnegative) &&        // condition 2?
      // ...and the following uses lambda_a
      is_infeasible_3 (p, lambda_a, is_nonnegative);          // condition 3?
  }
  case QP_UNBOUNDED: {
    std::vector<ET> ax_minus_b (m, et0);
    return 
      is_feasible (p, ax_minus_b, is_nonnegative) &&          // feasible?
      is_unbounded_1 (p) &&                                   // condition 1?
      is_unbounded_2 (p, is_nonnegative) &&                   // condition 2?
      is_unbounded_3 (p, is_linear);                          // condition 3?
  }
  default:
    return error ("solution in undefined state");
  }
}

// tests whether Ax ~ b is satisfied and computes d(Ax-b);
// precondition: ax_minus_b has length m and is zero
// ------------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::are_constraints_feasible 
(const Program& p, 
 typename std::vector<ET>& ax_minus_b)
{
  typedef typename Program::B_iterator B_column_iterator;
  typedef typename Program::R_iterator R_column_iterator;

  // fill ax_minus_b (length m)
  int m = p.get_m();
  B_column_iterator b = p.get_b();
  add_Az (p, variable_numerators_begin(), ax_minus_b);     // d(Ax)
  ET d = variables_common_denominator();
  for (int i=0; i<m; ++i, ++b)
    ax_minus_b[i] -= ET (*b) * d;                          // d(Ax-b)
  
  // now validate constraints 
  if (d <= et0) 
    return error ("common variable denominator is negative");
  R_column_iterator r = p.get_r();
  for (int i=0; i<m; ++i, ++r)
    switch (*r) {
    case SMALLER: // <=
      if (ax_minus_b[i] > et0) return error("inequality (<=) violated");
      break;
    case EQUAL:   // ==
      if (ax_minus_b[i] != et0) return error("inequality (==) violated");
      break;
    case LARGER:  // >=
      if (ax_minus_b[i] < et0) return error("inequality (>=) violated"); 
      break;
    }
  return true;
}

// tests whether Ax ~ b, x >= 0 is satisfied and computes d(Ax-b)
// precondition: ax_minus_b has length m and is zero
// --------------------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_feasible 
(const Program& p, 
 typename std::vector<ET>& ax_minus_b,
 Tag_true /*is_nonnegative*/)
{
  // test Ax ~ b
  if (!are_constraints_feasible (p, ax_minus_b)) return false;
  // test x >= 0  
  if (variables_common_denominator() <= et0) 
    return error ("common variable denominator is negative");
  for (Variable_numerator_iterator v = variable_numerators_begin();
       v < variable_numerators_end(); ++v)
    if (*v < et0) return error("bound (>= 0) violated");
  return true;
}

// tests whether Ax ~ b, l <= x <= u is satisfied and computes d(Ax-b)
// precondition: ax_minus_b has length m and is zero
// -------------------------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_feasible 
(const Program& p, 
 typename std::vector<ET>& ax_minus_b,
 Tag_false /*is_nonnegative*/)
{
  // test Ax ~ b
  if (!are_constraints_feasible (p, ax_minus_b)) return false; // (*)
  // test l <= x <= u   
  ET d = variables_common_denominator();
  if (d <= et0) 
    return error ("common variable denominator is negative");
  typedef typename Program::FL_iterator FL_column_iterator;
  typedef typename Program::L_iterator L_column_iterator;  
  typedef typename Program::FU_iterator FU_column_iterator;
  typedef typename Program::U_iterator U_column_iterator;
  FL_column_iterator fl = p.get_fl();
  FU_column_iterator fu = p.get_fu();
  L_column_iterator l = p.get_l();
  U_column_iterator u = p.get_u();  
  for (Variable_numerator_iterator v = variable_numerators_begin();
       v < variable_numerators_end(); ++v, ++fl, ++l, ++fu, ++u) {
    if (*fl && (*v < ET(*l) * d)) 
      return error("bound (>=l) violated");
    if (*fu && (*v > ET(*u) * d)) 
      return error("bound (<=u) violated");
  }
  return true;
}

// checks whether constraint type <= (>=) implies lambda >= (<=) 0
// --------------------------------------------------------------- 
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_optimal_1 
(const Program& p)
{
  if (variables_common_denominator() <= et0) 
    return error ("wrong sign of optimality certificate");
  typedef typename Program::R_iterator R_column_iterator;
  int m = p.get_m();
  R_column_iterator r = p.get_r();
  Optimality_certificate_numerator_iterator opt = 
    optimality_certificate_numerators_begin();
  for (int i=0; i<m; ++i, ++r, ++opt) {
    if (*r == SMALLER &&  *opt < et0)
      return error("constraint (<=) with lambda < 0");
    if (*r == LARGER &&  *opt > et0)
      return error("constraint (>=) with lambda > 0");
  }
  return true;
}

// checks whether constraint type <= (>=) implies lambda >= (<=) 0
// --------------------------------------------------------------- 
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_infeasible_1 
(const Program& p)
{ 
  int m = p.get_m();
  typedef typename Program::R_iterator R_column_iterator;
  R_column_iterator r = p.get_r();
  Infeasibility_certificate_iterator inf = infeasibility_certificate_begin();
  for (int i=0; i<m; ++i, ++r, ++inf) {
    if (*r == SMALLER &&  *inf < et0)
      return error("constraint (<=) with lambda < 0");
    if (*r == LARGER &&  *inf > et0)
      return error("constraint (>=) with lambda > 0");
  }
  return true;
}

// checks whether lambda and Ax-b are complementary
// ------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_optimal_2 
(const Program& p, 
 const typename std::vector<ET>& ax_minus_b)
{
  int m = p.get_m();  
  Optimality_certificate_numerator_iterator opt = 
    optimality_certificate_numerators_begin();
  for (int i=0; i<m; ++i, ++opt)
    if (*opt != et0 && ax_minus_b[i] != et0) 
      return error("lambda and Ax-b are not complementary");
  return true;					     
}

// checks whether (c^T + lambda^TA)_j >= (==) 0 if x_j = (>) 0
// -----------------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_optimal_3 
(const Program& p, std::vector<ET>& q,
 Tag_true /*is_linear*/, Tag_true /*is_nonnegative*/)
{
  int n = p.get_n();                                           // q = 0
  add_c (p, q);                                                // d c^T
  add_zA (p, optimality_certificate_numerators_begin(), q);    // d lambda^T A
  Variable_numerator_iterator v = variable_numerators_begin();
  for (int j=0; j<n; ++j, ++v) {
    if (q[j] < et0) 
      return error("some (c^T + lambda^TA)_j is negative");
    if (*v > et0 && q[j] > et0)
      return error("(c^T + lambda^TA) and x are not complementary");
  }
  return true;
}

// checks whether (c^T + lambda^TA + 2x^TD)_j >= (==) 0 if x_j = (>) 0
// -------------------------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_optimal_3 
(const Program& p, std::vector<ET>& q,
 Tag_false /*is_linear*/, Tag_true /*is_nonnegative*/)
{
  int n = p.get_n();                                          // q = 2Dx
  add_c (p, q);                                               // d * c^T
  add_zA (p, optimality_certificate_numerators_begin(), q);   // d * lambda^T A
  Variable_numerator_iterator v = variable_numerators_begin();
  for (int j=0; j<n; ++j, ++v) {
    if (q[j] < et0) 
      return error("some (c^T + lambda^TA + 2Dx)_j is negative");
    if (*v > et0 && q[j] > et0)
      return error("(c^T + lambda^TA + 2Dx) and x are not complementary");
  }
  return true;
}

// checks whether (c^T + lambda^TA )_j >= (==, <=) 0 if 
// x_j = l_j < u_j (l_j < x_j < u_j,  l_j < u_j = x_j)
// ---------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_optimal_3 
(const Program& p, std::vector<ET>& q, 
 Tag_true /*is_linear*/, Tag_false /*is_nonnegative*/)
{
  int n = p.get_n();                                          // q = 0
  add_c (p, q);                                               // d * c^T
  add_zA (p, optimality_certificate_numerators_begin(), q);   // d * lambda^T A
  typedef typename Program::FL_iterator FL_column_iterator;
  typedef typename Program::L_iterator L_column_iterator;  
  typedef typename Program::FU_iterator FU_column_iterator;
  typedef typename Program::U_iterator U_column_iterator;
  FL_column_iterator fl = p.get_fl();
  FU_column_iterator fu = p.get_fu();
  L_column_iterator l = p.get_l();
  U_column_iterator u = p.get_u();  
  Variable_numerator_iterator v = variable_numerators_begin();  
  ET d = variables_common_denominator();
  if (d <= et0)
    return error ("common variable denominator is negative");
  for (int j=0; j<n; ++j, ++v, ++fl, ++l, ++fu, ++u) {
    if (*fl && *v == ET(*l) * d && (!*fu || *l < *u) && q[j] < et0) 
      return error("x_j = l_j < u_j but (c^T + lambda^TA )_j < 0");
    if ( (!*fl || *v > ET(*l) * d) && (!*fu || *v < ET(*u) * d) && q[j] != et0)
      return error("l_j < x_j < u_j but (c^T + lambda^TA )_j != 0");
    if (*fu && *v == ET(*u) * d && (!*fl || *l < *u) && q[j] > et0) 
      return error("x_j = u_j > l_j but (c^T + lambda^TA )_j > 0");
  }
  return true;
}

// checks whether (c^T + lambda^TA +2x^*D)_j >= (==, <=) 0 if 
// x_j = l_j < u_j (l_j < x_j < u_j,  l_j < u_j = x_j)
// ----------------------------------------------------------
template <typename ET>
template <typename Program>
bool  Quadratic_program_solution<ET>::is_optimal_3 
(const Program& p, std::vector<ET>& q, 
 Tag_false /*is_linear*/, Tag_false /*is_nonnegative*/)
{
  int n = p.get_n();                                         // q = d * 2Dx
  add_c (p, q);                                              // d * c^T
  add_zA (p, optimality_certificate_numerators_begin(), q);  // d * lambda^T A
  typedef typename Program::FL_iterator FL_column_iterator;
  typedef typename Program::L_iterator L_column_iterator;  
  typedef typename Program::FU_iterator FU_column_iterator;
  typedef typename Program::U_iterator U_column_iterator;
  FL_column_iterator fl = p.get_fl();
  FU_column_iterator fu = p.get_fu();
  L_column_iterator l = p.get_l();
  U_column_iterator u = p.get_u();  
  Variable_numerator_iterator v = variable_numerators_begin();  
  ET d = variables_common_denominator();
  if (d <= et0)
    return error ("common variable denominator is negative");
  for (int j=0; j<n; ++j, ++v, ++fl, ++l, ++fu, ++u) {
    if (*fl && *v == ET(*l) * d && (!*fu || *l < *u) && q[j] < et0) 
      return error("x_j = l_j < u_j but (c^T + lambda^TA + 2Dx)_j < 0");
    if ( (!*fl || *v > ET(*l) * d) && (!*fu || *v < ET(*u) * d) && q[j] != et0)
      return error("l_j < x_j < u_j but (c^T + lambda^TA + 2Dx)_j != 0");
    if (*fu && *v == ET(*u) * d && (!*fl || *l < *u) && q[j] > et0) 
      return error("x_j = u_j > l_j but (c^T + lambda^TA + 2Dx)_j > 0");
  }
  return true;
}

// checks whether objective function value is c^T x + c_0, but
// leaves the vector q alone
// -------------------------------------------------------------- 
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_value_correct 
(const Program& p, std::vector<ET>& /*q*/, Tag_true /*is_linear*/)
{
  // check objective value c^T x + c_0
  ET d = variables_common_denominator();
  if (d <= et0)
    return error ("common variable denominator is negative");
  Variable_numerator_iterator v = variable_numerators_begin();
  ET obj = d  * ET(p.get_c0());                             // d * c_0
  typedef typename Program::C_iterator C_column_iterator;
  C_column_iterator c = p.get_c();
  int n = p.get_n();
  for (int j=0; j<n; ++j, ++v, ++c)
    obj += ET(*c) * *v;                                     // d * c^T x
  if (Quotient<ET>(obj, d) != objective_value())
    return error ("optimal objective value c^T x + c_0 incorrect");
  return true;
}

// checks whether objective function value is x^TDx + c^T x + c_0
// and fills the vector q with 2Dx
// -------------------------------------------------------------- 
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_value_correct 
(const Program& p, std::vector<ET>& q, Tag_false /*is_linear*/)
{
  ET d = variables_common_denominator();
  if (d <= et0)
    return error ("common variable denominator is negative");
  Variable_numerator_iterator v = variable_numerators_begin();  
  int n = p.get_n();
  add_two_Dz (p, v, q);                                     // 2d * Dx
  ET et2(2);
  ET obj = et2 * d * d * ET(p.get_c0());                    // 2d^2 * c_0
  typedef typename Program::C_iterator C_column_iterator;
  C_column_iterator c = p.get_c();  
  for (int j=0; j<n; ++j, ++v, ++c) {
    obj += et2 * d * ET(*c) * *v;                           // 2d^2 * c^T x
    obj += q[j] * *v;                                       // 2d^2 * x^TDx
  }  
  if (Quotient<ET>(obj, et2*d*d) != objective_value())
    return error ("optimal objective value x^TDx + c^T x + c_0 incorrect");
  return true;
}

// checks whether lambda^TA >= 0
// -----------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_infeasible_2 
(const Program& p, typename std::vector<ET>& lambda_a, 
 Tag_true /*is_nonnegative*/)
{
  // fill lambda_a
  add_zA (p, infeasibility_certificate_begin(), lambda_a); 
  int n = p.get_n();
  for (int j=0; j<n; ++j)
    if (lambda_a[j] < et0)
      return error ("(lambda^TA)_j < 0 for some j");
  return true;
}

// checks whether lambda^TA >= 0 (<= 0) if u_j = infty (l_j = -infty)
// ------------------------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_infeasible_2 
(const Program& p, typename std::vector<ET>& lambda_a, 
 Tag_false /*is_nonnegative*/)
{
  // fill lambda_a
  add_zA (p, infeasibility_certificate_begin(), lambda_a); 
  int n = p.get_n();
  typedef typename Program::FL_iterator FL_column_iterator;
  typedef typename Program::FU_iterator FU_column_iterator;
  FL_column_iterator fl = p.get_fl();
  FU_column_iterator fu = p.get_fu();
  for (int j=0; j<n; ++j, ++fl, ++fu) {
    if (!*fu && lambda_a[j] < et0)
      return error ("u_j = infty but (lambda^TA)_j < 0");
    if (!*fl && lambda_a[j] > et0)
      return error ("l_j = -infty but (lambda^TA)_j > 0");
  }
  return true;
}

// checks whether lambda^T b < 0
// ----------------------------- 
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_infeasible_3 
(const Program& p, const typename std::vector<ET>& /*lambda_a*/,
 Tag_true /*is_nonnegative*/)
{
  ET lambda_b(0);
  int m = p.get_m();
  typedef typename Program::B_iterator B_column_iterator;
  B_column_iterator b = p.get_b();
  Infeasibility_certificate_iterator inf = infeasibility_certificate_begin();
  for (int i=0; i<m; ++i, ++b, ++inf)
    lambda_b += ET(*b) * *inf;
  if (lambda_b >= et0)
    return error ("lambda_b >= 0");
  return true;
}

// checks whether lambda^T b < \sum_{j: lambda^TA_j < 0}lambda^TA_j u_j + 
//                             \sum_{j: lambda^TA_j > 0}lambda^TA_j l_j
// ---------------------------------------------------------------------- 
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_infeasible_3 
(const Program& p, const typename std::vector<ET>& lambda_a,
 Tag_false /*is_nonnegative*/)
{
  ET lambda_b(0);
  int m = p.get_m(); 
  typedef typename Program::B_iterator B_column_iterator;
  B_column_iterator b = p.get_b();
  Infeasibility_certificate_iterator inf = infeasibility_certificate_begin();
  for (int i=0; i<m; ++i, ++b, ++inf)
    lambda_b += ET(*b) * *inf;

  ET sum(0); 
  int n = p.get_n();  
  typedef typename Program::L_iterator L_column_iterator;  
  typedef typename Program::U_iterator U_column_iterator;  
  L_column_iterator l = p.get_l();
  U_column_iterator u = p.get_u();  
  for (int j=0; j<n; ++j, ++l, ++u) {
    if (lambda_a[j] < et0) sum += lambda_a[j] * ET(*u);
    if (lambda_a[j] > et0) sum += lambda_a[j] * ET(*l);
  }
  // now compare
  if (lambda_b >= sum)
    return error("lambda^T b >= sum");
  return true;
}
 
// checks whether (Aw)_i <= (>=, ==) 0 if i-th constraint is <= (>=, ==)
// ---------------------------------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_unbounded_1 
(const Program& p)
{
  int m = p.get_m();
  std::vector<ET> aw (m, et0);
  add_Az (p, unboundedness_certificate_begin(), aw);
  typedef typename Program::R_iterator R_column_iterator;
  R_column_iterator r = p.get_r();
  for (int i=0; i<m; ++i, ++r)
    switch (*r) {
    case SMALLER:
      if (aw[i] > et0)
	return error ("i-th constraint <= but (Aw)_i > 0");
      break;
    case EQUAL:
      if (aw[i] != et0)
	return error ("i-th constraint == but (Aw)_i != 0");
      break;
    case LARGER:
      if (aw[i] < et0)
	return error ("i-th constraint >= but (Aw)_i < 0");
      break;
    }
  return true;
}

// checks whether w >= 0
// ---------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_unbounded_2 
(const Program& p, Tag_true /*is_nonnegative*/)
{
  int n = p.get_n();
  Unboundedness_certificate_iterator unb = unboundedness_certificate_begin();
  for (int j=0; j<n; ++j, ++unb)
    if (*unb < et0)
      return error ("some w_j < 0");
  return true;
}

// checks whether w_j >= (<=) 0 if l_j (u_j) is finite
// --------------------------------------------------- 
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_unbounded_2 
(const Program& p, Tag_false /*is_nonnegative*/)
{
  int n = p.get_n();
  typedef typename Program::FL_iterator FL_column_iterator;
  typedef typename Program::FU_iterator FU_column_iterator;
  FL_column_iterator fl = p.get_fl();
  FU_column_iterator fu = p.get_fu();
  Unboundedness_certificate_iterator unb = unboundedness_certificate_begin();
  for (int j=0; j<n; ++j, ++unb, ++fl, ++fu) {
    if (*fl && *unb < et0)
      return error ("some l_j is finite but w_j < 0");
    if (*fu && *unb > et0)
      return error ("some u_j is finite but w_j > 0");
  }  
  return true;
}

// checks whether c^Tw < 0
// -----------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_unbounded_3 
(const Program& p, Tag_true /*is_linear*/)
{
  ET cw(0);
  int n = p.get_n();
  typedef typename Program::C_iterator C_column_iterator;
  C_column_iterator c = p.get_c();
  Unboundedness_certificate_iterator unb = unboundedness_certificate_begin();
  for (int j=0; j<n; ++j, ++c, ++unb) 
    cw += ET(*c) * *unb;
  if (cw >= et0)
    return error("c^Tw >= 0");
  return true;
} 

// checks whether w^TDw = 0 and (c^T + 2x^TD)w < 0
// -----------------------------------------------
template <typename ET>
template <typename Program>
bool Quadratic_program_solution<ET>::is_unbounded_3 
(const Program& p, Tag_false /*is_linear*/)
{
  int n = p.get_n();
  std::vector<ET> two_dw (n, et0);
  add_two_Dz(p, unboundedness_certificate_begin(), two_dw);        // 2Dw
  // check w^TDw = 0
  Unboundedness_certificate_iterator unb = unboundedness_certificate_begin();
  ET wdw(0);
  for (int j=0; j<n; ++j, ++unb) 
    wdw += two_dw[j] * *unb;
  if (wdw != et0)
    return error("w^TDw != 0");
  // check (c^T + 2x^TD)w = c^Tw + 2x^T Dw < 0
  typedef typename Program::C_iterator C_column_iterator;
  C_column_iterator c = p.get_c();
  unb = unboundedness_certificate_begin();
  Variable_numerator_iterator v = variable_numerators_begin();
  ET d = variables_common_denominator();
  if (d <= et0) 
    return error ("common variable denominator is negative");
  ET res(0);
  for (int j=0; j<n; ++j, ++c, ++unb, ++v) 
    res += ET(*c) * *unb * d + two_dw[j] * *v;
  if (res >= 0)
    return error("c^T + 2x^TD)w >= 0");
  return true;
}
  
// computes the product of A with the n-vector given by z
// and adds the result to v; precondition: v has size m
// ------------------------------------------------------
template <typename ET>
template <typename Program, typename Z_iterator >
void Quadratic_program_solution<ET>::add_Az 
(const Program& p, Z_iterator z, typename std::vector<ET>& v)
{
  // iterator types
  typedef typename Program::A_iterator
    A_matrix_iterator;
  typedef typename std::iterator_traits<A_matrix_iterator>::value_type
    A_column_iterator;
  // now compute the product
  A_matrix_iterator a = p.get_a();
  int n = p.get_n();
  for (int j=0; j<n; ++j, ++a, ++z) {
    if (!CGAL::is_zero(*z)) {
      // add A_j * z_j to az
      A_column_iterator a_j = *a;
      for (int i=0; i<p.get_m(); ++i, ++a_j)
  	v[i] += *z * ET(*a_j);
    }
  }
}

// computes the product of 2D with the n-vector given by z
// and adds the result to v; precondition: v has size n
// -------------------------------------------------------
template <typename ET>
template <typename Program, typename Z_iterator >
void Quadratic_program_solution<ET>::add_two_Dz 
(const Program& p, Z_iterator z, typename std::vector<ET>& v)
{
  // we compute the contribution from the enries of D on or below
  // the diagonal rowwise, and the remaining entries columnwise
  typedef typename Program::D_iterator
    D_matrix_iterator;
  typedef typename std::iterator_traits<D_matrix_iterator>::value_type
    D_row_iterator;
  int n = p.get_n();
  // make sure that we only handle the nonzero terms of z
  std::vector<int> z_indices;
  Z_iterator l = z;
  for (int i=0; i<n; ++i, ++l) 
    if (*l != 0) z_indices.push_back(i);

  // the rowwise contribution: on/below the diagonal
  D_matrix_iterator d = p.get_d();
  for (int i=0; i<n; ++i, ++d) {
    // handle row i
    D_row_iterator d_i = *d;  // row i on/below diagonal
    for (std::vector<int>::const_iterator 
	   k = z_indices.begin(); k < z_indices.end(); ++k) {
      int j = *k;
      if (j <= i) v[i] += *(z+j) * ET (*(d_i+j));
    }
  }
  
  // the columnwise contribution: above the diagonal
  d = p.get_d();
  for (std::vector<int>::const_iterator 
	 k = z_indices.begin(); k < z_indices.end(); ++k) {
    int j = *k;
    // handle column j
    D_row_iterator d_j = *(d+j); // column j above diagonal
    for (int i=0; i<j; ++i, ++d_j)
      v[i] += *(z+j) * ET (*d_j);
  }
}


// computes the product of the m-vector given by z with the matrix
// A and adds the result to v; precondition: v has size n 
// ----------------------------------------------------------------
template <typename ET>
template <typename Program, typename Z_iterator >
void  Quadratic_program_solution<ET>::add_zA 
(const Program& p, Z_iterator z, typename std::vector<ET>& v)
{
  typedef typename Program::A_iterator
    A_matrix_iterator;
  typedef typename std::iterator_traits<A_matrix_iterator>::value_type
    A_column_iterator;
  A_matrix_iterator a = p.get_a();
  int n = p.get_n();
  int m = p.get_m();
  // make sure that we only handle the nonzero terms of z
  std::vector<int> z_indices; // indices of nonzero entries
  Z_iterator l = z;
  for (int i=0; i<m; ++i, ++l) 
    if (*l != 0) z_indices.push_back(i);
  // now compute the product columnwise 
  for (int j=0; j<n; ++j, ++a) {
    // add z^T * A_j to v[j] 
    A_column_iterator a_j = *a;
    for (std::vector<int>::const_iterator 
	   k = z_indices.begin(); k < z_indices.end(); ++k)
      v[j] += *(z+*k) * ET(*(a_j+*k));
  }  
}

// adds d c^T to v; precondition: v has length n
// ---------------------------------------------
template <typename ET>
template <typename Program>
void  Quadratic_program_solution<ET>::add_c
(const Program& p, typename std::vector<ET>& v)
{
  int n = p.get_n();
  typedef typename Program::C_iterator C_column_iterator;
  C_column_iterator c = p.get_c();
  ET d = variables_common_denominator();
  for (int j=0; j<n; ++j, ++c)
    v[j] += ET(*c) * d;
}


} //namespace CGAL

#endif // CGAL_QP_SOLUTION_IMPL_H
