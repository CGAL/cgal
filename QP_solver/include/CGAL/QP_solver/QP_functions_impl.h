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
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>
#ifndef CGAL_QP_FUNCTIONS_IMPL_H
#define CGAL_QP_FUNCTIONS_IMPL_H

#include <CGAL/iterator.h>
#include <CGAL/QP_solver.h>
#include <CGAL/QP_models.h>

CGAL_BEGIN_NAMESPACE

template <typename R>
bool is_in_equational_form (const R& r) 
{	     
  typename R::R_iterator it = r.r();
  typename R::R_iterator end = it + r.m();
  for (; it < end; ++it)
    if (*it != CGAL::EQUAL) return false;
  return true;
}

template <class Ar, class ET>
bool has_linearly_independent_equations 
(const Ar& ar, const ET& dummy) {
  // we solve the following auxiliary LP, using exact type ET:
  // --------
  // min 0
  // A x r  0
  //   x >= 0
  // --------
  // Then A has linearly independent equations if and only if all 
  // artificials have left the basis after phase I; the QP_solver 
  // diagnostics tells us this
  //
  // auxiliary LP type	
  typedef typename 
    std::iterator_traits<typename Ar::C_iterator>::value_type C_value;
  typedef typename 
    std::iterator_traits<typename Ar::B_iterator>::value_type B_value;
  typedef Const_oneset_iterator <C_value>  C_iterator;
  typedef Const_oneset_iterator <B_value>  B_iterator;
  typedef Nonnegative_linear_program_from_iterators
    <typename Ar::A_iterator, B_iterator, 
    typename Ar::R_iterator, C_iterator> LP;

  //  auxiliary LP
  LP lp (ar.n(), ar.m(), ar.a(), B_iterator(0), ar.r(), C_iterator(0));

  //  solver Tags
  typedef QP_solver_impl::QP_tags<
    Tag_true,  // Is_linear
    Tag_true>  // Is_in_standard_form
  Tags;

  // solver type
  typedef QP_solver<LP, ET, Tags> Solver;

  // now solve auxiliary LP and compute predicate value
  Solver solver (lp);
  return !solver.diagnostics.redundant_equations;
}

// helper for MPS output: BOUNDS
template <typename Fllfuu>
void print_fllfuu_bounds 
(std::ostream& out, const Fllfuu& fllfuu,
 CGAL::Tag_true is_in_standard_form)
{
  // nop (default bounds are nonnegative)
}

// helper for MPS output: BOUNDS
template <typename Fllfuu>
void print_fllfuu_bounds 
(std::ostream& out, const Fllfuu& fllfuu, 
 CGAL::Tag_false is_in_standard_form)
{
  typename Fllfuu::FL_iterator fl = fllfuu.fl();
  typename Fllfuu::FU_iterator fu = fllfuu.fu();
  typename Fllfuu::L_iterator l = fllfuu.l();
  typename Fllfuu::U_iterator u = fllfuu.u();
  int n = fllfuu.n();
  out << "BOUNDS\n"; 
  for (int j=0; j<n; ++j, ++fl, ++l, ++fu, ++u) {
    if (!*fl || !CGAL::is_zero(*l)) 
      if (*fl)
	out << "  LO  BND  x" << j << "  " << *l << "\n";
      else
	out << "  MI  BND  x" << j << "\n";
    if (*fu)
      out << "  UP  BND  x" << j << "  " << *u << "\n";
  } 
} 

// helper for MPS output: QMATRIX
template <typename D>
void print_d_qmatrix 
(std::ostream& out, const D& d, CGAL::Tag_true is_linear)
{
  // nop
}

// helper for MPS output: QMATRIX
template <typename D>
void print_d_qmatrix 
(std::ostream& out, const D& d, CGAL::Tag_false is_linear)
{
  typename D::D_iterator it = d.d();
  int n = d.n();
  bool empty_D = true;
  for (int i=0; i<n; ++i, ++it) {
    for (int j=0; j<n; ++j)
      if (!CGAL::is_zero((*it)[j])) {
	if (empty_D) {
	  // first time we see a nonzero entry
	  out << "QMATRIX\n";
	  empty_D = false;
	}
	out << "  x" << i << "  x" << j << "  " << 2*(*it)[j] << "\n";
      }
  }
}

template <typename Abrc, 
	  typename Is_linear, typename Is_in_standard_form>
void print_linear_inequality_system
(std::ostream& out, 
 const Abrc &abrc,
 const std::string& problem_name,
 Is_linear is_linear, 
 Is_in_standard_form is_in_standard_form)
{
  // NAME:
  out << "NAME " << problem_name << "\n";

  int n = abrc.n();
  int m = abrc.m();
 
  // ROWS section: 
  typename Abrc::R_iterator r = abrc.r();
  out << "ROWS\n"
      << "  N obj\n";                       // for the objective function
  for (int i=0; i<m; ++i, ++r) {
    if (*r == CGAL::SMALLER)
      out << "  L";
    else if (*r == CGAL::EQUAL)
      out << "  E";
    else if (*r == CGAL::LARGER)
      out << "  G";
    else
      CGAL_qpe_assertion_msg(false, "incorrect row-type");
    out << " c" << i << "\n";               // row name is CI 
  }

  // COLUMNS section:
  typename Abrc::A_iterator a = abrc.a();
  typename Abrc::C_iterator c = abrc.c();
  out << "COLUMNS\n";
  for (int j=0; j<n; ++j, ++c, ++a) {
    if (!CGAL_NTS is_zero (*c))
      out << "  x" << j << "  obj  " << *c << "\n";
    for (int i=0; i<m; ++i) { 
      if (!CGAL_NTS is_zero ((*a)[i]))
	out << "  x" << j << "  c" << i << "  " << (*a)[i] << "\n";
    }
  }
 
  // RHS section:
  typename Abrc::B_iterator b = abrc.b();
  out << "RHS\n";
  if (!CGAL_NTS is_zero (abrc.c0()))
    out << "  rhs obj   -" << abrc.c0() << "\n";
  for (int i=0; i<m; ++i, ++b)  
    if (!CGAL_NTS is_zero (*b))
      out << "  rhs c" << i << "  " << *b << "\n";

  // BOUNDS section:
  print_fllfuu_bounds (out, abrc, is_in_standard_form); 

  // QMATRIX section:
  print_d_qmatrix (out, abrc, is_linear);
 
  // output end:
  out << "ENDATA\n";
}


CGAL_END_NAMESPACE

#endif // CGAL_QP_FUNCTIONS_IMPL_H
