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
 CGAL::Tag_true /*is_in_standard_form*/)
{
  // nop (default bounds are nonnegative)
}

// helper for MPS output: BOUNDS
template <typename Fllfuu>
void print_fllfuu_bounds 
(std::ostream& out, const Fllfuu& fllfuu, 
 CGAL::Tag_false /*is_in_standard_form*/)
{
  typename Fllfuu::FL_iterator fl = fllfuu.fl();
  typename Fllfuu::FU_iterator fu = fllfuu.fu();
  typename Fllfuu::L_iterator l = fllfuu.l();
  typename Fllfuu::U_iterator u = fllfuu.u();
  int n = fllfuu.n();
  out << "BOUNDS\n"; 
  for (int j=0; j<n; ++j, ++fl, ++l, ++fu, ++u) {
    if (!*fl || !CGAL::is_zero(*l)) {
      if (*fl)
	out << "  LO  BND  x" << j << "  " << *l << "\n";
      else
	out << "  MI  BND  x" << j << "\n";
    }
    if (*fu)
      out << "  UP  BND  x" << j << "  " << *u << "\n";
  } 
} 

// helper for MPS output: QMATRIX
template <typename D>
void print_d_qmatrix 
(std::ostream& out, const D& d, CGAL::Tag_true /*is_linear*/)
{
  // nop
}

// helper for MPS output: QMATRIX
template <typename D>
void print_d_qmatrix 
(std::ostream& out, const D& d, CGAL::Tag_false /*is_linear*/)
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

template <typename Quadratic_program1, typename Quadratic_program2>
bool are_equal_qp 
(const Quadratic_program1 &qp1, const Quadratic_program2 &qp2)
{
  if (qp1.n() != qp2.n()) return false;
  if (qp1.m() != qp2.m()) return false;
  int n = qp1.n();
  int m = qp1.m();
  // check A
  typename Quadratic_program1::A_iterator a1 = qp1.a();
  typename Quadratic_program2::A_iterator a2 = qp2.a();
  for (int j=0; j<n; ++j, ++a1, ++a2)
    for (int i=0; i<m; ++i) 
      if (*((*a1)+i) != *((*a2)+i)) return false;
  // check b
  typename Quadratic_program1::B_iterator b1 = qp1.b();
  typename Quadratic_program2::B_iterator b2 = qp2.b();
  for (int i=0; i<m; ++i, ++b1, ++b2)
    if (*b1 != *b2) return false;
  // check r
  typename Quadratic_program1::R_iterator r1 = qp1.r();
  typename Quadratic_program2::R_iterator r2 = qp2.r();
  for (int i=0; i<m; ++i, ++r1, ++r2)
    if (*r1 != *r2) return false;
  // check fl
  typename Quadratic_program1::FL_iterator fl1 = qp1.fl();
  typename Quadratic_program2::FL_iterator fl2 = qp2.fl();
  for (int j=0; j<n; ++j, ++fl1, ++fl2)
    if (*fl1 != *fl2) return false;
  // check l
  typename Quadratic_program1::L_iterator l1 = qp1.l();
  typename Quadratic_program2::L_iterator l2 = qp2.l();
  for (int j=0; j<n; ++j, ++l1, ++l2)
    if (*l1 != *l2) return false;
  // check fu
  typename Quadratic_program1::FU_iterator fu1 = qp1.fu();
  typename Quadratic_program2::FU_iterator fu2 = qp2.fu();
  for (int j=0; j<n; ++j, ++fu1, ++fu2)
    if (*fu1 != *fu2) return false;
  // check u
  typename Quadratic_program1::U_iterator u1 = qp1.u();
  typename Quadratic_program2::U_iterator u2 = qp2.u();
  for (int j=0; j<n; ++j, ++u1, ++u2)
    if (*u1 != *u2) return false;
  // check d
  typename Quadratic_program1::D_iterator d1 = qp1.d();
  typename Quadratic_program2::D_iterator d2 = qp2.d();
  for (int j=0; j<n; ++j, ++d1, ++d2)
    for (int i=0; i<n; ++i) 
      if (*((*d1)+i) != *((*d2)+i)) return false;
  // check c
  typename Quadratic_program1::C_iterator c1 = qp1.c();
  typename Quadratic_program2::C_iterator c2 = qp2.c();
  for (int j=0; j<n; ++j, ++c1, ++c2)
    if (*c1 != *c2) return false;
  // check c0
  typename Quadratic_program1::C_entry c01 = qp1.c0();
  typename Quadratic_program2::C_entry c02 = qp2.c0();
  if (c01 != c02) return false;
  return true;
}
CGAL_END_NAMESPACE

#endif // CGAL_QP_FUNCTIONS_IMPL_H
