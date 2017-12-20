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
#ifndef CGAL_QP_FUNCTIONS_IMPL_H
#define CGAL_QP_FUNCTIONS_IMPL_H

#include <CGAL/license/QP_solver.h>


#include <iostream>
#include <fstream>
#include <CGAL/iterator.h>
#include <CGAL/QP_solver/QP_solver.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_solution.h>

namespace CGAL {

namespace QP_functions_detail {
  // test whether the system is of the form A x == b (equations only)
  template <typename R>
  bool is_in_equational_form (const R& r) 
  {	     
    typename R::R_iterator it = r.get_r();
    typename R::R_iterator end = it + r.get_m();
    for (; it < end; ++it)
      if (*it != CGAL::EQUAL) return false;
    return true;
  }

  // test whether the row vectors of A that correpsond to equations 
  // are linearly independent; this is done using type ET. The value
  // type of LinearInequalitySystem must be convertible to ET
  template <class Ar, class ET>
  bool has_linearly_independent_equations 
  (const Ar& ar, const ET& /*dummy*/) {
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
    LP lp (ar.get_n(), ar.get_m(), ar.get_a(), B_iterator(0), ar.get_r(), C_iterator(0));

    //  solver Tags
    typedef QP_solver_impl::QP_tags<
      Tag_true,  // Is_linear
      Tag_true>  // Is_nonnegative
      Tags;

    // solver type
    typedef QP_solver<LP, ET, Tags> Solver;

    // now solve auxiliary LP and compute predicate value
    Solver solver (lp);
    return !solver.diagnostics.redundant_equations;
  }

  // helper for MPS output: BOUNDS
  template <typename P>
  void print_bounds 
  (std::ostream& , const P& ,
   CGAL::Tag_true /*is_nonnegative*/)
  {
    // nop (default bounds are nonnegative)
  }

  // helper for MPS output: BOUNDS
  template <typename P>
  void print_bounds 
  (std::ostream& out, const P& p, 
   CGAL::Tag_false /*is_nonnegative*/)
  {
    typename P::FL_iterator fl = p.get_fl();
    typename P::FU_iterator fu = p.get_fu();
    typename P::L_iterator l = p.get_l();
    typename P::U_iterator u = p.get_u();
    int n = p.get_n();
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

  // helper for MPS output: DMATRIX/QMATRIX
  template <typename P>
  void print_qmatrix 
  (std::ostream& , const P& , 
   CGAL::Tag_true /*is_linear*/)
  {
    // nop
  }

  // helper for MPS output: DMATRIX/QMATRIX
  template <typename P>
  void print_qmatrix 
  (std::ostream& out, const P& p, 
   CGAL::Tag_false /*is_linear*/)
  {
    typename P::D_iterator it = p.get_d();
    int n = p.get_n();
    bool empty_D = true;
    for (int i=0; i<n; ++i, ++it) {
      // handle only entries on/below diagonal
      for (int j=0; j<i+1; ++j)
	if (!CGAL::is_zero(*(*it + j))) {
	  if (empty_D) {
	    // first time we see a nonzero entry
	    out << "QMATRIX\n";
	    empty_D = false;
	  }
	  out << "  x" << i << "  x" << j << "  " << *(*it + j) << "\n";
	  // QMATRIX format prescribes symmetric matrix
	  if (i != j)
	    out << "  x" << j << "  x" << i << "  " << *(*it + j) << "\n";
	}
    }
  }

  // check whether the two qp's have the same data; this is the case iff
  // they agree in n, m, a, b, r, fl, l, fu, u, d, c, c0
  // PRE: qp1, qp2 have the same internal number type
  template <typename Quadratic_program1, typename Quadratic_program2>
  bool are_equal_qp 
  (const Quadratic_program1 &qp1, const Quadratic_program2 &qp2)
  {
    bool return_val = true;
    // check n
    if (qp1.get_n() != qp2.get_n()) {
      std::cerr << "Equality test fails with n: " 
		<< qp1.get_n() << " vs. " << qp2.get_n() << std::endl;
      return false; // wildly wrong, abort now
    }
    // check m
    if (qp1.get_m() != qp2.get_m()) {
      std::cerr << "Equality test fails with m: " 
		<< qp1.get_m() << " vs. " << qp2.get_m() << std::endl;
      return false; // wildly wrong, abort now
    }
    int n = qp1.get_n();
    int m = qp1.get_m();
    // check A
    typename Quadratic_program1::A_iterator a1 = qp1.get_a();
    typename Quadratic_program2::A_iterator a2 = qp2.get_a();
    for (int j=0; j<n; ++j, ++a1, ++a2)
      for (int i=0; i<m; ++i) 
	if (*((*a1)+i) != *((*a2)+i)) {
	  std::cerr << "Equality test fails with A[" 
		    << j << "][" << i << "]: "
		    << *((*a1)+i) << " vs. " <<  *((*a2)+i) << std::endl;
	  return_val = false;
	}
    // check b
    typename Quadratic_program1::B_iterator b1 = qp1.get_b();
    typename Quadratic_program2::B_iterator b2 = qp2.get_b();
    for (int i=0; i<m; ++i, ++b1, ++b2)
      if (*b1 != *b2) {
	std::cerr << "Equality test fails with b[" << i << "]: "
		  << *b1 << " vs. " <<  *b2 << std::endl;	
	return_val = false;
      }
    // check r
    typename Quadratic_program1::R_iterator r1 = qp1.get_r();
    typename Quadratic_program2::R_iterator r2 = qp2.get_r();
    for (int i=0; i<m; ++i, ++r1, ++r2)
      if (*r1 != *r2) {
	std::cerr << "Equality test fails with r[" << i << "]: "
		  << *r1 << " vs. " <<  *r2 << std::endl;	
	return_val = false;
      }
    // check fl, l
    typename Quadratic_program1::FL_iterator fl1 = qp1.get_fl();
    typename Quadratic_program2::FL_iterator fl2 = qp2.get_fl();
    typename Quadratic_program1::L_iterator l1 = qp1.get_l();
    typename Quadratic_program2::L_iterator l2 = qp2.get_l();
    for (int j=0; j<n; ++j, ++fl1, ++fl2, ++l1, ++l2) {
      if (*fl1 != *fl2) {
	std::cerr << "Equality test fails with fl[" << j << "]: "
		  << *fl1 << " vs. " <<  *fl2 << std::endl;	
	return_val = false;
      }
      if ((*fl1 == true) && (*l1 != *l2)) {
	std::cerr << "Equality test fails with l[" << j << "]: "
		  << *l1 << " vs. " <<  *l2 << std::endl;
	return_val = false;
      }
    }
    
    // check fu, u 
    typename Quadratic_program1::FU_iterator fu1 = qp1.get_fu();
    typename Quadratic_program2::FU_iterator fu2 = qp2.get_fu();
    typename Quadratic_program1::U_iterator u1 = qp1.get_u();
    typename Quadratic_program2::U_iterator u2 = qp2.get_u();
    for (int j=0; j<n; ++j, ++fu1, ++fu2, ++u1, ++u2) {
      if (*fu1 != *fu2) {
	std::cerr << "Equality test fails with fu[" << j << "]: "
		  << *fu1 << " vs. " <<  *fu2 << std::endl;
	return_val = false;
      }
      if ((*fu1 == true) && (*u1 != *u2)) {
	std::cerr << "Equality test fails with u[" << j << "]: "
		  << *u1 << " vs. " <<  *u2 << std::endl;
	return_val = false;
      }
    }
    // check d
    typename Quadratic_program1::D_iterator d1 = qp1.get_d();
    typename Quadratic_program2::D_iterator d2 = qp2.get_d();
    for (int i=0; i<n; ++i, ++d1, ++d2)
      for (int j=0; j<i+1; ++j)  // only access entries on/below diagonal
	if (*((*d1)+j) != *((*d2)+j)) {
	  std::cerr << "Equality test fails with D["
		    << i << "][" << j << "]: "
		    << *((*d1)+j) << " vs. " <<  *((*d2)+j) << std::endl; 
	  return_val = false;
	}
    // check c
    typename Quadratic_program1::C_iterator c1 = qp1.get_c();
    typename Quadratic_program2::C_iterator c2 = qp2.get_c();
    for (int j=0; j<n; ++j, ++c1, ++c2)
      if (*c1 != *c2) {
	std::cerr << "Equality test fails with c[" << j << "]: "
		  << *c1 << " vs. " <<  *c2 << std::endl;
	return_val = false;
      }
    // check c0
    typename Quadratic_program1::C_entry c01 = qp1.get_c0();
    typename Quadratic_program2::C_entry c02 = qp2.get_c0();
    if (c01 != c02) {
      std::cerr << "Equality test fails with c0: "
		<< c01 << " vs. " <<  c02 << std::endl;
      return_val = false;
    }
    return return_val;
  }

  template <typename P, typename Is_linear, typename Is_nonnegative>
  void print_program
  (std::ostream& out, const P& p, 
   const std::string& problem_name,
   Is_linear is_linear, 
   Is_nonnegative is_nonnegative)
  {
    // NAME:
    out << "NAME " << problem_name << "\n";

    int n = p.get_n();
    int m = p.get_m();
 
    // ROWS section: 
    typename P::R_iterator r = p.get_r();
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
    typename P::A_iterator a = p.get_a();
    typename P::C_iterator c = p.get_c();
    typedef 
      typename std::iterator_traits<typename P::C_iterator>::value_type IT;
    out << "COLUMNS\n";
    for (int j=0; j<n; ++j, ++c, ++a) {
      // make sure that variable appears here even if it has only
      // zero coefficients
      bool written = false;
      if (!CGAL_NTS is_zero (*c)) {
	out << "  x" << j << "  obj  " << *c << "\n";
	written = true;
      }
      for (int i=0; i<m; ++i) { 
	if (!CGAL_NTS is_zero (*((*a)+i))) {
	  out << "  x" << j << "  c" << i << "  " << *((*a)+i) << "\n";
	  written = true;
	}
      }
      if (!written)
	out << "  x" << j << "  obj  " << IT(0) << "\n";
    }
 
    // RHS section:
    typename P::B_iterator b = p.get_b();
    out << "RHS\n";
    if (!CGAL_NTS is_zero (p.get_c0()))
      out << "  rhs obj " << -p.get_c0() << "\n";
    for (int i=0; i<m; ++i, ++b)  
      if (!CGAL_NTS is_zero (*b))
	out << "  rhs c" << i << "  " << *b << "\n";

    // BOUNDS section:
    QP_functions_detail::print_bounds (out, p, is_nonnegative); 

    // QMATRIX section:
    QP_functions_detail::print_qmatrix (out, p, is_linear);
 
    // output end:
    out << "ENDATA\n";
  }

  template <typename Program, typename ET, 
	    typename Is_linear,typename Is_nonnegative >
  Quadratic_program_solution<ET> solve_program 
  (const Program &p, const ET&, 
   Is_linear, 
   Is_nonnegative,
   const Quadratic_program_options& options)
  { 
    typedef QP_solver<
      Program, ET, 
      QP_solver_impl::QP_tags<Is_linear, Is_nonnegative> >
      Solver;
    const Solver* s = new Solver(p, options);
    Quadratic_program_solution<ET> solution(s);
    if (options.get_auto_validation()) {
      // validate solution
      if (options.get_verbosity() > 0)
	std::cout << "Validating solution...\n";
      if (!solution.solves_program(p, Is_linear(), Is_nonnegative())) {
	// write log file
	std::ofstream out("QP_solver.log");
	out << "Error: Program solution is invalid\n"
	    << "  (error message: " << solution.get_error() << ")\n"
	    << "------------------\n"
	    << "Solution function:\n"
	    << "------------------\n";
	print_solution_function (out, Is_linear(), Is_nonnegative());
	out << "\n"
	    << "--------\n"
	    << "Program:\n" 
	    << "--------\n";
	print_program (out, p, "unsolved", Is_linear(), Is_nonnegative());
	out << "--------\n"
	    << "Options:\n"
	    << "--------\n"
	    << options << std::endl;
	// print warning
	std::cerr 
	  << "Error: Program solution is invalid "
	  << "(see QP_solver.log for details)" << std::endl;
      }
    }
    return solution;
      
  }
}

template <typename QuadraticProgram, typename ET>
Quadratic_program_solution<ET> solve_quadratic_program 
(const QuadraticProgram &qp, const ET& dummy, 
 const Quadratic_program_options& options)
{
  return QP_functions_detail::
    solve_program(qp, dummy, Tag_false(), Tag_false(), options);
}

template <typename NonnegativeQuadraticProgram, typename ET>
Quadratic_program_solution<ET> solve_nonnegative_quadratic_program 
(const NonnegativeQuadraticProgram &qp, const ET& dummy,
 const Quadratic_program_options& options)
{
  return QP_functions_detail::
    solve_program(qp, dummy, Tag_false(), Tag_true(), options);
}

template <typename LinearProgram, typename ET>
Quadratic_program_solution<ET> solve_linear_program 
(const LinearProgram &lp, const ET& dummy,
 const Quadratic_program_options& options)
{
  return QP_functions_detail::
    solve_program(lp, dummy, Tag_true(), Tag_false(), options);
}

template <typename NonnegativeLinearProgram, typename ET>
Quadratic_program_solution<ET> solve_nonnegative_linear_program 
(const NonnegativeLinearProgram &lp, const ET& dummy,
 const Quadratic_program_options& options)
{
  return QP_functions_detail::
    solve_program(lp, dummy, Tag_true(), Tag_true(), options);
}

} //namespace CGAL

#endif // CGAL_QP_FUNCTIONS_IMPL_H
