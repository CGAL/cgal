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
// $URL: svn+ssh://gaertner@scm.gforge.inria.fr/svn/cgal/trunk/QP_solver/include/CGAL/QP_functions.h $
// $Id: QP_functions.h 33922 2006-09-05 12:32:25Z gaertner $
// 
//
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>

#ifndef CGAL_QP_FUNCTIONS_H
#define CGAL_QP_FUNCTIONS_H

#include <iostream>
#include <string>
#include <CGAL/iterator.h>
#include <CGAL/QP_solver.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_solution.h>

CGAL_BEGIN_NAMESPACE

template <typename R>
bool is_in_equational_form (const R& r);
// test whether the system is of the form A x == b (equations only)

template <typename Ar, class ET>
bool has_linearly_independent_equations 
(const Ar& r, const ET& dummy);
// test whether the row vectors of A that correpsond to equations 
// are linearly independent; this is done using type ET. The value
// type of LinearInequalitySystem must be convertible to ET


template <typename Abrc, typename Is_linear, typename Is_in_standard_form>
void print_linear_inequality_system
(std::ostream& out, 
 const Abrc &abrc,
 const std::string& problem_name,
 Is_linear is_linear, 
 Is_in_standard_form is_in_standard_form);
// internal routine; writes lis to out in MPS format
// Is_linear == Tag_true / Tag_false:
//    lis is treated as LinearProgram / QuadraticProgram
// Is_in_standard_form == Tag_true / Tag_false
//    lis is treated as Nonnegative / Arbitrary  

template <typename QuadraticProgram>
void print_quadratic_program 
(std::ostream& out, const QuadraticProgram &qp,
 const std::string& problem_name = std::string("MY_MPS"))
// writes qp to out in MPS format
{
  print_linear_inequality_system 
    (out, qp, problem_name, CGAL::Tag_false(), CGAL::Tag_false());
}

template <typename LinearProgram>
void print_linear_program 
(std::ostream& out, const LinearProgram &lp,
 const std::string& problem_name = std::string("MY_MPS"))
// writes lp to out in MPS format
{
  print_linear_inequality_system 
    (out, lp, problem_name, CGAL::Tag_true(), CGAL::Tag_false());
}

template <typename NonnegativeQuadraticProgram>
void print_nonnegative_quadratic_program 
(std::ostream& out, const NonnegativeQuadraticProgram &qp,
 const std::string& problem_name = std::string("MY_MPS"))
// writes qp to out in MPS format
{
  print_linear_inequality_system 
    (out, qp, problem_name, CGAL::Tag_false(), CGAL::Tag_true());
}

template <typename NonnegativeLinearProgram>
void print_nonnegative_linear_program 
(std::ostream& out, const NonnegativeLinearProgram &lp,
 const std::string& problem_name = std::string("MY_MPS"))
// writes lp to out in MPS format
{
  print_linear_inequality_system 
    (out, lp,  problem_name, CGAL::Tag_true(), CGAL::Tag_true());
}

template <typename QuadraticProgram, typename ET>
QP_solution<ET> solve_quadratic_program 
(const QuadraticProgram &qp, const ET& dummy)
{
  typedef QP_solver<
    QuadraticProgram, ET, 
    QP_solver_impl::QP_tags<Tag_false, Tag_false, Tag_false, Tag_false> >
    Solver;
  const Solver* s = new Solver(qp);
  return QP_solution<ET>(s);
}

template <typename QuadraticProgram, typename ET>
QP_solution<ET> solve_nonnegative_quadratic_program 
(const QuadraticProgram &qp, const ET& dummy)
{
  typedef QP_solver<
    QuadraticProgram, ET, 
    QP_solver_impl::QP_tags<Tag_false, Tag_false, Tag_false, Tag_true> >
    Solver;
  const Solver* s = new Solver(qp);
  return QP_solution<ET>(s);
}

template <typename QuadraticProgram, typename ET>
QP_solution<ET> solve_linear_program 
(const QuadraticProgram &qp, const ET& dummy)
{
  typedef QP_solver<
    QuadraticProgram, ET, 
    QP_solver_impl::QP_tags<Tag_true, Tag_false, Tag_false, Tag_false> >
    Solver;
  const Solver* s = new Solver(qp);
  return QP_solution<ET>(s);
}

template <typename QuadraticProgram, typename ET>
QP_solution<ET> solve_nonnegative_linear_program 
(const QuadraticProgram &qp, const ET& dummy)
{
  typedef QP_solver<
    QuadraticProgram, ET, 
    QP_solver_impl::QP_tags<Tag_true, Tag_false, Tag_false, Tag_true> >
    Solver;
  const Solver* s = new Solver(qp);
  return QP_solution<ET>(s);
}
CGAL_END_NAMESPACE

#include <CGAL/QP_solver/QP_functions_impl.h>

#endif // CGAL_QP_FUNCTIONS_H
