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
// 
//
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>

#ifndef CGAL_QP_FUNCTIONS_H
#define CGAL_QP_FUNCTIONS_H

#include <CGAL/license/QP_solver.h>


#include <iostream>
#include <string>
#include <CGAL/QP_options.h>
#include <CGAL/QP_solution.h>

namespace CGAL {

namespace QP_functions_detail { 
  // internal routine; writes P to out in MPS format
  // Is_linear == Tag_true / Tag_false:
  //    p is treated as LinearProgram / QuadraticProgram
  // Is_nonnegative == Tag_true / Tag_false
  //    p is treated as Nonnegative / Arbitrary  
  // the dmatrix parameter specificies whether the quadratic matrix (if any)
  // is written in DMATRIX format (no multiplication by two, good for
  // cross-checking output, or in QMATRIX format (good for using other
  // solvers like CPLEX)
  template <typename P, typename Is_linear, typename Is_nonnegative>
  void print_program
  (std::ostream& out, 
   const P &p,
   const std::string& problem_name,
   Is_linear is_linear, 
   Is_nonnegative is_nonnegative); 

  // internal routine: solves a program, depending on the tags
  template <typename Program, typename ET, 
	    typename Is_linear,typename Is_nonnegative >
  Quadratic_program_solution<ET> solve_program 
  (const Program &p, const ET&, 
   Is_linear is_linear, 
   Is_nonnegative is_nonnegative,
   const Quadratic_program_options& options = Quadratic_program_options()); 

  // internal routines: prints name of solution function
  inline void print_solution_function 
  (std::ostream& out, 
   Tag_true /*is_linear*/, Tag_true /*is_nonnegative*/)
  {
    out << "solve_nonnegative_linear_program"; 
  }
  inline void print_solution_function 
  (std::ostream& out, 
   Tag_false /*is_linear*/, Tag_true /*is_nonnegative*/)
  {
    out << "solve_nonnegative_quadratic_program"; 
  }
  inline void print_solution_function 
  (std::ostream& out, 
   Tag_true /*is_linear*/, Tag_false /*is_nonnegative*/)
  {
    out << "solve_linear_program"; 
  }
  inline void print_solution_function 
  (std::ostream& out, 
   Tag_false /*is_linear*/, Tag_false /*is_nonnegative*/)
  {
    out << "solve_quadratic_program"; 
  }
				
  // internal routine:
  // test whether the system is of the form A x == b (equations only)
  template <typename R>
  bool is_in_equational_form (const R& r);

  // internal routine:
  // test whether the row vectors of A that correpsond to equations 
  // are linearly independent; this is done using type ET. The value
  // type of LinearInequalitySystem must be convertible to ET
  template <class Ar, class ET>
  bool has_linearly_independent_equations 
  (const Ar& ar, const ET& dummy);
}

template <typename QuadraticProgram>
void print_quadratic_program 
(std::ostream& out, const QuadraticProgram &qp, 
 const std::string& problem_name = std::string("MY_MPS"))
// writes qp to out in MPS format
{
  QP_functions_detail::print_program 
    (out, qp, problem_name, CGAL::Tag_false(), CGAL::Tag_false());
}

template <typename LinearProgram>
void print_linear_program 
(std::ostream& out, const LinearProgram &lp, 
 const std::string& problem_name = std::string("MY_MPS"))
// writes lp to out in MPS format
{
  QP_functions_detail::print_program 
    (out, lp, problem_name, CGAL::Tag_true(), CGAL::Tag_false());
}

template <typename NonnegativeQuadraticProgram>
void print_nonnegative_quadratic_program 
(std::ostream& out, const NonnegativeQuadraticProgram &qp,
 const std::string& problem_name = std::string("MY_MPS"))
// writes qp to out in MPS format
{
  QP_functions_detail::print_program 
    (out, qp, problem_name, CGAL::Tag_false(), CGAL::Tag_true());
}

template <typename NonnegativeLinearProgram>
void print_nonnegative_linear_program 
(std::ostream& out, const NonnegativeLinearProgram &lp,
 const std::string& problem_name = std::string("MY_MPS"))
// writes lp to out in MPS format
{
  QP_functions_detail::print_program 
    (out, lp, problem_name, CGAL::Tag_true(), CGAL::Tag_true());
}

template <typename QuadraticProgram, typename ET>
Quadratic_program_solution<ET> solve_quadratic_program 
(const QuadraticProgram &qp, const ET&,
const Quadratic_program_options& options = Quadratic_program_options());

template <typename NonnegativeQuadraticProgram, typename ET>
Quadratic_program_solution<ET> solve_nonnegative_quadratic_program 
(const NonnegativeQuadraticProgram &qp, const ET&,
const Quadratic_program_options& options = Quadratic_program_options());

template <typename LinearProgram, typename ET>
Quadratic_program_solution<ET> solve_linear_program 
(const LinearProgram &lp, const ET&, 
const Quadratic_program_options& options = Quadratic_program_options());

template <typename NonnegativeLinearProgram, typename ET>
Quadratic_program_solution<ET> solve_nonnegative_linear_program 
(const NonnegativeLinearProgram &lp, const ET&,
const Quadratic_program_options& options = Quadratic_program_options());

} //namespace CGAL

#include <CGAL/QP_solver/QP_functions_impl.h>

#endif // CGAL_QP_FUNCTIONS_H
