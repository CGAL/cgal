// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licenseges holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Kaspar Fischer
//               : Bernd Gaertner <gaertner@inf.ethz.ch>
//               : Sven Schoenherr
//               : Franz Wessendorp

#ifndef CGAL_QP_SOLUTION_H
#define CGAL_QP_SOLUTION_H

#include <CGAL/license/QP_solver.h>


#include <iostream>
#include <vector>
#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/function_objects.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/QP_solver/assertions.h>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {

// forward references
template <typename Q, typename ET, typename Tags>
class QP_solver;

namespace QP_solution_detail {
  template <typename ET>
  class Quotient_normalizer;
  
  template <typename ET>
  class Value_by_index;

  template <typename ET>
  class Unbounded_direction_by_index;

  template <typename ET>
  class Lambda_by_index;
}

// global status type
enum Quadratic_program_status 
  { 
    QP_UPDATE, 
    QP_INFEASIBLE, 
    QP_UNBOUNDED, 
    QP_OPTIMAL 
  };

// abstract base class of all QP-solvers
// -------------------------------------
template <class ET>
class QP_solver_base
{
public:
  // types
  typedef  CGAL::Creator_2< ET, ET, Quotient<ET> >
  U_Quotient_creator;  // unnormalized quotient creator ET x ET -> (ET, ET)

  typedef QP_solution_detail::Quotient_normalizer<ET> 
  Quotient_normalizer; // normalizer (ET, ET) -> (ET, ET)

  typedef boost::function1< Quotient<ET>, ET > 
  Quotient_maker;

  typedef std::vector<int> 
  Indices;

  typedef Indices::iterator    
  Index_mutable_iterator;

  typedef Indices::const_iterator    
  Index_const_iterator;

  typedef typename QP_solution_detail::Value_by_index<ET> Value_by_index;

  typedef typename boost::transform_iterator
  <Value_by_index, boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t> >
  Variable_numerator_iterator;

  typedef boost::transform_iterator
  <Quotient_maker, Variable_numerator_iterator>
  Variable_value_iterator;

  typedef typename QP_solution_detail::Unbounded_direction_by_index<ET> 
  Unbounded_direction_by_index;

  typedef boost::transform_iterator
  <Unbounded_direction_by_index, boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t> >
  Unbounded_direction_iterator;

  typedef typename QP_solution_detail::Lambda_by_index<ET> 
  Lambda_by_index;
  
  typedef boost::transform_iterator
  <Lambda_by_index, boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t> >
  Lambda_numerator_iterator;

  typedef boost::transform_iterator
  <Quotient_maker,Lambda_numerator_iterator>
  Lambda_iterator;

public:

  // virtual access functions to solution that will 
  // be overridden by QP_solver below

  // Solution
  // --------
  virtual ET solution_numerator() const = 0;
  virtual ET solution_denominator() const = 0;
  Quotient<ET> solution( ) const
  { 
    // workaround to please Boost 1.33.1: 
    ET n = solution_numerator();
    ET d = solution_denominator();
    return 
      boost::bind 
      (Quotient_normalizer(), boost::bind
       (U_Quotient_creator(), _1, _2))
      (n, d);
      // (solution_numerator(), solution_denominator());
  }
  virtual Quadratic_program_status status() const = 0;
  virtual int iterations() const = 0;

  // Variable values
  // ---------------
  virtual ET variable_numerator_value (int i) const = 0; 
  virtual const ET& variables_common_denominator( ) const = 0;
  virtual int number_of_variables() const = 0;

  // value type ET
  Variable_numerator_iterator
  original_variables_numerator_begin( ) const
  { return Variable_numerator_iterator 
      (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0), 
       Value_by_index(this));}
				  
    
  Variable_numerator_iterator
  original_variables_numerator_end  ( ) const
  { return Variable_numerator_iterator 
      (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(number_of_variables()) , 
       Value_by_index(this));} 

  // value type Quotient<ET>   
  Variable_value_iterator
  original_variable_values_begin( ) const
  { return Variable_value_iterator
      (original_variables_numerator_begin(),
       boost::bind 
       (boost::bind 
	(Quotient_normalizer(), boost::bind
	 (U_Quotient_creator(), _1, _2)), _1, variables_common_denominator()));
  }
    
  Variable_value_iterator
  original_variable_values_end  ( ) const
  { return Variable_value_iterator
      (original_variables_numerator_end(),
       boost::bind 
       (boost::bind 
	(Quotient_normalizer(), boost::bind
	 (U_Quotient_creator(), _1, _2)), _1, variables_common_denominator()));
  }
    
  // Basic variables and constraints
  // -------------------------------
  virtual Index_const_iterator 
  basic_original_variable_indices_begin() const = 0;
  virtual Index_const_iterator 
  basic_original_variable_indices_end() const = 0;
  virtual int number_of_basic_original_variables() const = 0;
  virtual Index_const_iterator 
  basic_constraint_indices_begin() const = 0;
  virtual Index_const_iterator 
  basic_constraint_indices_end() const = 0;
  virtual int number_of_basic_constraints() const = 0;

  // Unboundedness
  // -------------
  virtual ET unbounded_direction_value(int i) const = 0;

  Unbounded_direction_iterator unbounded_direction_begin() const 
  { return Unbounded_direction_iterator 
      (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0), 
       Unbounded_direction_by_index(this));}

  // Returns the past-the-end iterator corresponding to
  // unbounded_direction_begin().
  Unbounded_direction_iterator unbounded_direction_end() const
  { return Unbounded_direction_iterator 
      (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(number_of_variables()), 
       Unbounded_direction_by_index(this));}


  // Optimality
  // ----------
  virtual ET lambda_numerator(int i) const = 0;
  virtual int number_of_constraints() const = 0;

  // value type ET
  Lambda_numerator_iterator 
  lambda_numerator_begin() const 
  { return Lambda_numerator_iterator 
      (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0), 
       Lambda_by_index(this));}

  Lambda_numerator_iterator 
  lambda_numerator_end() const
  { return Lambda_numerator_iterator 
      (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(number_of_constraints()), 
       Lambda_by_index(this));}

  // value type Quotient<ET>
  Lambda_iterator
  lambda_begin() const
  {
    return Lambda_iterator
     (lambda_numerator_begin(),
      boost::bind 
      (boost::bind 
       (Quotient_normalizer(), boost::bind
	(U_Quotient_creator(), _1, _2)), _1, variables_common_denominator()));
  }

  Lambda_iterator
  lambda_end() const
  {
    return Lambda_iterator
     (lambda_numerator_end(),
      boost::bind 
      (boost::bind 
       (Quotient_normalizer(), boost::bind
	(U_Quotient_creator(), _1, _2)), _1, variables_common_denominator()));
  }

  // destruction
  // -----------
  virtual ~QP_solver_base() {}
};


// Quadratic_program_solution class: a handle for QP_solver_base<ET>
// ----------------------------------------------------------------- 
template <class ET_>
class Quadratic_program_solution: Handle_for<const QP_solver_base<ET_>*> 
{
public:
  typedef ET_ ET;
  // interface types
  // ===============

  // variable values / indices
  // -------------------------
  typedef typename QP_solver_base<ET>::Variable_value_iterator
  Variable_value_iterator;

  typedef typename QP_solver_base<ET>::Variable_numerator_iterator
  Variable_numerator_iterator;

  typedef typename QP_solver_base<ET>::Index_const_iterator
  Index_iterator;

  // certificates
  // ------------
  typedef typename QP_solver_base<ET>::Unbounded_direction_iterator
  Unboundedness_certificate_iterator;
  
  typedef 
  typename QP_solver_base<ET>::Lambda_numerator_iterator
  Optimality_certificate_numerator_iterator;
  
  typedef typename QP_solver_base<ET>::Lambda_iterator
  Optimality_certificate_iterator;

  typedef typename QP_solver_base<ET>::Lambda_numerator_iterator
  Infeasibility_certificate_iterator;

  // methods
  // -------
  Quadratic_program_solution ()
    : Handle_for<const QP_solver_base<ET>*>(), et0(0)
  {
    *(this->ptr()) = 0; // unitialized solution
  }

  Quadratic_program_solution (const QP_solver_base<ET>* s)
    : Handle_for<const QP_solver_base<ET>*>(s), et0(0)
  {}

  Quadratic_program_solution& 
  operator= (const Quadratic_program_solution& sol)
  {
    if (this != &sol) {
      // delete the old solver if necessary
      if (!this->is_shared()) delete *(this->ptr());
      this->Handle_for<const QP_solver_base<ET>*>::operator=(sol);
    }
    return *this;
  }

  ~Quadratic_program_solution()
  {
    if (!this->is_shared()) delete *(this->ptr());
  }

private:
  const QP_solver_base<ET>* solver() const 
  {
    return *(this->Ptr());
  }

public:
  bool is_void() const 
  {
    return solver() == 0;
  }

  Quotient<ET> objective_value() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->solution();
  }

  ET objective_value_numerator() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->solution_numerator();
  }

  ET objective_value_denominator() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->solution_denominator();
  }

  Quadratic_program_status status() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->status();
  }

  bool is_optimal() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return status() == QP_OPTIMAL;
  }

  bool is_infeasible() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return status() == QP_INFEASIBLE;
  }

  bool is_unbounded() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return status() == QP_UNBOUNDED;
  }

  int number_of_iterations() const
  { 
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->iterations();
  }

  Variable_value_iterator variable_values_begin() const
  {
   CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
   return solver()->original_variable_values_begin();
  }

  Variable_value_iterator variable_values_end() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->original_variable_values_end();
  }

  Variable_numerator_iterator variable_numerators_begin() const
  {
   CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
   return solver()->original_variables_numerator_begin();
  }

  Variable_numerator_iterator variable_numerators_end() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->original_variables_numerator_end();
  }

  const ET& variables_common_denominator() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->variables_common_denominator();
  }

  Index_iterator basic_variable_indices_begin() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->basic_original_variable_indices_begin();
  }

  Index_iterator basic_variable_indices_end() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->basic_original_variable_indices_end();
  }

  int number_of_basic_variables() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->number_of_basic_original_variables();
  }

  Index_iterator basic_constraint_indices_begin() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->basic_constraint_indices_begin();
  }

  Index_iterator basic_constraint_indices_end() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->basic_constraint_indices_end();
  }

  int number_of_basic_constraints() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    return solver()->number_of_basic_constraints();
  }

  Optimality_certificate_numerator_iterator 
  optimality_certificate_numerators_begin() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    CGAL_qpe_assertion(status() == QP_OPTIMAL);
    return solver()->lambda_numerator_begin();
  }

  Optimality_certificate_numerator_iterator 
  optimality_certificate_numerators_end() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    CGAL_qpe_assertion(status() == QP_OPTIMAL);
    return solver()->lambda_numerator_end();
  }

  Optimality_certificate_iterator 
  optimality_certificate_begin() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    CGAL_qpe_assertion(status() == QP_OPTIMAL);
    return solver()->lambda_begin();
  }

  Optimality_certificate_iterator 
  optimality_certificate_end() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    CGAL_qpe_assertion(status() == QP_OPTIMAL);
    return solver()->lambda_end();
  }

  // infeasibility
  // -------------
  Infeasibility_certificate_iterator 
  infeasibility_certificate_begin() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    CGAL_qpe_assertion(status() == QP_INFEASIBLE);
    return solver()->lambda_numerator_begin();
  }

  Infeasibility_certificate_iterator 
  infeasibility_certificate_end() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    CGAL_qpe_assertion(status() == QP_INFEASIBLE);
    return solver()->lambda_numerator_end();
  }
  
  // unboundedness
  // -------------
  Unboundedness_certificate_iterator unboundedness_certificate_begin() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    CGAL_qpe_assertion(status() == QP_UNBOUNDED);
    return solver()->unbounded_direction_begin();
  }

  Unboundedness_certificate_iterator unboundedness_certificate_end() const
  {
    CGAL_qpe_assertion_msg(!is_void(), "Solution not initialized");
    CGAL_qpe_assertion(status() == QP_UNBOUNDED);
    return solver()->unbounded_direction_end();
  }

private:
  ET et0; // 0

  // validity
  // --------
  
  // error message returned by failing validation
  std::string err_msg;

  // the error message is set by the following function
  bool error (const std::string& message) 
  {
    err_msg = message;
    return false;
  }

public:
  bool is_valid() const
  {
    return err_msg.empty();
  }

  const std::string& get_error() const
  {
    return err_msg;
  }

  // these four methods use the certificates to validate the solution
  // of all four program types; in case this fails, the solution becomes
  // invalid (and this explains why the methods are non-const)
  template <class QuadraticProgram>
  bool solves_quadratic_program 
  (const QuadraticProgram& qp)
  { return solves_program(qp, Tag_false(), Tag_false()); }

  template <class LinearProgram>
  bool solves_linear_program 
  (const LinearProgram& lp)
  { return solves_program(lp, Tag_true(), Tag_false()); }

  template <class NonnegativeQuadraticProgram>
  bool solves_nonnegative_quadratic_program 
  (const NonnegativeQuadraticProgram& qp)
  { return solves_program(qp, Tag_false(), Tag_true()); }

  template <class NonnegativeLinearProgram>
  bool solves_nonnegative_linear_program 
  (const NonnegativeLinearProgram& lp)
  { return solves_program(lp, Tag_true(), Tag_true()); }

  // helper used by all four validation methods above; see
  // QP_solver/QP_solution_impl.h for its implementation 
  template <class Program, typename Is_linear, typename Is_nonnegative>
  bool solves_program (const Program& p, 
		       Is_linear is_linear, Is_nonnegative is_nonnegative);

private:
  // helpers used by the previous method 
  template <typename Program>
  bool is_feasible (const Program& p, 
		    typename std::vector<ET>& ax_minus_b,
		    Tag_true /*is_nonnegative*/);
  template <typename Program>
  bool is_feasible (const Program& p, 
		    typename std::vector<ET>& ax_minus_b,
		    Tag_false /*is_nonnegative*/);

  template <typename Program>
  bool is_optimal_1 (const Program& p);

  template <typename Program>
  bool is_optimal_2 (const Program& p, 
		     const typename std::vector<ET>& ax_minus_b);

  template <typename Program>
  bool is_optimal_3 (const Program& p, typename std::vector<ET>& two_Dx,
		     Tag_true /*is_linear*/, Tag_true /*is_nonnegative*/);
  template <typename Program>
  bool is_optimal_3 (const Program& p, typename std::vector<ET>& two_Dx,
		     Tag_false /*is_linear*/, Tag_true /*is_nonnegative*/);
  template <typename Program>
  bool is_optimal_3 (const Program& p, typename std::vector<ET>& two_Dx,
		     Tag_true /*is_linear*/, Tag_false /*is_nonnegative*/);
  template <typename Program>
  bool is_optimal_3 (const Program& p, typename std::vector<ET>& two_Dx,
		     Tag_false /*is_linear*/, Tag_false /*is_nonnegative*/);

  template <typename Program>
  bool is_infeasible_1 (const Program& p);

  template <typename Program>
  bool is_infeasible_2 (const Program& p, 
			typename std::vector<ET>& lambda_a,
			Tag_true /*is_nonnegative*/);
  template <typename Program>
  bool is_infeasible_2 (const Program& p, 
			typename std::vector<ET>& lambda_a,
			Tag_false /*is_nonnegative*/);

  template <typename Program>
  bool is_infeasible_3 (const Program& p, 
			const typename std::vector<ET>& /*lambda_a*/,
			Tag_true /*is_nonnegative*/); 
  template <typename Program>
  bool is_infeasible_3 (const Program& p, 
			const typename std::vector<ET>& lambda_a,
			Tag_false /*is_nonnegative*/);
 
  template <typename Program>
  bool is_unbounded_1 (const Program& p);
 
  template <typename Program>
  bool is_unbounded_2 (const Program& p, Tag_true /*is_nonnegative*/);
  template <typename Program>
  bool is_unbounded_2 (const Program& p, Tag_false /*is_nonnegative*/);

  template <typename Program>
  bool is_unbounded_3 (const Program& p, Tag_true /*is_linear*/);
  template <typename Program>
  bool is_unbounded_3 (const Program& p, Tag_false /*is_linear*/);

  template <typename Program>
  bool is_value_correct 
  (const Program& p, typename std::vector<ET>& /*two_Dx*/, 
   Tag_true /*is_linear*/); 
  
  template <typename Program>
  bool is_value_correct 
  (const Program& p, typename std::vector<ET>& two_Dx,
   Tag_false /*is_linear*/); 

  template <typename Program>
  bool are_constraints_feasible 
  (const Program& p, typename std::vector<ET>& ax);

  template <typename Program>
  bool are_bounds_feasible (const Program& p,  Tag_true /*is_nonnegative*/);
  template <typename Program>
  bool are_bounds_feasible (const Program& p,  Tag_false /*is_nonnegative*/);
  
  template <typename Program, typename Z_iterator >
  void add_Az 
  (const Program& p, Z_iterator z, typename std::vector<ET>& v);

  template <typename Program, typename Z_iterator >
  void add_two_Dz 
  (const Program& p, Z_iterator z, typename std::vector<ET>& v);
  
  template <typename Program, typename Z_iterator >
  void add_zA 
  (const Program& p, Z_iterator z, typename std::vector<ET>& v);

  template <typename Program>
  void add_c
  (const Program& p, typename std::vector<ET>& v);

}; 

// output
template <typename ET>
std::ostream& operator<<
  (std::ostream& o, const Quadratic_program_solution<ET>& s)
{
  o << "status:          ";
  switch (s.status()) {
  case QP_INFEASIBLE:
    return o << "INFEASIBLE\n";
  case QP_UNBOUNDED:
    return o << "UNBOUNDED\n";
  case QP_OPTIMAL:
    o << "OPTIMAL\n";
    break;
  default:
    CGAL_qpe_assertion(false);
  }
  o << "objective value: " << s.objective_value() << "\n";
  o << "variable values:\n";
  int j=0;
  for ( typename Quadratic_program_solution<ET>::Variable_value_iterator 
	  it = s.variable_values_begin(); 
	it < s.variable_values_end(); ++it, ++j)
    o << "  " << j << ": " << *it << "\n";
  return o; 
}

// Details
namespace QP_solution_detail {
  // Quotient_normalizer
  // -------------------
  template < typename ET>
  class Quotient_normalizer {
  public:
    typedef CGAL::Quotient<ET> result_type;
   
  private:
      typedef CGAL::Algebraic_structure_traits<ET> AST;
      typedef typename AST::Algebraic_category Category; 
    
  public:
      typedef CGAL::Boolean_tag<
      CGAL::is_same_or_derived<CGAL::Unique_factorization_domain_tag,Category>::value> 
      Has_gcd;
    
      typedef CGAL::Boolean_tag<
      CGAL::is_same_or_derived<CGAL::Integral_domain_tag,Category>::value> 
      Has_exact_division;

    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_true /*has_gcd*/,
     Tag_true /*has_exact_division*/) const
    {
      if (CGAL::is_zero (q.numerator()))
	return CGAL::Quotient<ET>(ET(0), ET(1));
      ET gcd = CGAL::gcd (q.numerator(), q.denominator());
      return CGAL::Quotient<ET> 
	(CGAL::integral_division (q.numerator(), gcd),
	 CGAL::integral_division (q.denominator(), gcd));
    }  

    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_true /*has_gcd*/,
     Tag_false /*has_exact_division*/) const
    {
      return q;
    }
  
    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_false /*has_gcd*/,
     Tag_true /*has_exact_division*/) const
    {
      return q;
    }

    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_false /*has_gcd*/,
     Tag_false /*has_exact_division*/) const
    {
      return q;
    }

    CGAL::Quotient<ET> operator() (const CGAL::Quotient<ET>& q) const
    {
      return normalize (q, Has_gcd(), Has_exact_division());
    }
  };

  // Value_by_index
  // --------------
  template < typename ET>
  class Value_by_index : public std::unary_function< std::size_t, ET>
  {
  public:
    typedef QP_solver_base<ET> QP;
    typedef ET result_type;

    Value_by_index(const QP* solver)
      : s (solver)
    {}

    // returns value * denominator 
    result_type operator () ( std::size_t i) const
    {
      return s->variable_numerator_value(static_cast<int>(i));
    }
    
    const QP* s;
  };

  // Unbounded_direction_by_index
  // ----------------------------
  template < typename ET>
  class Unbounded_direction_by_index : public std::unary_function< std::size_t, ET>
  {
  public:
    typedef QP_solver_base<ET> QP;
    typedef ET result_type;

    Unbounded_direction_by_index(const QP* solver)
      : s (solver)
    {}

    result_type operator () ( std::size_t i) const
    {
      return s->unbounded_direction_value(static_cast<int>(i));
    }
      
    const QP* s;
  };

  // Lambda_by_index
  // ---------------
  template < typename ET>
  class Lambda_by_index : public std::unary_function< std::size_t, ET>
  {
  public:
    typedef QP_solver_base<ET> QP;
    typedef ET result_type;

    Lambda_by_index(const QP* solver)
      : s (solver)
    {}

    result_type operator () ( std::size_t i) const
    {
      return s->lambda_numerator(static_cast<int>(i));
    }
      
    const QP* s;
  };
}
} //namespace CGAL

#include <CGAL/QP_solver/QP_solution_impl.h>

#endif// CGAL_QP_SOLUTION_H
