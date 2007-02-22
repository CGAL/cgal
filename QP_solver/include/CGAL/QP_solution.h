// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>
//               : Bernd Gaertner <gaertner@inf.ethz.ch>
//               : Sven Schoenherr <sven@inf.ethz.ch>
//               : Franz Wessendorp < fransw@inf.ethz.ch>

#ifndef CGAL_QP_SOLUTION_H
#define CGAL_QP_SOLUTION_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/functional.h>
#include <CGAL/function_objects.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/QP_solver/iterator.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

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
  class Optimality_certificate_by_index;
}

// global status type
enum QP_status { QP_UPDATE, QP_INFEASIBLE, QP_UNBOUNDED, QP_OPTIMAL };

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

  typedef typename CGAL::Compose<Quotient_normalizer, 
				 U_Quotient_creator>::Type
  Quotient_creator; // normalized quotient creator (ET, ET)

  typedef  typename CGAL::Bind<Quotient_creator, ET, 2>::Type
  Quotient_maker; // normalized quotient creator (ET, const) 

  typedef std::vector<int> 
  Indices;

  typedef Indices::iterator    
  Index_mutable_iterator;

  typedef Indices::const_iterator    
  Index_const_iterator;

  typedef Transform_diff_const_iterator<int, Identity <int> >
  Original_index_const_iterator; 

  typedef typename QP_solution_detail::Value_by_index<ET> Value_by_index;

  typedef Transform_diff_const_iterator<int, Value_by_index>
  Variable_numerator_iterator;

  typedef  Join_input_iterator_1< Variable_numerator_iterator,
				  Quotient_maker >
  Variable_value_iterator;

  typedef typename QP_solution_detail::Unbounded_direction_by_index<ET> 
  Unbounded_direction_by_index;

  typedef Transform_diff_const_iterator<int, Unbounded_direction_by_index>
  Unbounded_direction_iterator;

  typedef typename QP_solution_detail::Optimality_certificate_by_index<ET> 
  Optimality_certificate_by_index;
  
  typedef Transform_diff_const_iterator<int, Optimality_certificate_by_index>
  Optimality_certificate_numerator_iterator;

  typedef Join_input_iterator_1<Optimality_certificate_numerator_iterator,
				Quotient_maker >
  Optimality_certificate_iterator;

public:

  // virtual access functions to solution that will 
  // be overridden by QP_solver below

  // Solution
  // --------
  virtual Quotient<ET> solution() const = 0;
  virtual QP_status status() const = 0;

  // Variable values
  // ---------------
  virtual ET variable_numerator_value (int i) const = 0; 
  virtual const ET& variables_common_denominator( ) const = 0;

  // value type ET
  virtual Variable_numerator_iterator 
  original_variables_numerator_begin() const = 0;
  virtual Variable_numerator_iterator 
  original_variables_numerator_end() const = 0;

  // value type Quotient<ET>
  virtual Variable_value_iterator original_variable_values_begin() const = 0;
  virtual Variable_value_iterator original_variable_values_end() const = 0;

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
  virtual Unbounded_direction_iterator unbounded_direction_begin() const = 0;
  virtual Unbounded_direction_iterator unbounded_direction_end() const = 0;

  // Optimality
  // ----------
  virtual ET optimality_certificate_numerator(int i) const = 0;

  // value type ET
  virtual Optimality_certificate_numerator_iterator 
  optimality_certificate_numerator_begin() const = 0;
  virtual Optimality_certificate_numerator_iterator 
  optimality_certificate_numerator_end() const = 0;

  // value type Quotient<ET>
  virtual Optimality_certificate_iterator 
  optimality_certificate_begin() const = 0;
  virtual Optimality_certificate_iterator 
  optimality_certificate_end() const = 0;

  // Infeasibility
  // -------------

  // Validity
  // --------
  virtual bool is_valid() const = 0;

  // destruction
  // -----------
  virtual ~QP_solver_base() {}
};


// Quadratic_program_solution class: a handle for QP_solver_base<ET>
// -------------------------------------------------- 
template <class ET>
class Quadratic_program_solution: Handle_for<const QP_solver_base<ET>*> 
{
public:
  // interface types
  typedef typename QP_solver_base<ET>::Variable_value_iterator
  Variable_value_iterator;

  typedef typename QP_solver_base<ET>::Variable_numerator_iterator
  Variable_numerator_iterator;

  typedef typename QP_solver_base<ET>::Index_const_iterator
  Index_iterator;

  typedef typename QP_solver_base<ET>::Unbounded_direction_iterator
  Unbounded_direction_iterator;

  // methods
  Quadratic_program_solution ()
    : Handle_for<const QP_solver_base<ET>*>()
  {
    *(this->ptr()) = 0; // unitialized solution
  }

  Quadratic_program_solution (const QP_solver_base<ET>* s)
    : Handle_for<const QP_solver_base<ET>*>(s)
  {}

  Quotient<ET> solution() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->solution();
  }

  QP_status status() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->status();
  }

  Variable_value_iterator variable_values_begin() const
  {
   CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
   return (*(this->Ptr()))->original_variable_values_begin();
  }

  Variable_value_iterator variable_values_end() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->original_variable_values_end();
  }

  Variable_numerator_iterator variable_numerators_begin() const
  {
   CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
   return (*(this->Ptr()))->original_variables_numerator_begin();
  }

  Variable_numerator_iterator variable_numerators_end() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->original_variables_numerator_end();
  }

  ET variable_numerator_value (int i) const 
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->variable_numerator_value(i);
  }

  const ET& variables_common_denominator() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->variables_common_denominator();
  }

  Unbounded_direction_iterator unbounded_direction_begin() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->unbounded_direction_begin();
  }

  Unbounded_direction_iterator unbounded_direction_end() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->unbounded_direction_end();
  }

  Index_iterator basic_variable_indices_begin() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->basic_original_variable_indices_begin();
  }

  Index_iterator basic_variable_indices_end() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->basic_original_variable_indices_end();
  }

  int number_of_basic_variables() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->number_of_basic_original_variables();
  }

  Index_iterator basic_constraint_indices_begin() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->basic_constraint_indices_begin();
  }

  Index_iterator basic_constraint_indices_end() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->basic_constraint_indices_end();
  }

  int number_of_basic_constraints() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->number_of_basic_constraints();
  }

  bool is_valid() const
  {
    CGAL_qpe_precondition_msg(*(this->ptr()) != 0, "Solution not initialized");
    return (*(this->Ptr()))->is_valid();
  }

  ~Quadratic_program_solution()
  {
    if (!this->is_shared()) delete *(this->ptr());
  }
}; 

// Details
namespace QP_solution_detail {
  // Quotient_normalizer
  // -------------------
  template < typename ET>
  class Quotient_normalizer {
  public:
    typedef CGAL::Quotient<ET> result_type;
    typedef CGAL::Arity_tag<1> Arity;
    
    typedef 
    CGAL::Boolean_tag<CGALi::Is_unique_factorization_domain<ET>::value> 
    Has_gcd;
    
    typedef 
    CGAL::Boolean_tag<CGALi::Is_integral_domain<ET>::value> 
    Has_exact_division;

    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_true /*has_gcd*/,
     Tag_true /*has_exact_division*/) const
    {
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
  class Value_by_index : public std::unary_function< int, ET>
  {
  public:
    typedef QP_solver_base<ET> QP;
    typedef ET result_type;

    Value_by_index(const QP* solver)
      : s (solver)
    {}

    // returns value * denominator 
    result_type operator () ( int i) const
    {
      return s->variable_numerator_value(i);
    }
    
    const QP* s;
  };

  // Unbounded_direction_by_index
  // ----------------------------
  template < typename ET>
  class Unbounded_direction_by_index : public std::unary_function< int, ET>
  {
  public:
    typedef QP_solver_base<ET> QP;
    typedef ET result_type;

    Unbounded_direction_by_index(const QP* solver)
      : s (solver)
    {}

    result_type operator () ( int i) const
    {
      return s->unbounded_direction_value(i);
    }
      
    const QP* s;
  };

  // Optimality_certificate_by_index
  // -------------------------------
  template < typename ET>
  class Optimality_certificate_by_index : public std::unary_function< int, ET>
  {
  public:
    typedef QP_solver_base<ET> QP;
    typedef ET result_type;

    Optimality_certificate_by_index(const QP* solver)
      : s (solver)
    {}

    result_type operator () ( int i) const
    {
      return s->optimality_certificate_numerator(i);
    }
      
    const QP* s;
  };
}
CGAL_END_NAMESPACE

#endif// CGAL_QP_SOLUTION_H
