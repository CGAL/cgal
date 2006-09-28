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
// $URL: svn+ssh://gaertner@scm.gforge.inria.fr/svn/cgal/trunk/QP_solver/include/CGAL/QP_solution.h $
// $Id: QP_solution.h 34216 2006-09-13 21:10:11Z gaertner $
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

  typedef Join_input_iterator_1< Original_index_const_iterator,Value_by_index >
  Variable_numerator_iterator;

  typedef  Join_input_iterator_1< Variable_numerator_iterator,
				  Quotient_maker >
  Variable_value_iterator;

public:

  // virtual access functions to solution that will 
  // be overridden by QP_solver below
  virtual Quotient<ET> solution() const = 0;
  virtual QP_status status() const = 0;
  virtual ET variable_value (int i) const = 0; 
  virtual Variable_value_iterator original_variable_values_begin() const = 0;
  virtual Variable_value_iterator original_variable_values_end() const = 0;
  virtual Index_const_iterator 
  basic_original_variable_indices_begin() const = 0;
  virtual Index_const_iterator 
  basic_original_variable_indices_end() const = 0;
  virtual int number_of_basic_original_variables() const = 0;
  virtual Index_const_iterator basic_constraint_indices_begin() const = 0;
  virtual Index_const_iterator basic_constraint_indices_end() const = 0;
  virtual int number_of_basic_constraints() const = 0;
  virtual bool is_valid() const = 0;

  // destruction (overridden by QP_solver)
  virtual ~QP_solver_base() {}
};


// QP_solution class: a handle for QP_solver_base<ET>
// -------------------------------------------------- 
template <class ET>
class QP_solution: Handle_for<const QP_solver_base<ET>*> 
{
private:
  ET variable_value (int i) const
  {
    // gives just the numerator; rename!
    return (*(this->Ptr()))->variable_value(i);
  }
public:
  // interface types
  typedef typename QP_solver_base<ET>::Variable_value_iterator
  Variable_value_iterator;

  typedef typename QP_solver_base<ET>::Index_const_iterator
  Index_iterator;

  // methods
  QP_solution ()
    : Handle_for<const QP_solver_base<ET>*>()
  {
    *(this->ptr()) = 0; // unitialized solution
  }

  QP_solution (const QP_solver_base<ET>* s)
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

  ~QP_solution()
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
    typedef typename ET::Has_gcd Has_gcd; 
    typedef typename ET::Has_exact_division Has_exact_division;

    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_true has_gcd,
     Tag_true has_exact_division) const
    {
      ET gcd = CGAL::gcd (q.numerator(), q.denominator());
      return CGAL::Quotient<ET> 
	(CGAL::exact_division (q.numerator(), gcd),
	 CGAL::exact_division (q.denominator(), gcd));
    }  

    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_true has_gcd,
     Tag_false has_exact_division) const
    {
      return q;
    }
  
    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_false has_gcd,
     Tag_true has_exact_division) const
    {
      return q;
    }

    CGAL::Quotient<ET> normalize 
    (const CGAL::Quotient<ET>& q, 
     Tag_false has_gcd,
     Tag_false has_exact_division) const
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
      return s->variable_value(i);
    }
    const QP* s;
  };

}
CGAL_END_NAMESPACE

#endif// CGAL_QP_SOLUTION_H
