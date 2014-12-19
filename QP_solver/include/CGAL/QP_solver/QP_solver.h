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

#ifndef CGAL_QP_SOLVER_H
#define CGAL_QP_SOLVER_H

#include <CGAL/iterator.h>
#include <CGAL/QP_solver/basic.h>
#include <CGAL/QP_solver/functors.h>
#include <CGAL/QP_options.h>
#include <CGAL/QP_solution.h>
#include <CGAL/QP_solver/QP_basis_inverse.h>
#include <CGAL/QP_solver/QP_pricing_strategy.h>
#include <CGAL/QP_solver/QP_full_exact_pricing.h>
#include <CGAL/QP_solver/QP_partial_exact_pricing.h>
#include <CGAL/QP_solver/QP_full_filtered_pricing.h>
#include <CGAL/QP_solver/QP_partial_filtered_pricing.h>
#include <CGAL/QP_solver/QP_exact_bland_pricing.h>

#include <CGAL/algorithm.h>

#include <CGAL/IO/Verbose_ostream.h>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/iterator/transform_iterator.hpp>

#include <vector>
#include <numeric>
#include <algorithm>

namespace CGAL {

// ==================
// class declarations
// ==================
template < typename Q, typename ET, typename Tags  >
class QP_solver;

template <class ET>
class QP_solution; 

namespace QP_solver_impl {   // namespace for implemenation details
  // --------------
  // Tags generator
  // --------------
  template < typename Linear, 
	     typename Nonnegative >
  struct QP_tags {
    typedef Linear                Is_linear;
    typedef Nonnegative           Is_nonnegative;
  };

  template < class Q, class Is_linear >
  struct D_selector {};

  template <class Q>
  struct D_selector<Q, Tag_false> // quadratic
  {
    typedef typename Q::D_iterator D_iterator;
  };

  template <class Q>
  struct D_selector<Q, Tag_true> // linear
  {
    // dummy type, not used
    typedef int** D_iterator;
  };

  template < class Q, class Is_nonnegative >
  struct Bd_selector {};

  template < class Q >
  struct Bd_selector<Q, Tag_false> // nonstandard form
  {
    typedef typename Q::FL_iterator FL_iterator;
    typedef typename Q::L_iterator L_iterator;
    typedef typename Q::FU_iterator FU_iterator;
    typedef typename Q::U_iterator U_iterator;
  };

  template < class Q >
  struct Bd_selector<Q, Tag_true> // standard form
  {
    // dummy types, not used
    typedef int* FL_iterator;
    typedef int* L_iterator;
    typedef int* FU_iterator;
    typedef int* U_iterator;
  };

  // only allow filtered pricing if NT = double
  template <typename Q, typename ET, typename Tags, typename NT>
  struct Filtered_pricing_strategy_selector
  {
    typedef QP_full_exact_pricing<Q, ET, Tags> FF;
    typedef QP_partial_exact_pricing<Q, ET, Tags> PF;
  };

  template <typename Q, typename ET, typename Tags>
  struct Filtered_pricing_strategy_selector<Q, ET, Tags, double> 
  {
    typedef QP_full_filtered_pricing<Q, ET, Tags> FF;
    typedef QP_partial_filtered_pricing<Q, ET, Tags> PF;
  };

} // end of namespace for implementation details

// ================
// class interfaces
// ================


template < typename Q, typename ET, typename Tags >
class QP_solver : public QP_solver_base<ET> {

public: // public types
  typedef  QP_solver<Q, ET, Tags> Self;
  typedef  QP_solver_base<ET> Base;
  
  // types from the QP
  typedef  typename Q::A_iterator   A_iterator;
  typedef  typename Q::B_iterator   B_iterator;
  typedef  typename Q::C_iterator   C_iterator;
  typedef  CGAL::Comparison_result Row_type;
  typedef  typename Q::R_iterator Row_type_iterator;
  
  // the remaining types might not be present in the qp, so the
  // following selectors generate dummy types for them 
  typedef  typename QP_solver_impl::
  D_selector<Q, typename Tags::Is_linear>::
  D_iterator D_iterator;
  typedef typename QP_solver_impl::
  Bd_selector<Q, typename Tags::Is_nonnegative>::
  L_iterator L_iterator;
  typedef typename QP_solver_impl::
  Bd_selector<Q, typename Tags::Is_nonnegative>::
  U_iterator U_iterator;
  typedef typename QP_solver_impl::
  Bd_selector<Q, typename Tags::Is_nonnegative>::
  FL_iterator FL_iterator;
  typedef typename QP_solver_impl::
  Bd_selector<Q, typename Tags::Is_nonnegative>::
  FU_iterator FU_iterator;

  // types from the Tags
  typedef  typename Tags::Is_linear    Is_linear;
  typedef  typename Tags::Is_nonnegative Is_nonnegative;

  // friends
  template <class Q_, class ET_>
  friend bool has_linearly_independent_equations 
  (const Q_& qp, const ET_& dummy);

private: // private types

  // types of original problem:
  typedef  typename std::iterator_traits<A_iterator>::value_type  A_column;
  typedef  typename std::iterator_traits<D_iterator>::value_type  D_row;
  
  typedef  typename std::iterator_traits<A_column  >::value_type  A_entry;
  typedef  typename std::iterator_traits<B_iterator>::value_type  B_entry;
  typedef  typename std::iterator_traits<C_iterator>::value_type  C_entry;
  typedef  typename std::iterator_traits<D_row     >::value_type  D_entry;
  typedef  typename std::iterator_traits<L_iterator>::value_type  L_entry;
  typedef  typename std::iterator_traits<U_iterator>::value_type  U_entry;
  
  // slack columns:
  //
  // The following two types are used to (conceptually) add to the matrix A
  // additional columns that model the constraints "x_s>=0" for the slack
  // variables x_s.  Of course, we do not store the column (which is just
  // plus/minus a unit vector), but maintain a pair (int,bool): the first
  // entry says in which column the +-1 is and the second entry of the pair
  // says whether it is +1 (false) or -1 (true).
  typedef  std::pair<int,bool>        Slack_column;
  typedef  std::vector<Slack_column>  A_slack;

  // artificial columns
  //
  // Artificial columns that are (conceptually) added to the matrix A are
  // handled exactly like slack columns (see above).
  typedef  std::pair<int,bool>        Art_column;
  typedef  std::vector<Art_column>    A_art;

  // special artificial column:
  //
  // Also for the special artificial variable we (conceptually) add a column
  // to A. This column contains only +-1's (but it may contain several nonzero
  // entries).
  typedef  std::vector<A_entry>       S_art;
  
  // auxiliary objective vector (i.e., the objective vector for phase I):
  typedef  std::vector<C_entry>       C_aux;

public: // export some additional types:
  
  typedef  typename Base::Indices     Indices; 
  typedef  typename Base::Index_mutable_iterator   Index_iterator;
  typedef  typename Base::Index_const_iterator     Index_const_iterator;
 
  // For problems in nonstandard form we also export the following type, which
  // for an original variable will say whether it sits at is lower, upper, at
  // its lower and upper (fixed) bound, or at zero, or whether the variable is
  // basic:
  enum  Bound_index  { LOWER, ZERO, UPPER, FIXED, BASIC };

private:
  typedef  std::vector<Bound_index>    Bound_index_values;
  typedef  typename Bound_index_values::iterator
  Bound_index_value_iterator;
  typedef  typename Bound_index_values::const_iterator
  Bound_index_value_const_iterator;
  
  // values (variables' numerators):
  typedef  std::vector<ET>            Values;
  typedef  typename Values::iterator  Value_iterator;
  typedef  typename Values::const_iterator
  Value_const_iterator;
    
  // access values by basic index functor:
  typedef  CGAL::Value_by_basic_index<Value_const_iterator>
  Value_by_basic_index;

  // access to original problem by basic variable/constraint index:
  typedef  QP_vector_accessor<A_column, false, false >  A_by_index_accessor;
  typedef  boost::transform_iterator
  < A_by_index_accessor,Index_const_iterator >
  A_by_index_iterator;

  // todo kf: following can be removed once we have all these (outdated)
  // accessors removed:
  typedef  QP_vector_accessor< B_iterator, false, false >
  B_by_index_accessor;
  typedef  boost::transform_iterator
  < B_by_index_accessor, Index_const_iterator >
  B_by_index_iterator;

  typedef  QP_vector_accessor< C_iterator, false, false >
  C_by_index_accessor;
  typedef  boost::transform_iterator
  <C_by_index_accessor, Index_const_iterator >
  C_by_index_iterator;

  typedef  QP_matrix_accessor< A_iterator, false, true, false, false>
  A_accessor;
  typedef  boost::function1<typename A_accessor::result_type, int>
  A_row_by_index_accessor;
  typedef  boost::transform_iterator 
  < A_row_by_index_accessor, Index_iterator >
  A_row_by_index_iterator;

  // Access to the matrix D sometimes converts to ET, and 
  // sometimes retruns the original input type
  typedef  QP_matrix_pairwise_accessor< D_iterator, ET >
  D_pairwise_accessor;
  typedef boost::transform_iterator 
  < D_pairwise_accessor, Index_const_iterator>
  D_pairwise_iterator;
  typedef  QP_matrix_pairwise_accessor< D_iterator, D_entry >
  D_pairwise_accessor_input_type;
  typedef  boost::transform_iterator
  < D_pairwise_accessor_input_type, Index_const_iterator >
  D_pairwise_iterator_input_type;

  // access to special artificial column by basic constraint index:
  typedef  QP_vector_accessor< typename S_art::const_iterator, false, false>
  S_by_index_accessor;
  typedef  boost::transform_iterator
  < S_by_index_accessor, Index_iterator >
  S_by_index_iterator;
  
public:
    
  typedef  typename A_slack::const_iterator
  A_slack_iterator;

  typedef  typename A_art::const_iterator
  A_artificial_iterator;
    
  typedef  typename C_aux::const_iterator
  C_auxiliary_iterator;

  typedef typename Base::Variable_numerator_iterator
  Variable_numerator_iterator;

  typedef  Index_const_iterator       Basic_variable_index_iterator;
  typedef  Value_const_iterator       Basic_variable_numerator_iterator;
  typedef  Index_const_iterator       Basic_constraint_index_iterator;
        
  typedef  QP_pricing_strategy<Q, ET, Tags>  Pricing_strategy;

private:
  // compile time tag for symbolic perturbation, should be moved into traits
  // class when symbolic perturbation is to be implemented
  Tag_false                is_perturbed;
    
  // some constants
  const ET                 et0, et1, et2;

  // verbose output streams
  mutable Verbose_ostream  vout;      // used for any  diagnostic output
  mutable Verbose_ostream  vout1;     // used for some diagnostic output
  mutable Verbose_ostream  vout2;     // used for more diagnostic output
  mutable Verbose_ostream  vout3;     // used for full diagnostic output
  mutable Verbose_ostream  vout4;     // used for output of basis inverse
  mutable Verbose_ostream  vout5;     // used for output of validity tests
    
  // pricing strategy
  Pricing_strategy*        strategyP;

  // given QP
  int                      qp_n;      // number of variables
  int                      qp_m;      // number of constraints
    
  // min x^T D x + c^T x + c0
  A_iterator               qp_A;      // constraint matrix
  B_iterator               qp_b;      // right-hand-side vector
  C_iterator               qp_c;      // objective vector
  C_entry                  qp_c0;     // constant term in objective function
  // attention: qp_D represents *twice* the matrix D
  D_iterator               qp_D;      // objective matrix
  Row_type_iterator        qp_r;      // row-types of constraints
  FL_iterator              qp_fl;     // lower bound finiteness vector
  L_iterator               qp_l;      // lower bound vector
  FU_iterator              qp_fu;     // upper bound finiteness vector
  U_iterator               qp_u;      // upper bound vector

  A_slack                  slack_A;   // slack part of constraint matrix

  // auxiliary problem    
  A_art                    art_A;     // artificial part of constraint matrix
  // Note: in phase I there is an
  // additional "fake" column attached
  // to this "matrix", see init_basis()

  S_art                    art_s;     // special artificial column for slacks
  int                      art_s_i;   // art_s_i>=0  -> index of special
  //                artificial column
  // art_s_i==-1 -> no sp. art. col
  // art_s_i==-2 -> sp. art. col removed
  //                after it left basis 
  int                      art_basic; // number of basic artificial variables
  C_aux                    aux_c;     // objective function for phase I
  // initially has the same size as A_art

  Indices                  B_O;       // basis (original variables)
  // Note: the size of B_O is always
  // correct, i.e., equals the number of
  // basic original variables, plus (in
  // phase I) the number of basic
  // artificial variables.

  Indices                  B_S;       // basis (   slack variables)
    
  Indices                  C;         // basic constraints ( C = E+S_N )
  // Note: the size of C is always
  // correct, i.e., corresponds to the
  // size of the (conceptual) set
  // $E\cup S_N$.

  Indices                  S_B;       // nonbasic constraints ( S_B '=' B_S)
    
  QP_basis_inverse<ET,Is_linear>
  inv_M_B;   // inverse of basis matrix

  const ET&                d;         // reference to `inv_M_B.denominator()'
    
  Values                   x_B_O;     // basic variables (original)
  // Note: x_B_O is only enlarged,
  // so its size need not be |B|.

  Values                   x_B_S;     // basic variables (slack)
  Values                   lambda;    // lambda (from KKT conditions)
  Bound_index_values       x_O_v_i;   // bounds value index vector
  // the following vectors are updated
  // with each update in order to avoid
  // evaluating a matrix vector
  // multiplication
  Values                   r_C;       // r_C = A_{C,N_O}x_{N_O}
  // Note: r_C.size() == C.size().

  Values                   r_S_B;     // r_S_B = A_{S_B,N_O}x_{N_O}

  // The following to variables are initialized (if used at all) in
  // transition().  They are not used in case Is_linear or
  // Is_nonnegative is set to Tag_true.
  Values                   r_B_O;     // r_B_O = 2D_{B_O,N_O}x_{N_O}
  Values                   w;         // w = 2D_{O, N_O}x_{N_O}
    
  int                      m_phase;   // phase of the Simplex method
  Quadratic_program_status                   m_status;  // status of last pivot step
  int                      m_pivots;  // number of pivot steps
    
  bool                     is_phaseI; // flag indicating phase I
  bool                     is_phaseII;// flag indicating phase II
  bool                     is_RTS_transition; // flag indicating transition
  // from Ratio Test Step1 to Ratio
  // Test Step2                                           
  const bool               is_LP;     // flag indicating a linear program
  const bool               is_QP;     // flag indicating a quadratic program

  // the following flag indicates whether the program is in equational form
  // AND still has all its equations; this is given in phase I for any
  // program in equational form, but it may change if redundant constraints
  // get removed from the basis. If no_ineq == true, the program is treated
  // in a more efficient manner, since in that case we need no bookkeeping 
  // for basic constraints
  bool                     no_ineq;   
  bool                     has_ineq;  // !no_ineq

  const bool               is_nonnegative; // standard form, from Tag

  // additional variables
  int                      l;         // minimum of 'qp_n+e+1' and 'qp_m'
  // Note: this is an upper bound for
  // the size of the reduced basis in
  // phase I (in phase II, the reduced
  // basis size can be arbitrarily
  // large)
    
  int 		     e;         // number of equality constraints
    
  // Given a variable number i, in_B[i] is -1 iff x_i is not in the current
  // basis.  If the number in_B[i] is >=0, it is the basis heading of x_i.
  Indices                  in_B;      // variable   in basis, -1 if non-basic

  // Given a number i in {0,...,qp_m-1} of a constraint, 
  Indices                  in_C;      // constraint in basis, -1 if non-basic
  // Note: in_C is only maintained if
  // there are inequality constraints.

  Values                   b_C;       // exact version of `qp_b'
  // restricted to basic constraints C
  Values                   minus_c_B; // exact version of `-qp_c'
  // restricted to basic variables B_O
  // Note: minus_c_B is only enlarged,
  // so its size need not be |B|.

  Values                   A_Cj;      // exact version of j-th column of A
  // restricted to basic constraints C
  Values                   two_D_Bj;  // exact version of twice the j-th
  // column of D restricted to B_O
  // Note: tmp_x_2 is only enlarged,
  // so its size need not be |B|.
    
  int                      j;         // index of entering variable `x_j'
    
  int                      i;         // index of leaving variable `x_i'
  ET                       x_i;       // numerator of leaving variable `x_i'
  ET                       q_i;       // corresponding `q_i'
  int                      direction; // indicates whether the current
  // entering variable x_j is increased
  // or decreased
  Bound_index              ratio_test_bound_index;  // indicates for leaving
  // original variables which bound
  // was hit with upper bounding

  ET                       mu;        //   numerator of `t_j'
  ET                       nu;        // denominator of `t_j'
    
  Values                   q_lambda;  // length dependent on C
  Values                   q_x_O;     // used in the ratio test & update
  // Note: q_x_O is only enlarged,
  // so its size need not be |B|.

  Values                   q_x_S;     // 

  Values                   tmp_l;     // temporary vector of size l
  Values                   tmp_x;     // temporary vector of s. >= B_O.size()
  // Note: tmp_x is only enlarged,
  // so its size need not be |B|.

  Values                   tmp_l_2;   // temporary vector of size l
  Values                   tmp_x_2;   // temporary vector of s. >= B_O.size()
  // Note: tmp_x_2 is only enlarged,
  // so its size need not be |B|.
  // Diagnostics
  struct Diagnostics {
    bool redundant_equations;
  };
    
  Diagnostics              diagnostics;

public:

  /*
   * Note: Some member functions below are suffixed with '_'.
   * They are member templates and their declaration is "hidden",
   * because they are also implemented in the class interface.
   * This is a workaround for M$-VC++, which otherwise fails to
   * instantiate them correctly.
   */

  // creation & initialization
  // -------------------------
  // creation
  QP_solver(const Q& qp, 
	    const Quadratic_program_options& options = 
	    Quadratic_program_options());

  virtual ~QP_solver()
  {
    if (strategyP != static_cast<Pricing_strategy*>(0))
      delete strategyP;
  }

	      
private:
  // set-up of QP
  void set( const Q& qp); 
  void set_D (const Q& qp, Tag_true is_linear);
  void set_D (const Q& qp, Tag_false is_linear);
	         
  // set-up of explicit bounds
  void set_explicit_bounds(const Q& qp); 
  void set_explicit_bounds(const Q& qp, Tag_true /*is_nonnegative*/); 
  void set_explicit_bounds(const Q& qp, Tag_false /*is_nonnegative*/);

  // initialization (of phase I)
  void  init( );
    
  // initialization (of phase II)
  /*
    template < class InputIterator >
    void  init( InputIterator  basic_variables_first,
    InputIterator  basic_variables_beyond);
  */
    
  // operations
  // ----------
  // pivot step
  Quadratic_program_status  pivot( )
  { CGAL_qpe_assertion( phase() > 0);
  CGAL_qpe_assertion( phase() < 3);
  pivot_step();
  return status(); }

  // solve QP
  Quadratic_program_status  solve( )
  { CGAL_qpe_assertion( phase() > 0);
  while ( phase() < 3) { pivot_step(); }
  return status(); }

public:

  // access
  // ------
  // access to QP
  int  number_of_variables  ( ) const { return qp_n; }
  int  number_of_constraints( ) const { return qp_m; }
    
  A_iterator  a_begin( ) const { return qp_A;      }
  A_iterator  a_end  ( ) const { return qp_A+qp_n; }
    
  B_iterator  b_begin( ) const { return qp_b;      }
  B_iterator  b_end  ( ) const { return qp_b+qp_m; }
    
  C_iterator  c_begin( ) const { return qp_c;      }
  C_iterator  c_end  ( ) const { return qp_c+qp_n; }

  C_entry     c_0    ( ) const { return qp_c0;}
    
  D_iterator  d_begin( ) const { return qp_D;      }
  D_iterator  d_end  ( ) const { return qp_D+qp_n; }
    
  Row_type_iterator  row_type_begin( ) const { return qp_r;      }
  Row_type_iterator  row_type_end  ( ) const { return qp_r+qp_m; }

  // access to current status
  int     phase     ( ) const { return m_phase;  }
  Quadratic_program_status  status    ( ) const { return m_status; }
  int     iterations( ) const { return m_pivots; }
    
  // access to common denominator
  const ET& variables_common_denominator( ) const 
  { 
    CGAL_qpe_assertion (d > 0);
    return d; 
  }

  // access to current solution
  ET  solution_numerator( ) const;

  // access to current solution
  ET  solution_denominator( ) const { return et2*d*d; }
    
  // access to original variables
  int  number_of_original_variables( ) const { return qp_n; }
    
  // access to slack variables
  int  number_of_slack_variables( ) const { return static_cast<int>(slack_A.size()); }

  // access to artificial variables
  int  number_of_artificial_variables( ) const { return static_cast<int>(art_A.size()); }
    
  C_auxiliary_iterator
  c_auxiliary_value_iterator_begin( ) const { return aux_c.begin(); }
  C_auxiliary_iterator
  c_auxiliary_value_iterator_end( ) const {return aux_c.end(); }

  // access to basic variables
  int  number_of_basic_variables( ) const { return static_cast<int>(B_O.size()+B_S.size()); }
  int  number_of_basic_original_variables( ) const { return static_cast<int>(B_O.size()); }
  int  number_of_basic_slack_variables( ) const { return static_cast<int>(B_S.size()); }

  Basic_variable_index_iterator
  basic_original_variable_indices_begin( ) const { return B_O.begin(); }
  Basic_variable_index_iterator
  basic_original_variable_indices_end  ( ) const { return B_O.end(); }
    
  Basic_variable_numerator_iterator
  basic_original_variables_numerator_begin( ) const { return x_B_O.begin(); }
  Basic_variable_numerator_iterator
  basic_original_variables_numerator_end  ( ) const { return x_B_O.begin()
							+ B_O.size(); }

public: // only the pricing strategies (including user-defined ones
        // need access to this) -- make them friends?

  // access to working variables
  int  number_of_working_variables( ) const { return static_cast<int>(in_B.size()); }
  
  bool is_basic( int j) const
  { 
    CGAL_qpe_assertion(j >= 0);
    CGAL_qpe_assertion(j < number_of_working_variables());
    return (in_B[ j] >= 0);
  }
  
  bool is_original(int j) const
  {
    CGAL_qpe_assertion(j >= 0);
    CGAL_qpe_assertion(j < number_of_working_variables());
    return (j < qp_n);    
  }
    
  bool phaseI( ) const {return is_phaseI;}
  
  bool is_artificial(int k) const;

  int get_l() const;

  // Returns w[j] for an original variable x_j.
  ET w_j_numerator(int j) const
  { 
    CGAL_qpe_assertion((0 <= j) && (j < qp_n) && is_phaseII);
    return w[j];
  }
  
  Bound_index nonbasic_original_variable_bound_index(int i) const
    // Returns on which bound the nonbasic variable x_i is currently
    // sitting:
    //
    // - LOWER: the variable is sitting on its lower bound.
    // - UPPER: the variable is sitting on its upper bound.
    // - FIXED: the variable is sitting on its lower and upper bound.
    // - ZERO: the variable has value zero and is sitting on its lower
    //   bound, its upper bound, or betweeen the two bounds.
    //
    // Note: in the latter case you can call state_of_zero_nonbasic_variable()
    // to find out which bound is active, if any.
  {
    CGAL_assertion(!check_tag(Is_nonnegative()) &&
		   !is_basic(i) && i < qp_n);
    if (x_O_v_i[i] == BASIC) {
      CGAL_qpe_assertion(false);
    }
    return x_O_v_i[i];  
  };
  
  int state_of_zero_nonbasic_variable(int i) const
    // Returns -1 if the original variable x_i equals its lower bound,
    // 0 if it lies strictly between its lower and upper bound, and 1 if
    // it coincides with its upper bound.
    // 
    // See also the documentation of nonbasic_original_variable_bound_index()
    // above.
  {
    CGAL_assertion(!check_tag(Is_nonnegative()) &&
		   !is_basic(i) && i < qp_n && x_O_v_i[i] == ZERO);
    if (*(qp_fl+i) && CGAL::is_zero(*(qp_l+i)))
      return -1;
    if (*(qp_fu+i) && CGAL::is_zero(*(qp_u+i)))
      return 1;
    return 0;
  }
  
private:
  // miscellaneous
  // -------------
  // setting the pricing strategy:
  void  set_pricing_strategy ( Quadratic_program_pricing_strategy strategy);

  // diagnostic output
  void  set_verbosity( int verbose = 0, std::ostream& stream = std::cout);


public:
  // access to indices of basic constraints
  int  number_of_basic_constraints( ) const { return static_cast<int>(C.size()); }

  Basic_constraint_index_iterator
  basic_constraint_indices_begin( ) const { return C.begin(); }
  Basic_constraint_index_iterator
  basic_constraint_indices_end  ( ) const { return C.end(); }

  // helper functions
  template < class RndAccIt1, class RndAccIt2, class NT >  
  NT  mu_j_( int j, RndAccIt1 lambda_it, RndAccIt2 x_it, const NT& dd) const;

  ET  dual_variable( int i)
  {
    for ( int j = 0; j < qp_m; ++j) {
      tmp_x[ j] = inv_M_B.entry( j, i);
    }
    return std::inner_product( tmp_x.begin(), tmp_x.begin()+qp_m,
			       minus_c_B.begin(), et0);
  }

public:
  // public access to compressed lambda (used in filtered base)
  Value_const_iterator get_lambda_begin() const
  {
    return lambda.begin();
  }
  Value_const_iterator get_lambda_end() const
  {
    return lambda.begin() + C.size();
  }
  

private:    

  // private member functions
  // ------------------------
  // initialization
  void  init_basis( );
  void  init_basis__slack_variables( int s_i, Tag_true  has_no_inequalities);
  void  init_basis__slack_variables( int s_i, Tag_false has_no_inequalities);
  void  init_basis__slack_variables( int s_i, bool has_no_inequalities) {
    if (has_no_inequalities)
      init_basis__slack_variables (s_i, Tag_true());
    else 
      init_basis__slack_variables (s_i, Tag_false());
  }

  void  init_basis__constraints    ( int s_i, Tag_true  has_no_inequalities);
  void  init_basis__constraints    ( int s_i, Tag_false has_no_inequalities);
  void  init_basis__constraints    ( int s_i, bool has_no_inequalities) {
    if (has_no_inequalities)
      init_basis__constraints (s_i, Tag_true());
    else 
      init_basis__constraints (s_i, Tag_false());
  }

  void  init_x_O_v_i();
  void  init_r_C(Tag_true  /*is_nonnegative*/);
  void  init_r_C(Tag_false /*is_nonnegative*/);
  void  init_r_S_B(Tag_true  /*is_nonnegative*/);
  void  init_r_S_B(Tag_false /*is_nonnegative*/);

  void  init_r_B_O();
  void  init_w();


  void  init_solution( );
  void  init_solution__b_C( Tag_true  has_no_inequalities);
  void  init_solution__b_C( Tag_false has_no_inequalities); 
  void  init_solution__b_C( bool has_no_inequalities) {
    if (has_no_inequalities)
      init_solution__b_C (Tag_true());
    else
      init_solution__b_C (Tag_false());
  }

  void  init_additional_data_members( );
    
  // function needed for set up of auxiliary problem for symbolic perturbation
  int  signed_leading_exponent( int row);
  // This is a variant of set_up_auxiliary_problem for symbolic perturbation
  // for the perturbed case
  void  set_up_auxiliary_problemI( Tag_true is_perturbed);

  void  set_up_auxiliary_problem();

  // transition (to phase II)
  void  transition( );
  void  transition( Tag_true  is_linear);
  void  transition( Tag_false is_linear);
    
  // pivot step
  void  pivot_step( );

  // pricing
  void  pricing( );

  template < class NT, class It >
  void  mu_j__linear_part_( NT& mu_j, int j, It lambda_it,
			    Tag_true  has_no_inequalities) const;
  template < class NT, class It >
  void  mu_j__linear_part_( NT& mu_j, int j, It lambda_it,
			    Tag_false has_no_inequalities) const;
  template < class NT, class It >
  void  mu_j__linear_part_( NT& mu_j, int j, It lambda_it,
			    bool has_no_inequalities) const {
    if (has_no_inequalities)
      mu_j__linear_part_ (mu_j, j, lambda_it, Tag_true());
    else
      mu_j__linear_part_ (mu_j, j, lambda_it, Tag_false());
  }


//   template < class NT, class It >
//   void  mu_j__quadratic_part_( NT& mu_j, int j, It x_it,
// 			       Tag_true  is_linear) const;
//   template < class NT, class It >
//   void  mu_j__quadratic_part_( NT& mu_j, int j, It x_it,
// 			       Tag_false is_linear) const;
//   template < class NT, class It >
//   void  mu_j__quadratic_part_( NT& mu_j, int j, It x_it,
// 			       Tag_false is_linear,
// 			       Tag_true  is_symmetric) const;
//   template < class NT, class It >
//   void  mu_j__quadratic_part_( NT& mu_j, int j, It x_it,
// 			       Tag_false is_linear,
// 			       Tag_false is_symmetric) const;

  template < class NT, class It >
  void  mu_j__slack_or_artificial_( NT& mu_j, int j, It lambda_it, 
				    const NT& dd,
				    Tag_true  has_no_inequalities) const;
  template < class NT, class It >
  void  mu_j__slack_or_artificial_( NT& mu_j, int j, It lambda_it, 
				    const NT& dd,
				    Tag_false has_no_inequalities) const;

  template < class NT, class It >
  void  mu_j__slack_or_artificial_( NT& mu_j, int j, It lambda_it,
				    const NT& dd,
				    bool has_no_inequalities) const {
    if (has_no_inequalities)
      mu_j__slack_or_artificial_ (mu_j, j, lambda_it, dd, Tag_true());
    else
      mu_j__slack_or_artificial_ (mu_j, j, lambda_it, dd, Tag_false());
  }

  // ratio test
  void  ratio_test_init( );
  void  ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j,
			       Tag_true  has_no_inequalities);
  void  ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j,
			       Tag_false has_no_inequalities);
  void  ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j,
			       bool has_no_inequalities) {
    if (has_no_inequalities) 
      ratio_test_init__A_Cj (A_Cj_it, j, Tag_true());
    else
      ratio_test_init__A_Cj (A_Cj_it, j, Tag_false());
  }

  void  ratio_test_init__2_D_Bj( Value_iterator two_D_Bj_it, int j,
				 Tag_true  is_linear);
  void  ratio_test_init__2_D_Bj( Value_iterator two_D_Bj_it, int j,
				 Tag_false is_linear);
  void  ratio_test_init__2_D_Bj( Value_iterator two_D_Bj_it, int j,
				 Tag_false is_linear,
				 Tag_true  has_no_inequalities);
  void  ratio_test_init__2_D_Bj( Value_iterator two_D_Bj_it, int j,
				 Tag_false is_linear,
				 Tag_false has_no_inequalities);
  void  ratio_test_init__2_D_Bj( Value_iterator two_D_Bj_it, int j,
				 Tag_false is_linear,
				 bool has_no_inequalities) {
    if (has_no_inequalities)
      ratio_test_init__2_D_Bj( two_D_Bj_it, j, is_linear, Tag_true());
    else
      ratio_test_init__2_D_Bj( two_D_Bj_it, j, is_linear, Tag_false());
  }


  void  ratio_test_1( );
  void  ratio_test_1__q_x_O( Tag_true  is_linear);
  void  ratio_test_1__q_x_O( Tag_false is_linear);
  void  ratio_test_1__q_x_S( Tag_true  has_no_inequalities);
  void  ratio_test_1__q_x_S( Tag_false has_no_inequalities);
  void  ratio_test_1__q_x_S( bool has_no_inequalities) {
    if (has_no_inequalities)
      ratio_test_1__q_x_S (Tag_true());
    else
      ratio_test_1__q_x_S (Tag_false());
  }

  void  ratio_test_1__t_min_j(Tag_true  /*is_nonnegative*/);  
  void  ratio_test_1__t_min_j(Tag_false /*is_nonnegative*/);
    
  void  ratio_test_1__t_i( Index_iterator i_it, Index_iterator end_it,
			   Value_iterator x_it, Value_iterator   q_it,
			   Tag_true  no_check);
  void  ratio_test_1__t_i( Index_iterator i_it, Index_iterator end_it,
			   Value_iterator x_it, Value_iterator   q_it,
			   Tag_false  no_check);
    
  // replaces the above two functions
  void  ratio_test_1__t_min_B(Tag_true has_no_inequalities );
  void  ratio_test_1__t_min_B(Tag_false has_no_inequalities ); 
  void  ratio_test_1__t_min_B(bool has_no_inequalities ) {
    if (has_no_inequalities)
      ratio_test_1__t_min_B (Tag_true());
    else
      ratio_test_1__t_min_B (Tag_false());
  }
   
  void  ratio_test_1_B_O__t_i(Index_iterator i_it, Index_iterator end_it,
			      Value_iterator x_it, Value_iterator q_it,
			      Tag_true  /*is_nonnegative*/);
  void  ratio_test_1_B_O__t_i(Index_iterator i_it, Index_iterator end_it,
			      Value_iterator x_it, Value_iterator q_it,
			      Tag_false /*is_nonnegative*/);
  void  ratio_test_1_B_S__t_i(Index_iterator i_it, Index_iterator end_it,
			      Value_iterator x_it, Value_iterator q_it,
			      Tag_true  /*is_nonnegative*/);
  void  ratio_test_1_B_S__t_i(Index_iterator i_it, Index_iterator end_it,
			      Value_iterator x_it, Value_iterator q_it,
			      Tag_false /*is_nonnegative*/);
			     
  void  test_implicit_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
				     int& i_min, ET& d_min, ET& q_min);
  void  test_implicit_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
				     int& i_min, ET& d_min, ET& q_min);
  void  test_explicit_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
				     int& i_min, ET& d_min, ET& q_min);
  void  test_explicit_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
				     int& i_min, ET& d_min, ET& q_min);
  void  test_mixed_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
				  int& i_min, ET& d_min, ET& q_min);
  void  test_mixed_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
				  int& i_min, ET& d_min, ET& q_min);    
                                    
  void  ratio_test_1__t_j( Tag_true  is_linear);
  void  ratio_test_1__t_j( Tag_false is_linear);

  void  ratio_test_2( Tag_true  is_linear);
  void  ratio_test_2( Tag_false is_linear);
  void  ratio_test_2__p( Tag_true  has_no_inequalities);
  void  ratio_test_2__p( Tag_false has_no_inequalities);   
  void  ratio_test_2__p( bool has_no_inequalities) {
    if (has_no_inequalities)
      ratio_test_2__p (Tag_true());
    else
      ratio_test_2__p (Tag_false());
  }

  // update
  void  update_1( );
  void  update_1( Tag_true  is_linear);
  void  update_1( Tag_false is_linear);

  void  update_2( Tag_true  is_linear);
  void  update_2( Tag_false is_linear);

  void  replace_variable( );
  void  replace_variable( Tag_true  has_no_inequalities);
  void  replace_variable( Tag_false has_no_inequalities);
  void  replace_variable( bool has_no_inequalities) {
    if (has_no_inequalities)
      replace_variable (Tag_true());
    else
      replace_variable (Tag_false());
  }

  void  replace_variable_original_original( );
  // update of the vector r
  void  replace_variable_original_original_upd_r(Tag_true
						 /*is_nonnegative*/);
  void  replace_variable_original_original_upd_r(Tag_false
						 /*is_nonnegative*/);

  void  replace_variable_original_slack( );
  // update of the vector r
  void  replace_variable_original_slack_upd_r(Tag_true /*is_nonnegative*/);
  void  replace_variable_original_slack_upd_r(Tag_false /*is_nonnegative*/);

  void  replace_variable_slack_original( );
  // update of the vector r
  void  replace_variable_slack_original_upd_r(Tag_true /*is_nonnegative*/);
  void  replace_variable_slack_original_upd_r(Tag_false /*is_nonnegative*/);
    
  void  replace_variable_slack_slack( );
  // update of the vector r
  void  replace_variable_slack_slack_upd_r(Tag_true /*is_nonnegative*/);
  void  replace_variable_slack_slack_upd_r(Tag_false /*is_nonnegative*/);
    
  void  remove_artificial_variable_and_constraint( );
  // update of the vector r
  void  remove_artificial_variable_and_constraint_upd_r(Tag_true
							/*is_nonnegative*/);
  void  remove_artificial_variable_and_constraint_upd_r(Tag_false
							/*is_nonnegative*/);    
    
  void  expel_artificial_variables_from_basis( );
    
  // update that occurs only with upper bounding in ratio test step 1
  void  enter_and_leave_variable( );

  void  enter_variable( );
  // update of the vectors w and r
  void  enter_variable_original_upd_w_r(Tag_true /*is_nonnegative*/);
  void  enter_variable_original_upd_w_r(Tag_false /*is_nonnegative*/);
  void  enter_variable_slack_upd_w_r(Tag_true /*is_nonnegative*/);
  void  enter_variable_slack_upd_w_r(Tag_false /*is_nonnegative*/);
    
  void  leave_variable( );
  // update of the vectors w and r
  void  leave_variable_original_upd_w_r(Tag_true /*is_nonnegative*/);    
  void  leave_variable_original_upd_w_r(Tag_false /*is_nonnegative*/);
  void  leave_variable_slack_upd_w_r(Tag_true /*is_nonnegative*/);
  void  leave_variable_slack_upd_w_r(Tag_false /*is_nonnegative*/);
    
  void  z_replace_variable( );
  void  z_replace_variable( Tag_true has_no_inequalities);
  void  z_replace_variable( Tag_false has_no_inequalities);
  void  z_replace_variable( bool has_no_inequalities) {
    if (has_no_inequalities) 
      z_replace_variable (Tag_true());
    else
      z_replace_variable (Tag_false());
  }
    
  void  z_replace_variable_original_by_original( );
  // update of the vectors w and r
  void  z_replace_variable_original_by_original_upd_w_r(Tag_true 
							/*is_nonnegative*/);
  void  z_replace_variable_original_by_original_upd_w_r(Tag_false 
							/*is_nonnegative*/);
    
  void  z_replace_variable_original_by_slack( );
  // update of the vectors w and r    
  void  z_replace_variable_original_by_slack_upd_w_r(Tag_true 
						     /*is_nonnegative*/);
  void  z_replace_variable_original_by_slack_upd_w_r(Tag_false
						     /*is_nonnegative*/);
    
  void  z_replace_variable_slack_by_original( );
  // update of the vectors w and r
  void  z_replace_variable_slack_by_original_upd_w_r(Tag_true
						     /*is_nonnegative*/);
  void  z_replace_variable_slack_by_original_upd_w_r(Tag_false
						     /*is_nonnegative*/);
    
  void  z_replace_variable_slack_by_slack( );
  // update of the vectors w and r
  void  z_replace_variable_slack_by_slack_upd_w_r(Tag_true
						  /*is_nonnegative*/);
  void  z_replace_variable_slack_by_slack_upd_w_r(Tag_false
						  /*is_nonnegative*/);
    
  // update of the parts r_C and r_S_B
  void  update_r_C_r_S_B__j(ET& x_j);
  void  update_r_C_r_S_B__j_i(ET& x_j, ET& x_i);
  void  update_r_C_r_S_B__i(ET& x_i);
    
  // update of w and r_B_O 
  void  update_w_r_B_O__j(ET& x_j);
  void  update_w_r_B_O__j_i(ET& x_j, ET& x_i);
  void  update_w_r_B_O__i(ET& x_i);
    
    
  bool  basis_matrix_stays_regular( );

  // current solution
  void  compute_solution(Tag_true  /*is_nonnegative*/);
  void  compute_solution(Tag_false /*is_nonnegative*/);

  void  compute__x_B_S( Tag_false  has_no_inequalities,
			Tag_false /*is_nonnegative*/);
  void  compute__x_B_S( Tag_false  has_no_inequalities,
			Tag_true  /*is_nonnegative*/);
  void  compute__x_B_S( Tag_true  has_no_inequalities,
			Tag_false /*is_nonnegative*/);
  void  compute__x_B_S( Tag_true  has_no_inequalities,
			Tag_true  /*is_nonnegative*/);
  void  compute__x_B_S( bool  has_no_inequalities,
			Tag_true  is_nonnegative) {
    if (has_no_inequalities)
      compute__x_B_S (Tag_true(), is_nonnegative);
    else
      compute__x_B_S (Tag_false(), is_nonnegative);
  }
    
  void  compute__x_B_S( bool  has_no_inequalities,
			Tag_false  is_nonnegative) {
    if (has_no_inequalities)
      compute__x_B_S (Tag_true(), is_nonnegative);
    else
      compute__x_B_S (Tag_false(), is_nonnegative);
  }  

  void  multiply__A_S_BxB_O( Value_iterator in, Value_iterator out) const;
    
  ET    multiply__A_ixO(int row) const;
  void  multiply__A_CxN_O(Value_iterator out) const;
  bool  check_r_C(Tag_true  /*is_nonnegative*/) const;
  bool  check_r_C(Tag_false /*is_nonnegative*/) const;
    
  void  multiply__A_S_BxN_O(Value_iterator out) const;
  bool  check_r_S_B(Tag_true  /*is_nonnegative*/) const;
  bool  check_r_S_B(Tag_false /*is_nonnegative*/) const;
    
  void  multiply__2D_B_OxN_O(Value_iterator out) const;
  bool  check_r_B_O(Tag_true  /*is_nonnegative*/) const;
  bool  check_r_B_O(Tag_false /*is_nonnegative*/) const;
        
  void  multiply__2D_OxN_O(Value_iterator out) const;
  bool  check_w(Tag_true  /*is_nonnegative*/) const;
  bool  check_w(Tag_false /*is_nonnegative*/) const;
    
  // utility routines for QP's in nonstandard form:
  ET original_variable_value_under_bounds(int i) const;
  ET nonbasic_original_variable_value (int i) const;

public: 
  // for original variables
  ET variable_numerator_value(int i) const;
  ET unbounded_direction_value(int i) const;
  ET lambda_numerator(int i) const
  {
    // we use the vector lambda which conforms to C (basic constraints)
    CGAL_qpe_assertion (i >= 0);
    CGAL_qpe_assertion (i <= qp_m);
    if (no_ineq)
      return lambda[i];
    else {
      int k = in_C[i];     // position of i in C
      if (k != -1) 
	return lambda[k];
      else 
	return et0;
    }   
}

private:
  // check basis inverse
  bool  check_basis_inverse( );
  bool  check_basis_inverse( Tag_true  is_linear);
  bool  check_basis_inverse( Tag_false is_linear);

  // diagnostic output
  void  print_program ( ) const;
  void  print_basis   ( ) const;
  void  print_solution( ) const;
  void  print_ratio_1_original(int k, const ET& x_k, const ET& q_k);
  void  print_ratio_1_slack(int k, const ET& x_k, const ET& q_k);

  const char*  variable_type( int k) const;
    
  // ensure container size
  template <class Container>
  void ensure_size(Container& c, typename Container::size_type desired_size) {
    typedef typename Container::value_type Value_type;
    for (typename Container::size_type i=c.size(); i < desired_size; ++i) {
      c.push_back(Value_type());
    }
  }
    
private:

private:  // (inefficient) access to bounds of variables:
  // Given an index of an original or slack variable, returns whether
  // or not the variable has a finite lower bound.
  bool has_finite_lower_bound(int i) const;

  // Given an index of an original or slack variable, returns whether
  // or not the variable has a finite upper bound.
  bool has_finite_upper_bound(int i) const;

  // Given an index of an original or slack variable, returns its
  // lower bound.
  ET lower_bound(int i) const;

  // Given an index of an original variable, returns its upper bound.
  ET upper_bound(int i) const;

  struct Bnd { // (inefficient) utility class representing a possibly
    // infinite bound
    enum Kind { MINUS_INF=-1, FINITE=0, PLUS_INF=1 };
    const Kind kind;      // whether the bound is finite or not
    const ET value;       // bound's value in case it is finite

    Bnd(bool is_upper, bool is_finite, const ET& value) 
      : kind(is_upper? (is_finite? FINITE : PLUS_INF) :
	     (is_finite? FINITE : MINUS_INF)),
	value(value) {}
    Bnd(Kind kind, const ET& value) : kind(kind), value(value) {}
    
    bool operator==(const ET& v) const { return kind == FINITE && value == v; }
    bool operator==(const Bnd& b) const {
      return kind == b.kind && (kind != FINITE || value == b.value);
    }
    bool operator!=(const Bnd& b) const { return !(*this == b); }
    bool operator<(const ET& v) const { return kind == FINITE && value < v; }
    bool operator<(const Bnd& b) const {
      return kind < b.kind ||
	(kind == b.kind && kind == FINITE && value < b.value);
    }
    bool operator<=(const Bnd& b) const { return *this < b || *this == b; }
    bool operator>(const ET& v) const { return kind == FINITE && value > v; }
    bool operator>(const Bnd& b)  const { return !(*this <= b); }
    bool operator>=(const Bnd& b) const { return !(*this < b); }
    
    Bnd operator*(const ET& f) const { return Bnd(kind, value*f); }
  };

  // Given an index of an original, slack, or artificial variable,
  // return its lower bound.
  Bnd lower_bnd(int i) const;

  // Given an index of an original, slack, or artificial variable,
  // return its upper bound.
  Bnd upper_bnd(int i) const;

private:
  bool is_value_correct() const;
  // ----------------------------------------------------------------------------

  // ===============================
  // class implementation (template)
  // ===============================

public:

  // pricing
  // -------
  // The solver provides three methods to compute mu_j; the first
  // two below take additional information (which the pricing
  // strategy either provides in exact- or NT-form), and the third
  // simply does the exact computation. (Note: internally, we use
  // the third version, too, see ratio_test_1__t_j().)

  // computation of mu_j with standard form
  template < class RndAccIt1, class RndAccIt2, class NT >  
  NT
  mu_j( int j, RndAccIt1 lambda_it, RndAccIt2 x_it, const NT& dd) const
  {
    NT  mu_j;

    if ( j < qp_n) {                                // original variable

      // [c_j +] A_Cj^T * lambda_C
      mu_j = ( is_phaseI ? NT( 0) : dd * NT(*(qp_c+ j)));
      mu_j__linear_part( mu_j, j, lambda_it, no_ineq);

      // ... + 2 D_Bj^T * x_B
      mu_j__quadratic_part( mu_j, j, x_it, Is_linear());

    } else {                                        // slack or artificial

      mu_j__slack_or_artificial( mu_j, j, lambda_it, dd,
				 no_ineq);

    }

    return mu_j;
  }
    
  // computation of mu_j with upper bounding
  template < class RndAccIt1, class RndAccIt2, class NT >  
  NT
  mu_j( int j, RndAccIt1 lambda_it, RndAccIt2 x_it, const NT& w_j,
	const NT& dd) const
  {
    NT  mu_j;

    if ( j < qp_n) {                                // original variable

      // [c_j +] A_Cj^T * lambda_C
      mu_j = ( is_phaseI ? NT( 0) : dd * NT(*(qp_c+ j)));
      mu_j__linear_part( mu_j, j, lambda_it, no_ineq);

      // ... + 2 D_Bj^T * x_B + 2 D_Nj x_N
      mu_j__quadratic_part( mu_j, j, x_it, w_j, dd, Is_linear());

    } else {                                        // slack or artificial

      mu_j__slack_or_artificial( mu_j, j, lambda_it, dd,
				 no_ineq);

    }

    return mu_j;
  }

  // computation of mu_j (exact, both for upper bounding and standard form)
  ET
  mu_j( int j) const
  {
    CGAL_qpe_assertion(!is_basic(j));
    
    if (!check_tag(Is_nonnegative()) &&
	!check_tag(Is_linear()) &&
	!is_phaseI && is_original(j)) {
      return mu_j(j,
		  lambda.begin(),
		  basic_original_variables_numerator_begin(),
		  w_j_numerator(j),
		  variables_common_denominator());
    } else {
      return mu_j(j,
		  lambda.begin(),
		  basic_original_variables_numerator_begin(),
		  variables_common_denominator());
    }
  }

private:

  // pricing (private helper functions)
  // ----------------------------------
  template < class NT, class It > inline                      // no ineq.
  void
  mu_j__linear_part( NT& mu_j, int j, It lambda_it, Tag_true) const
  {
    mu_j += inv_M_B.inner_product_l( lambda_it, *(qp_A+ j));
  }

  template < class NT, class It > inline                      // has ineq.
  void
  mu_j__linear_part( NT& mu_j, int j, It lambda_it, Tag_false) const
  {
    mu_j += inv_M_B.inner_product_l
      ( lambda_it,
	A_by_index_iterator( C.begin(),
			     A_by_index_accessor( *(qp_A + j))));
  }

  template < class NT, class It > inline                     
  void
  mu_j__linear_part( NT& mu_j, int j, It lambda_it, 
		     bool has_no_inequalities) const {
    if (has_no_inequalities) 
      mu_j__linear_part (mu_j, j, lambda_it, Tag_true());
    else
      mu_j__linear_part (mu_j, j, lambda_it, Tag_false());     
  }

  template < class NT, class It > inline          // LP case, standard form
  void
  mu_j__quadratic_part( NT&, int, It, Tag_true) const
  {
    // nop
  }
    
  template < class NT, class It > inline          // LP case, upper bounded
  void
  mu_j__quadratic_part( NT&, int, It, const NT& /*w_j*/, const NT& /*dd*/,
			Tag_true) const
  {
    // nop
  }    

  template < class NT, class It > inline          // QP case, standard form
  void
  mu_j__quadratic_part( NT& mu_j, int j, It x_it, Tag_false) const
  {
    if ( is_phaseII) {
      // 2 D_Bj^T * x_B
      mu_j += inv_M_B.inner_product_x
	( x_it,
	  D_pairwise_iterator_input_type( B_O.begin(),
			       D_pairwise_accessor_input_type(qp_D, j)));
    }
  }

  template < class NT, class It > inline          // QP case, upper bounded
  void
  mu_j__quadratic_part( NT& mu_j, int j, It x_it, const NT& w_j,
			const NT& dd, Tag_false) const
  {
    if ( is_phaseII) {
      mu_j += dd * w_j;
      // 2 D_Bj^T * x_B
      mu_j += inv_M_B.inner_product_x
	( x_it,
	  D_pairwise_iterator_input_type( B_O.begin(),
			       D_pairwise_accessor_input_type(qp_D, j)));
    }
  }


  template < class NT, class It >  inline                     // no ineq.
  void
  mu_j__slack_or_artificial( NT& mu_j, int j, It lambda_it, const NT& dd, Tag_true) const
  {
    j -= qp_n;
    // artificial variable
    // A_j^T * lambda
    mu_j = lambda_it[ j];
    if ( art_A[ j].second) mu_j = -mu_j;

    // c_j + ...
    mu_j += dd*NT(aux_c[ j]);

  }

  template < class NT, class It >  inline                     // has ineq.
  void
  mu_j__slack_or_artificial( NT& mu_j, int j, It lambda_it, const NT& dd, Tag_false) const
  {
    j -= qp_n;

    if ( j < static_cast<int>(slack_A.size())) {                 // slack variable

      // A_Cj^T * lambda_C
      mu_j = lambda_it[ in_C[ slack_A[ j].first]];
      if ( slack_A[ j].second) mu_j = -mu_j;

    } else {                                        // artificial variable
      j -= static_cast<int>(slack_A.size());

      // A_Cj^T * lambda_C
      mu_j = lambda_it[ in_C[ art_A[ j].first]];
      if ( art_A[ j].second) mu_j = -mu_j;

      // c_j + ...
      mu_j += dd*NT(aux_c[ j]);
    }
  }

  template < class NT, class It >  inline
  void
  mu_j__slack_or_artificial( NT& mu_j, int j, It lambda_it, 
			     const NT& dd, bool has_no_inequalities) const {
    if (has_no_inequalities)
      mu_j__slack_or_artificial (mu_j, j, lambda_it, dd, Tag_true());
    else
      mu_j__slack_or_artificial (mu_j, j, lambda_it, dd, Tag_false());
  }

  
  
};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// initialization
// --------------

// transition
// ----------
template < class Q, typename ET, typename Tags >  inline                                 // QP case
void  QP_solver<Q, ET, Tags>::
transition( Tag_false)
{
  typedef  Creator_2< D_iterator, int, 
    D_pairwise_accessor >  D_transition_creator_accessor;

  typedef  Creator_2< Index_iterator, D_pairwise_accessor,
    D_pairwise_iterator >  D_transition_creator_iterator;

  // initialization of vector w and vector r_B_O:
  if (!check_tag(Is_nonnegative())) {
    init_w();                      
    init_r_B_O();
  }

  // here is what we need in the transition: an iterator that steps through 
  // the basic indices, where dereferencing
  // yields an iterator through the corresponding row of D, restricted 
  // to the basic indices. This means that we select the principal minor of D 
  // corresponding to the current basis.
 
  // To realize this, we transform B_O.begin() via the function h where
  //   h(i) = D_pairwise_iterator
  //           (B_O.begin(), 
  //            D_pairwise_accessor(qp_D, i))


  inv_M_B.transition 
    (boost::make_transform_iterator 
     (B_O.begin(),
      boost::bind 
      (D_transition_creator_iterator(), B_O.begin(), 
       boost::bind (D_transition_creator_accessor(), qp_D, _1))));
}

template < typename Q, typename ET, typename Tags >  inline                                 // LP case
void  QP_solver<Q, ET, Tags>::
transition( Tag_true)
{
  inv_M_B.transition();
}

// ratio test
// ----------
template < typename Q, typename ET, typename Tags > inline                                  // LP case
void  QP_solver<Q, ET, Tags>::
ratio_test_init__2_D_Bj( Value_iterator, int, Tag_true)
{
  // nop
}

template < typename Q, typename ET, typename Tags > inline                                  // QP case
void  QP_solver<Q, ET, Tags>::
ratio_test_init__2_D_Bj( Value_iterator two_D_Bj_it, int j_, Tag_false)
{
  if ( is_phaseII) {
    ratio_test_init__2_D_Bj( two_D_Bj_it, j_,
			     Tag_false(), no_ineq);
  }
}

template < typename Q, typename ET, typename Tags > inline                                  // QP, no ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_init__2_D_Bj( Value_iterator two_D_Bj_it, int j_, Tag_false,
			 Tag_true )
{
  // store exact version of `2 D_{B_O,j}'
  D_pairwise_accessor  d_accessor( qp_D, j_);
  std::copy( D_pairwise_iterator( B_O.begin(), d_accessor),
	     D_pairwise_iterator( B_O.end  (), d_accessor),
	     two_D_Bj_it);
}

template < typename Q, typename ET, typename Tags > inline                                  // QP, has ineq
void  QP_solver<Q, ET, Tags>::
ratio_test_init__2_D_Bj( Value_iterator two_D_Bj_it, int j_, Tag_false,
			 Tag_false)
{
  // store exact version of `2 D_{B_O,j}'
  if ( j_ < qp_n) {                               // original variable
    ratio_test_init__2_D_Bj( two_D_Bj_it, j_, Tag_false(), Tag_true());
  } else {                                        // slack variable
    std::fill_n( two_D_Bj_it, B_O.size(), et0);
  }
}

template < typename Q, typename ET, typename Tags >  inline                                 // LP case
void  QP_solver<Q, ET, Tags>::
ratio_test_1__q_x_O( Tag_true)
{
  inv_M_B.multiply_x( A_Cj.begin(), q_x_O.begin());
}

template < typename Q, typename ET, typename Tags >  inline                                 // QP case
void  QP_solver<Q, ET, Tags>::
ratio_test_1__q_x_O( Tag_false)
{
  if ( is_phaseI) {                                   // phase I
    inv_M_B.multiply_x(     A_Cj.begin(),    q_x_O.begin());
  } else {                                            // phase II
    inv_M_B.multiply  (     A_Cj.begin(), two_D_Bj.begin(),
			    q_lambda.begin(),    q_x_O.begin());
  }
}

template < typename Q, typename ET, typename Tags >  inline                                 // no ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_1__q_x_S( Tag_true)
{
  // nop
}

template < typename Q, typename ET, typename Tags >  inline                                 // has ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_1__q_x_S( Tag_false)
{
  // A_S_BxB_O * q_x_O
  multiply__A_S_BxB_O( q_x_O.begin(), q_x_S.begin());

  // ( A_S_BxB_O * q_x_O) - A_S_Bxj
  if ( j < qp_n) {
    std::transform( q_x_S.begin(),
		    q_x_S.begin()+S_B.size(),
		    A_by_index_iterator( S_B.begin(),
					 A_by_index_accessor( *(qp_A + j))),
		    q_x_S.begin(),
		    compose2_2( std::minus<ET>(),
				Identity<ET>(),
				std::bind1st( std::multiplies<ET>(), d)));
  }

  // q_x_S = -+ ( A_S_BxB_O * q_x_O - A_S_Bxj)
  Value_iterator  q_it = q_x_S.begin();
  Index_iterator  i_it;
  for ( i_it = B_S.begin(); i_it != B_S.end(); ++i_it, ++q_it) {
    if ( ! slack_A[ *i_it - qp_n].second) *q_it = -(*q_it);
  }
}

template < typename Q, typename ET, typename Tags >  inline                                 // no check
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_i( Index_iterator, Index_iterator,
		   Value_iterator, Value_iterator, Tag_true)
{
  // nop
}

template < typename Q, typename ET, typename Tags >  inline                                 // check
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_i( Index_iterator i_it, Index_iterator end_it,
		   Value_iterator x_it, Value_iterator   q_it, Tag_false)
{
  // check `t_i's
  for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it) {
    // BLAND rule: In case the ratios are the same, only update if the new index
    // is smaller. The special artificial variable is always made to leave first.
    if ( (*q_it > et0) && (
                           (( *x_it * q_i) < ( x_i * *q_it)) ||
                           ( (*i_it < i) && (i != art_s_i) && (( *x_it * q_i) == ( x_i * *q_it)) )
                           )
        ) {
      i = *i_it; x_i = *x_it; q_i = *q_it;
    }
  }
}

template < typename Q, typename ET, typename Tags >  inline                                 // LP case
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_j( Tag_true)
{
  // nop
}

template < typename Q, typename ET, typename Tags >  inline                                 // QP case
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_j( Tag_false)
{
  if ( is_phaseII) {

    // compute `nu' and `mu_j' 
    mu = mu_j(j);
    nu = inv_M_B.inner_product(     A_Cj.begin(), two_D_Bj.begin(),
				    q_lambda.begin(),    q_x_O.begin());
    if ( j < qp_n) {                                // original variable
      nu -= d*ET( (*(qp_D + j))[ j]);
    }
    CGAL_qpe_assertion_msg(nu <= et0,
			   "nu <= et0 violated -- is your D matrix positive semidefinite?");

    // check `t_j'
    CGAL_qpe_assertion(mu != et0);
    // bg: formula below compares abs values, assuming mu < 0
    if ( ( nu < et0) && ( ( (mu < et0 ? mu : -mu) * q_i) > ( x_i * nu))) {
      i = -1; q_i = et1;
    }
  }
}

template < typename Q, typename ET, typename Tags >  inline                                 // LP case
void  QP_solver<Q, ET, Tags>::
ratio_test_2( Tag_true)
{
  // nop
}

template < typename Q, typename ET, typename Tags >  inline                                 // no ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_2__p( Tag_true)
{
  // get column index of entering variable in basis
  int  col = in_B[ j];
 
  CGAL_qpe_assertion( col >= 0);
  col += l;

  // get (last) column of `M_B^{-1}' (Note: `p_...' is stored in `q_...')
  Value_iterator  it;
  int             row;
  unsigned int    k;
  for (   k = 0,            row = 0,   it = q_lambda.begin();
	  k < C.size();
	  ++k,              ++row,     ++it                   ) {
    *it = inv_M_B.entry( row, col);
  }
  for (   k = 0,            row = l,   it = q_x_O.begin();
	  k < B_O.size();
	  ++k,              ++row,     ++it                   ) {
    *it = inv_M_B.entry( row, col);
  }
}

template < typename Q, typename ET, typename Tags >  inline                                 // has ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_2__p( Tag_false)
{
  Value_iterator  v_it;
  Index_iterator  i_it;

  // compute 'p_lambda' and 'p_x_O' (Note: `p_...' is stored in `q_...')
  // -------------------------------------------------------------------

  // type of entering variable
  if ( j < qp_n) {                                        // original

    // use 'no_ineq' variant
    ratio_test_2__p( Tag_true());

  } else {                                                // slack

    j -= qp_n;

    // get column A_{S_j,B_O}^T (i.e. row of A_{S_B,B_O})
    int             row  = slack_A[ j].first;
    bool            sign = slack_A[ j].second;

    for (   i_it =  B_O.begin(),   v_it = tmp_x.begin();
	    i_it != B_O.end();
	    ++i_it,                ++v_it                ) {
      *v_it = ( sign ? 
		*((*(qp_A+ *i_it))+ row) : - (*((*(qp_A + *i_it))+ row)));
    }

    // compute  ( p_l | p_x_O )^T = M_B^{-1} * ( 0 | A_{S_j,B_O} )^T
    std::fill_n( tmp_l.begin(), C.size(), et0);
    inv_M_B.multiply( tmp_l     .begin(), tmp_x  .begin(),
		      q_lambda.begin(),   q_x_O.begin());

    j += qp_n;
  }

  // compute 'p_x_S'
  // ---------------
  // A_S_BxB_O * p_x_O
  multiply__A_S_BxB_O( q_x_O.begin(), q_x_S.begin());

  // p_x_S = +- ( A_S_BxB_O * p_x_O)
  for (   i_it =  B_S.begin(),   v_it = q_x_S.begin();
	  i_it != B_S.end();
	  ++i_it,                ++v_it                ) {
    if ( ! slack_A[ *i_it - qp_n].second) *v_it = -(*v_it);
  }
}

// update
// ------
template < typename Q, typename ET, typename Tags >  inline                                 // LP case
void  QP_solver<Q, ET, Tags>::
update_1( Tag_true)
{
  // replace leaving with entering variable
  if ((i == j) && (i >= 0)) {
    enter_and_leave_variable();
  } else {
    replace_variable();
  }
}

template < typename Q, typename ET, typename Tags >  inline                                 // QP case
void  QP_solver<Q, ET, Tags>::
update_1( Tag_false)
{
  if ( is_phaseI) {                                   // phase I

    // replace leaving with entering variable
    if ((i == j) && (i >= 0)) {
      enter_and_leave_variable();
    } else {
      replace_variable();
    }

  } else {                                            // phase II
        
    if ((i == j) && (i >= 0)) {
      enter_and_leave_variable();
    } else {

      if ( ( i >= 0) && basis_matrix_stays_regular()) {

	// leave variable from basis, if
	// - some leaving variable was found  and
	// - basis matrix stays regular
	leave_variable();

      } else {

	// enter variable into basis, if
	// - no leaving variable was found  or
	// - basis matrix would become singular when variable i leaves

	if ( i < 0 ) {
	  enter_variable();
	} else {
	  z_replace_variable();
	}
      }
    }
  }
}

template < typename Q, typename ET, typename Tags >  inline                                 // LP case
void  QP_solver<Q, ET, Tags>::
update_2( Tag_true)
{
  // nop
}

template < typename Q, typename ET, typename Tags >  inline                                 // no ineq.
void  QP_solver<Q, ET, Tags>::
replace_variable( Tag_true)
{
  replace_variable_original_original();
  strategyP->leaving_basis( i);
}

template < typename Q, typename ET, typename Tags >  inline                                 // has ineq.
void  QP_solver<Q, ET, Tags>::
replace_variable( Tag_false)
{
  // determine type of variables
  bool  enter_original = ( (j < qp_n) || (j >= static_cast<int>( qp_n+slack_A.size())));
  bool  leave_original = ( (i < qp_n) || (i >= static_cast<int>( qp_n+slack_A.size())));

  // update basis & basis inverse
  if ( leave_original) {
    if ( enter_original) {                              // orig  <--> orig
      replace_variable_original_original();
    } else {                                            // slack <--> orig
      replace_variable_slack_original();
    }

    // special artificial variable removed?
    if ( is_phaseI && ( i == art_s_i)) {
      // remove the fake column - it corresponds
      // to the special artificial variable which is
      // (like all artificial variables) not needed
      // anymore once it leaves the basis. Note:
      // regular artificial variables are only removed
      // from the problem after phase I
      // art_s_i == -1 -> there is no special artificial variable
      // art_s_i == -2 -> there was a special artificial variable, 
      // but has been removed  
      art_s_i = -2;
      art_A.pop_back();
      CGAL_qpe_assertion(in_B[in_B.size()-1] == -1); // really removed?
      in_B.pop_back();
      // BG: shouldn't the pricing strategy be notfied also here?
    } else {
      strategyP->leaving_basis( i);
    }
  } else {
    if ( enter_original) {                              // orig  <--> slack
      replace_variable_original_slack();
    } else {                                            // slack <--> slack
      replace_variable_slack_slack();
    }
    strategyP->leaving_basis( i);
  }
}

template < typename Q, typename ET, typename Tags >  inline
bool  QP_solver<Q, ET, Tags>::
basis_matrix_stays_regular()
{
  CGAL_qpe_assertion( is_phaseII);
  int new_row, k;
    
  if ( has_ineq && (i >= qp_n)) {	// slack variable
    new_row = slack_A[ i-qp_n].first;
    A_row_by_index_accessor  a_accessor =
      boost::bind (A_accessor( qp_A, 0, qp_n), _1, new_row);
    std::copy( A_row_by_index_iterator( B_O.begin(), a_accessor),
	       A_row_by_index_iterator( B_O.end  (), a_accessor),
	       tmp_x.begin());	   
    inv_M_B.multiply( tmp_x.begin(),                        // dummy (not used)
		      tmp_x.begin(), tmp_l_2.begin(), tmp_x_2.begin(),
		      Tag_false(),                                 // QP
		      Tag_false());                             // ignore 1st argument
    return ( -inv_M_B.inner_product_x( tmp_x_2.begin(), tmp_x.begin()) != et0);

	
  } else {						// check original variable
    k = l+in_B[ i];
    return ( inv_M_B.entry( k, k) != et0);
  }

  /* ToDo: check, if really not needed in 'update_1':
     - basis has already minimal size  or
     || ( B_O.size()==C.size()) 
  */
}

// current solution
// ----------------
template < typename Q, typename ET, typename Tags >  inline             // no inequalities, upper bounded
void  QP_solver<Q, ET, Tags>::
compute__x_B_S( Tag_true  /*has_equalities_only_and_full_rank*/,
                Tag_false /*is_nonnegative*/)
{
  // nop
}

template < typename Q, typename ET, typename Tags >  inline             // no inequalities, standard form
void  QP_solver<Q, ET, Tags>::
compute__x_B_S( Tag_true /*has_equalities_only_and_full_rank*/,
                Tag_true /*is_nonnegative*/)
{
  // nop
}


template < typename Q, typename ET, typename Tags >  inline             // has inequalities, upper bounded
void  QP_solver<Q, ET, Tags>::
compute__x_B_S( Tag_false /*has_equalities_only_and_full_rank*/,
                Tag_false /*is_nonnegative*/)
{
  // A_S_BxB_O * x_B_O
  multiply__A_S_BxB_O( x_B_O.begin(), x_B_S.begin());

  // b_S_B - ( A_S_BxB_O * x_B_O)
  B_by_index_accessor  b_accessor( qp_b);
  std::transform( B_by_index_iterator( S_B.begin(), b_accessor),
		  B_by_index_iterator( S_B.end  (), b_accessor),
		  x_B_S.begin(),
		  x_B_S.begin(),
		  compose2_2( std::minus<ET>(),
			      std::bind1st( std::multiplies<ET>(), d),
			      Identity<ET>()));
				
  // b_S_B - ( A_S_BxB_O * x_B_O) - r_S_B
  std::transform(x_B_S.begin(), x_B_S.begin()+S_B.size(),
		 r_S_B.begin(), x_B_S.begin(),
		 compose2_2(std::minus<ET>(),
			    Identity<ET>(),
			    std::bind1st( std::multiplies<ET>(), d)));
                        

  // x_B_S = +- ( b_S_B - A_S_BxB_O * x_B_O)
  Value_iterator  x_it = x_B_S.begin();
  Index_iterator  i_it;
  for ( i_it = B_S.begin(); i_it != B_S.end(); ++i_it, ++x_it) {
    if ( slack_A[ *i_it - qp_n].second) *x_it = -(*x_it);
  }
       
}



template < typename Q, typename ET, typename Tags >  inline             // has inequalities, standard form
void  QP_solver<Q, ET, Tags>::
compute__x_B_S( Tag_false /*has_equalities_only_and_full_rank*/,
                Tag_true  /*is_nonnegative*/)
{
  // A_S_BxB_O * x_B_O
  multiply__A_S_BxB_O( x_B_O.begin(), x_B_S.begin());

  // b_S_B - ( A_S_BxB_O * x_B_O)
  B_by_index_accessor  b_accessor( qp_b);
  std::transform( B_by_index_iterator( S_B.begin(), b_accessor),
		  B_by_index_iterator( S_B.end  (), b_accessor),
		  x_B_S.begin(),
		  x_B_S.begin(),
		  compose2_2( std::minus<ET>(),
			      std::bind1st( std::multiplies<ET>(), d),
			      Identity<ET>()));

  // x_B_S = +- ( b_S_B - A_S_BxB_O * x_B_O)
  Value_iterator  x_it = x_B_S.begin();
  Index_iterator  i_it;
  for ( i_it = B_S.begin(); i_it != B_S.end(); ++i_it, ++x_it) {
    if ( slack_A[ *i_it - qp_n].second) *x_it = -(*x_it);
  }
       
}

} //namespace CGAL

#include <CGAL/QP_solver/Unbounded_direction.h>
#include <CGAL/QP_solver/QP_solver_nonstandardform_impl.h>
#include <CGAL/QP_solver/QP_solver_bounds_impl.h>
#include <CGAL/QP_solver/QP_solver_impl.h>

#endif // CGAL_QP_SOLVER_H

// ===== EOF ==================================================================
