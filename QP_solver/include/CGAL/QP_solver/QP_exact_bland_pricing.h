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
// $URL: $
// $Id:  $
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp <fransw@inf.ethz.ch>
//                 Kaspar Fischer <fischerk@inf.ethz.ch>

#ifndef CGAL_QP_EXACT_BLAND_PRICING_H
#define CGAL_QP_EXACT_BLAND_PRICING_H

// includes
#include <CGAL/QP_solver/QP_pricing_strategy.h>

CGAL_BEGIN_NAMESPACE

// =================
// class declaration
// =================
template < class Rep_ >
class QP_exact_bland_pricing;

// ===============
// class interface
// ===============
template < class Rep_ >
class QP_exact_bland_pricing : public QP_pricing_strategy<Rep_> {

  // self
  typedef  Rep_                         Rep;
  typedef  QP_pricing_strategy<Rep>    Base;
  typedef  QP_exact_bland_pricing<Rep>  Self;

  // types from the base class
  typedef  typename Base::ET                      ET;
  typedef  typename Base::Is_in_standard_form     Is_in_standard_form;
  typedef  typename CGAL::QP_solver<Rep>          QP_solver;

 public:

  // creation
  QP_exact_bland_pricing();

  // operations
  int  pricing(int& direction );
    
  // cleanup
  ~QP_exact_bland_pricing() { };
    
 private:
  int pricing_helper(int& direction, Tag_true  is_in_standard_form);
  int pricing_helper(int& direction, Tag_false is_in_standard_form);
};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// construction
template < class Rep_ >  inline
QP_exact_bland_pricing<Rep_>::
QP_exact_bland_pricing()
  : QP_pricing_strategy<Rep_>("full exact")
{ }
    
// operations
template < class Rep_ >
int  QP_exact_bland_pricing<Rep_>::
pricing(int& direction )
{
  return (pricing_helper(direction, Is_in_standard_form()));
}

template < class Rep_ >
int  QP_exact_bland_pricing<Rep_>::
pricing_helper(int& direction, Tag_true is_in_standard_form)
{
  // get properties of quadratic program:
  int  w = this->solver().number_of_working_variables();

  // loop over all non-basic variables:
  int  j;
  ET   mu;
  for (j = 0; j < w; ++j) {

    // variable non-basic?
    if (!this->solver().is_basic(j)) {
	
      // don't price artificial variables:
      if (this->solver().is_artificial(j)) {
	CGAL_qpe_debug { 
	  this->vout() << "mu_" << j << ": artificial [ not priced ]"
		       << std::endl;
	}
	continue;
      }

      // compute mu_j:
      mu = this->mu_j(j);

      CGAL_qpe_debug {
	this->vout() << "mu_" << j << ": " << mu << std::endl;
      }

      // negative? 
      if (mu < this->et0) return j;
    }
  }
  this->vout() << std::endl;

  // fallback option
  return -1;
    
}

template < class Rep_ >
int  QP_exact_bland_pricing<Rep_>::
pricing_helper(int& direction, Tag_false is_in_standard_form)
{
  typedef typename QP_solver::Bound_index Bound_index;
    
  // get properties of quadratic program:
  int  w = this->solver().number_of_working_variables();

  // loop over all non-basic variables:
  int  j;
  // 
  for (j = 0; j < w; ++j) {

    // variable non-basic?
    if (!this->solver().is_basic(j)) {
	
      // don't price artificial variables:
      if (this->solver().is_artificial(j)) {
	CGAL_qpe_debug { 
	  this->vout() << "mu_" << j << ": artificial [ not priced ]"
		       << std::endl;
	}
	continue;
      }
            
      // original variable:
      if (this->solver().is_original(j)) {
	const Bound_index bnd_ind =
	  this->solver().nonbasic_original_variable_bound_index(j);
	switch (bnd_ind) {
	case QP_solver::LOWER:
	  {
	    // compute mu_j
	    const ET mu = this->mu_j(j);
	    
	    CGAL_qpe_debug { 
	      this->vout() << "mu_" << j << ": " << mu
			   << " LOWER" << std::endl;
	    }
	    
	    if (mu < this->et0) {
	      direction = 1;
	      return j;
	    }
	    break;
	  }
	case QP_solver::ZERO:
	  {
	    // compute mu_j
	    const ET mu = this->mu_j(j);

	    // determine whether the variable is on lower or upper bound, or
	    // somewhere it the middle:
	    //
	    // Note: it cannot be both on the lower and upper bound (as it is
	    // not FIXED).
	    const int where =
	      this->solver().state_of_zero_nonbasic_variable(j);
	    
	    CGAL_qpe_debug { 
	      this->vout() << "mu_" << j << ": " << mu
			   << " ZERO " 
			   << (where == -1? "(LOWER)" : 
			       (where == 0? "(MIDDLE)" : "(UPPER)"))
			   << std::endl;
	    }
	    
	    if (where >= 0 &&       // middle or on upper bound?
		mu > this->et0) {
	      direction = -1;
	      return j;
	    }
	    if (where <= 0 &&       // middle or on lower bound?
		mu < this->et0) {
	      direction = 1;
	      return j;                            
	    }
	    break;
	  }
	case QP_solver::UPPER:
	  {
	    // compute mu_j
	    const ET mu = this->mu_j(j);
	    
	    CGAL_qpe_debug { 
	      this->vout() << "mu_" << j << ": " << mu
			   << " UPPER" << std::endl;
	    }
	    
	    if (mu > this->et0) {
	      direction = -1;
	      return j;
	    }                    
	    break;
	  }
	case QP_solver::FIXED:
	  CGAL_qpe_debug {
	    this->vout() << "Fixed variable " << j << std::endl;
	  }
	  break;
	case QP_solver::BASIC:
	  CGAL_qpe_assertion(false);
	  break;
	}  
      } else {                                    // slack variable
	// compute mu_j
	const ET mu = this->mu_j(j);

	CGAL_qpe_debug {
	  this->vout() << "mu_" << j << ": " << mu 
		       << " LOWER (slack)" << std::endl;
	}

	if (mu < this->et0) {
	  direction = 1;
	  return j;
	}

      }
    }
  }
  this->vout() << std::endl;

  // fallback option
  return -1;
    
}


CGAL_END_NAMESPACE

#endif // CGAL_QP_EXACT_BLAND_PRICING_H

// ===== EOF ==================================================================
