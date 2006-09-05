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
    
  // get properties of quadratic program:
  int  w = this->solver().number_of_working_variables();

  // loop over all non-basic variables:
  int  j,  min_j  = -1;
  ET min_mu = this->et0;
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

      const ET mu = this->mu_j(j);
      // from pricing strategy base class
      price_dantzig (j, mu, this->et0, min_j, min_mu, direction); 
      if (min_j >= 0) return j;
    }
  }
  this->vout() << std::endl;

  // fallback option
  return -1;
    
}


CGAL_END_NAMESPACE

#endif // CGAL_QP_EXACT_BLAND_PRICING_H

// ===== EOF ==================================================================
