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
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp
//                 Kaspar Fischer

#ifndef CGAL_QP_FULL_EXACT_PRICING_H
#define CGAL_QP_FULL_EXACT_PRICING_H

#include <CGAL/license/QP_solver.h>


// includes
#include <CGAL/QP_solver/QP_pricing_strategy.h>

namespace CGAL {

// =================
// class declaration
// =================
template < typename Q, typename ET, typename Tags >
class QP_full_exact_pricing;

// ===============
// class interface
// ===============
template < typename Q, typename ET, typename Tags >
class QP_full_exact_pricing : public QP_pricing_strategy<Q,ET,Tags> {

  // self
  typedef  QP_pricing_strategy<Q,ET,Tags>    Base;
  typedef  QP_full_exact_pricing<Q,ET,Tags>  Self;

  // types from the base class
  typedef  typename Tags::Is_nonnegative     Is_nonnegative;
  typedef  typename CGAL::QP_solver<Q,ET,Tags>    QP_solver;

  using Base::price_dantzig;

 public:

  // creation
  QP_full_exact_pricing();

  // operations
  int  pricing(int& direction );
    
  // cleanup
  ~QP_full_exact_pricing() { };
    
 private:
  int pricing_helper(int& direction, Tag_true  /*is_in_standard_form*/);
  int pricing_helper(int& direction, Tag_false /*is_in_standard_form*/);
};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// construction
template < typename Q, typename ET, typename Tags >  inline
QP_full_exact_pricing<Q,ET,Tags>::
QP_full_exact_pricing()
  : QP_pricing_strategy<Q,ET,Tags>("full exact")
{ }
    
// operations
template < typename Q, typename ET, typename Tags >
int  QP_full_exact_pricing<Q,ET,Tags>::
pricing(int& direction )
{
  return (pricing_helper(direction, Is_nonnegative()));
}

template < typename Q, typename ET, typename Tags >
int  QP_full_exact_pricing<Q,ET,Tags>::
pricing_helper(int& /*direction*/, Tag_true /*is_in_standard_form*/)
{
  // get properties of quadratic program:
  int  w = this->solver().number_of_working_variables();

  // loop over all non-basic variables:
  int  j,  min_j  = -1;
  ET   mu, min_mu = this->et0;
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

      // new minimum?
      if (mu < min_mu) { min_j = j; min_mu = mu; }
    }
  }
  CGAL_qpe_debug { 
    this->vout() << std::endl;
  }

  // return index of entering variable:
  return min_j;
    
}

template < typename Q, typename ET, typename Tags >
int  QP_full_exact_pricing<Q,ET,Tags>::
pricing_helper(int& direction, Tag_false /*is_in_standard_form*/)
{
    
  // get properties of quadratic program:
  int  w = this->solver().number_of_working_variables();

  // loop over all non-basic variables:
  int  j,  min_j  = -1;
  // 
  ET   min_mu = this->et0;     // Note: for mu_j > 0 we will compare -mu_j and
			       // min_mu.
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
    }
  }
  CGAL_qpe_debug { 
    this->vout() << std::endl;
  }

  // return index of entering variable
  return min_j;
    
}


} //namespace CGAL

#endif // CGAL_QP_FULL_EXACT_PRICING_H

// ===== EOF ==================================================================
