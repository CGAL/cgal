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
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp 
//                 Kaspar Fischer

#ifndef CGAL_QP_FULL_FILTERED_PRICING_H
#define CGAL_QP_FULL_FILTERED_PRICING_H

// includes
#include <CGAL/QP_solver/QP__filtered_base.h>

namespace CGAL {

// =================
// class declaration
// =================
template < typename Q, typename ET, typename Tags, class NT_ = double, class ET2NT_ =
    To_double<ET> >
class QP_full_filtered_pricing;

// ===============
// class interface
// ===============
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
class QP_full_filtered_pricing : public QP__filtered_base<Q,ET,Tags,NT_,ET2NT_> {

    // self
    typedef  QP_pricing_strategy<Q, ET, Tags>       Base;
    typedef  typename Tags::Is_nonnegative     Is_nonnegative;
    typedef  QP__filtered_base<Q, ET, Tags, NT_, ET2NT_>         Filtered_base;
    typedef  QP_full_filtered_pricing<Q, ET, Tags, NT_, ET2NT_>  Self;

  public:

    // number type
    typedef  NT_                        NT;
    typedef  ET2NT_                     ET2NT;

    // creation
    QP_full_filtered_pricing( ET2NT et2nt = ET2NT());

    // operations
    int  pricing(int& direction );
    
    // cleanup
    ~QP_full_filtered_pricing() { };

  private:
    int pricing_helper(int& direction, Tag_true  is_in_standard_form);
    int pricing_helper(int& direction, Tag_false is_in_standard_form);

};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// construction
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >  inline
QP_full_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
QP_full_filtered_pricing( ET2NT et2nt)
    : Base( "full filtered"),
      Filtered_base( et2nt)
{ }
    
// operations
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
int  QP_full_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
pricing (int& direction) 
{
  return (pricing_helper(direction, Is_nonnegative()));
}

template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
int  QP_full_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
pricing_helper(int& /*direction*/, Tag_true ) // standard form
{
    // get properties of quadratic program
    int  w = this->solver().number_of_working_variables();

    // initialize filtered computation
    this->init_NT();

    // loop over all non-basic variables
    int  j,  min_j  = -1;
    NT   mu, min_mu = this->nt0;
    for ( j = 0; j < w; ++j) {

	// variable non-basic?
	if ( ! this->solver().is_basic( j)) {
	
	    // don't price artificial variables
	    if (this->solver().is_artificial( j)) continue;


	    // compute mu_j
	    mu = this->mu_j_NT( j);

	    CGAL_qpe_debug {
		this->vout() << "mu_" << j << " [NT]: " << mu << std::endl;
	    }

	    // new minimum?
	    if ( mu < min_mu) { min_j = j; min_mu = mu; }
	}
    }

    // exact check of entering variable
    if ( min_j >= 0) {
	if ( this->mu_j( min_j) >= this->et0) {

	    // exact check failed!
	    CGAL_qpe_debug {
		this->vout() << "--> exact check of entering variable failed!"
		       << std::endl;
	    }
	    
	    min_j  = -1;
	    min_mu = this->nt0;
	}
    }

    // certify non-existance of entering variable, if necessary
    if ( min_j < 0) {

	// update row and column maxima
	this->update_maxima();

	// loop over all non-basic variables again
	for ( j = 0; j < w; ++j) {

	    // variable non-basic?
	    if ( ! this->solver().is_basic( j)) {
	    
	        // don't price artificial variables
	        if (this->solver().is_artificial( j)) continue;


		// certify 'mu_j >= 0'
		if ( ! this->certify_mu_j_NT( j)) {

		    // entering variable missed by inexact arithmetic
		    min_j = j;
		    break;
		}
	    }
	}
    }
    CGAL_qpe_debug { 
      this->vout() << std::endl;
    }
    // return index of entering variable
    return min_j;
}

template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
int  QP_full_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
pricing_helper(int& direction, Tag_false ) // bounds for variables
{
    // get properties of quadratic program
    int  w = this->solver().number_of_working_variables();

    // initialize filtered computation
    this->init_NT();

    // loop over all non-basic variables
    int  j,  min_j  = -1;
    NT   mu, min_mu = this->nt0;
    for ( j = 0; j < w; ++j) {

	// variable non-basic?
	if ( ! this->solver().is_basic( j)) {
	
	    // don't price artificial variables
	    if (this->solver().is_artificial( j)) continue;


	    // compute mu_j
	    mu = this->mu_j_NT( j);

	    CGAL_qpe_debug {
		this->vout() << "mu_" << j << " [NT]: " << mu << std::endl;
	    }
	    // from pricing strategy base class
	    this->price_dantzig (j, mu, this->nt0, min_j, min_mu, direction);
	}
    }

    if ( min_j >= 0) {
        // exact check; do we really have an entering variable
	if (!this->is_improving(min_j, this->mu_j( min_j), this->et0)) {

	    // exact check failed!
	    CGAL_qpe_debug {
		this->vout() << "--> exact check of entering variable failed!"
		       << std::endl;
	    }
	    
	    min_j  = -1;
	    min_mu = this->nt0;
	}
    }

    if ( min_j == -1) {
        // try to certify non-existence of entering variable, based on
        // error bounds
	this->update_maxima();

	// loop over all non-basic variables again
	for ( j = 0; j < w; ++j) {

	    // variable non-basic?
	    if ( ! this->solver().is_basic( j)) {
	    
	        // don't price artificial variables
	        if (this->solver().is_artificial( j)) continue;

		// certify that j is not improving
		if ( ! this->certify_mu_j_NT( j)) {

		    // entering variable missed by inexact arithmetic
		    min_j = j;
		    break;
		}
	    }
	}
    }
    CGAL_qpe_debug { 
      this->vout() << std::endl;
    }

    // return index of entering variable
    return min_j;
}

} //namespace CGAL

#endif // CGAL_QP_FULL_FILTERED_PRICING_H

// ===== EOF ==================================================================
