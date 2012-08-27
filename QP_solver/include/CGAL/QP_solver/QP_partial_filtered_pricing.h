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

#ifndef CGAL_QP_PARTIAL_FILTERED_PRICING_H
#define CGAL_QP_PARTIAL_FILTERED_PRICING_H

// MSVC detection
#include <boost/config.hpp>

// includes
#include <CGAL/QP_solver/QP__partial_base.h>
#include <CGAL/QP_solver/QP__filtered_base.h>


// MSVC complains about inheritance through dominance when only one
// base implements virtual functions from the top of the diamond.
#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4250)
#endif

namespace CGAL {

// =================
// class declaration
// =================
template < typename Q, typename ET, typename Tags, class NT_ = double, class ET2NT_ =
    To_double<ET> >
class QP_partial_filtered_pricing;

// ===============
// class interface
// ===============
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
class QP_partial_filtered_pricing
    : public QP__partial_base <Q,ET,Tags>,
      public QP__filtered_base<Q,ET,Tags,NT_,ET2NT_> {

    // self
    typedef  QP_pricing_strategy<Q,ET,Tags>          Base;
    typedef  QP__partial_base<Q,ET,Tags>             Partial_base;
    typedef  QP__filtered_base<Q,ET,Tags, NT_, ET2NT_>  Filtered_base;
    typedef  QP_partial_filtered_pricing<Q, ET, Tags, NT_, ET2NT_>  Self;

    // types from the base class
    typedef  typename Tags::Is_nonnegative           Is_nonnegative;
    typedef  typename Partial_base::Index_iterator        Index_iterator;
    typedef  typename Partial_base::Index_const_iterator  Index_const_iterator;

    using Base::price_dantzig;
    using Base::is_improving;

  public:

    // number type
    typedef  NT_                        NT;
    typedef  ET2NT_                     ET2NT;

    // creation
    QP_partial_filtered_pricing( bool     randomize = false,
				  Random&  random    = default_random,
				  ET2NT    et2nt     = ET2NT());

    // operations
    int  pricing(int& direction );

    void  init( );
    void  transition( );
    
    
    // cleanup
    ~QP_partial_filtered_pricing() {};

  private:
  int pricing_helper(int& direction, Tag_true  /*is_in_standard_form*/);
  int pricing_helper(int& direction, Tag_false /*is_in_standard_form*/);
};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// construction
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >  inline
QP_partial_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
QP_partial_filtered_pricing( bool randomize, Random& random, ET2NT et2nt)
    : Base( "partial filtered"),
      Partial_base( randomize, random),
      Filtered_base( et2nt)
{ }
    
// operations
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >  inline
void  QP_partial_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
init( )
{
     Partial_base::init();
    Filtered_base::init();
}

template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >  inline
void  QP_partial_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
transition( )
{
     Partial_base::transition();
    Filtered_base::transition();
}


template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
int  QP_partial_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
pricing(int& direction ) 
{
  return (pricing_helper(direction, Is_nonnegative()));
}

template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
int  QP_partial_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
pricing_helper(int& /*direction*/, Tag_true /*is_in_standard_form*/ )
{
    // initialize filtered computation
    this->init_NT();

    // loop over all active non-basic variables
    CGAL_qpe_debug {
	this->vout() << "active variables:" << std::endl;
    }

    Index_const_iterator  it, min_it;
    NT                            mu, min_mu = this->nt0;
    for ( it = this->active_set_begin(); it != this->active_set_end(); ++it) {

	// don't price artificial variables
	if (this->solver().is_artificial( *it) ||
	    this->solver().is_basic( *it))  // added by kf
	  continue;

	// compute mu_j
	mu = this->mu_j_NT( *it);

	CGAL_qpe_debug {
	    this->vout() << "  mu_" << *it << " [NT]: " << mu << std::endl;
	}

	// new minimum?
	if ( mu < min_mu) { min_it = it; min_mu = mu; }
    }

    // exact check of entering variable
    if ( min_mu < this->nt0) {
	if ( this->mu_j( *min_it) >= this->et0) {

	    // exact check failed!
	    CGAL_qpe_debug {
		this->vout() << "--> exact check of entering variable failed!"
		       << std::endl;
	    }

	    // reject entering variable
	    min_mu = this->nt0;
	}
    } else {
	CGAL_qpe_debug {
	    this->vout() << "--> no entering variable found yet" << std::endl;
	}
    }

    // no entering variable found so far?
    if ( ( min_mu == this->nt0) &&
         ( this->inactive_set_begin() < this->inactive_set_end())) {

	// loop over all inactive non-basic variables
	CGAL_qpe_debug {
	    this->vout() << "inactive variables:" << std::endl;
	}
	Index_const_iterator  active_it;
	for ( it = this->inactive_set_begin(); it != this->inactive_set_end(); ++it) {

	    // don't price artificial variables
	    if (this->solver().is_artificial( *it)) continue;

	    // compute mu_j
	    mu = this->mu_j_NT( *it);

	    CGAL_qpe_debug {
		this->vout() << "  mu_" << *it << " [NT]: " << mu << std::endl;
	    }

	    // candidate for entering?
	    if ( mu < this->nt0) {

		// make variable active
		active_it = it;
		this->activating( active_it);

		// new minimum?
		if ( mu < min_mu) { min_it = active_it; min_mu = mu; }
	    }
	}

	// exact check of entering variable
	if ( min_mu < this->nt0) {
	    if ( this->mu_j( *min_it) >= this->et0) {

		// exact check failed!
		CGAL_qpe_debug {
		    this->vout() << "--> exact check of entering variable failed!"
			   << std::endl;
		}

		// reject entering variable
		min_mu = this->nt0;
	    }
	} else {
	    CGAL_qpe_debug {
		this->vout() << "--> still no entering variable found" << std::endl;
	    }
	}
    }

    // certify non-existance of entering variable, if necessary
    if ( min_mu == this->nt0) {

	// update row and column maxima
	this->update_maxima();

	// loop over all non-basic variables again
	for ( it = this->active_set_begin(); it != this->inactive_set_end(); ++it) {

	    // don't price artificial variables
	    if (this->solver().is_artificial( *it)) continue;

	    // certify 'mu_j >= 0'
	    if ( ! this->certify_mu_j_NT( *it)) {

		// entering variable missed by inexact arithmetic
		min_it = it;
		break;
	    }
	}
    }
    CGAL_qpe_debug { 
      this->vout() << std::endl;
    }

    // return index of entering variable, if any
    if ( min_mu < this->nt0) {
	int  j = *min_it;
        this->entering_basis( min_it);
	return j;
    }

    // no entering variable found
    return -1;
}
template < typename Q, typename ET, typename Tags, class NT_, class ET2NT_ >
int  QP_partial_filtered_pricing<Q,ET,Tags,NT_,ET2NT_>::
pricing_helper(int& direction, Tag_false /*is_in_standard_form*/ )
{
    // initialize filtered computation
    this->init_NT();

    // loop over all active non-basic variables
    CGAL_qpe_debug {
	this->vout() << "active variables:" << std::endl;
    }

    Index_const_iterator  it, min_it;    
    int min_j = -1;
    NT  mu, min_mu = this->nt0;
    for ( it = this->active_set_begin(); it != this->active_set_end(); ++it) {

	// don't price artificial variables
	if (this->solver().is_artificial( *it) ||
	    this->solver().is_basic( *it))  // added by kf
	  continue;

	// compute mu_j
	mu = this->mu_j_NT( *it);

	if (price_dantzig (*it, mu, this->nt0, min_j, min_mu, direction))
	  min_it = it;
    }

    if ( min_j >= 0 ) {
        // exact check; do we really have an entering variable
	if ( !this->is_improving(min_j, this->mu_j( min_j), this->et0)) {

	    // exact check failed!
	    CGAL_qpe_debug {
		this->vout() << "--> exact check of entering variable failed!"
		       << std::endl;
	    }

	    // reject entering variable
	    min_j = -1;
	    min_mu = this->nt0;
	}
    } else {
	CGAL_qpe_debug {
	    this->vout() << "--> no entering variable found yet" << std::endl;
	}
    }

    // no entering variable found so far?
    if ( ( min_j == -1) &&
         ( this->inactive_set_begin() < this->inactive_set_end())) {

	// loop over all inactive non-basic variables
	CGAL_qpe_debug {
	    this->vout() << "inactive variables:" << std::endl;
	}
	Index_const_iterator  active_it;
	for ( it = this->inactive_set_begin(); 
	      it != this->inactive_set_end(); ++it) {

	    // don't price artificial variables
	    if (this->solver().is_artificial( *it)) continue;

	    // compute mu_j
	    mu = this->mu_j_NT( *it);

	    CGAL_qpe_debug {
		this->vout() << "  mu_" << *it << " [NT]: " << mu << std::endl;
	    }

	    // candidate for entering?
	    if (is_improving(*it, mu, this->nt0)) {

		// make variable active
		active_it = it;
		this->activating( active_it);
		
		// new minimum
		if (price_dantzig (*active_it, mu, this->nt0, 
			            min_j, min_mu, direction))
		  min_it = active_it;
	    }
	}

	if ( min_j >= 0) {	
	    // exact check of entering variable
	    if (!this->is_improving(min_j, this->mu_j( min_j), this->et0)) {

		// exact check failed!
		CGAL_qpe_debug {
		    this->vout() << 
		      "--> exact check of entering variable failed!"
		      << std::endl;
		}

		// reject entering variable
		min_j = -1;
		min_mu = this->nt0;
	    }
	} else {
	    CGAL_qpe_debug {
		this->vout() << 
		  "--> still no entering variable found" 
		  << std::endl;
	    }
	}
    }

    // certify non-existance of entering variable, if necessary
    if ( min_j == -1) {

	// update row and column maxima
	this->update_maxima();

	// loop over all non-basic variables again
	for ( it = this->active_set_begin(); 
	      it != this->inactive_set_end(); ++it) {

	    // don't price artificial variables
	    if (this->solver().is_artificial( *it)) continue;

	    if ( ! this->certify_mu_j_NT( *it)) {

		// entering variable missed by inexact arithmetic
	      min_j = *it;
	      min_it = it;
	      break;
	    }
	}
    }
    CGAL_qpe_debug { 
      this->vout() << std::endl;
    }

    // return index of entering variable, if any
    if ( min_j >= 0) {
      CGAL_qpe_assertion(min_j == *min_it);
      this->entering_basis( min_it);
      return min_j;
    }

    // no entering variable found
    return -1;
}

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_QP_PARTIAL_FILTERED_PRICING_H

// ===== EOF ==================================================================
