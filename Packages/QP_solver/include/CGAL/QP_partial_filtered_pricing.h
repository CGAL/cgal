// ============================================================================
//
// Copyright (c) 1997-2004 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/QP_engine/QPE_partial_filtered_pricing.h
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Pricing Strategy with partial filtered pricing
// ============================================================================

#ifndef CGAL_QPE_PARTIAL_FILTERED_PRICING_H
#define CGAL_QPE_PARTIAL_FILTERED_PRICING_H

// includes
#include <CGAL/QP_engine/QPE__partial_base.h>
#include <CGAL/QP_engine/QPE__filtered_base.h>

CGAL_BEGIN_NAMESPACE

// =================
// class declaration
// =================
template < class Rep_, class NT_ = double, class ET2NT_ = To_Double >
class QPE_partial_filtered_pricing;

// ===============
// class interface
// ===============
template < class Rep_, class NT_, class ET2NT_ >
class QPE_partial_filtered_pricing
    : public QPE__partial_base <Rep_>,
      public QPE__filtered_base<Rep_,NT_,ET2NT_> {

    // self
    typedef  Rep_                               Rep;
    typedef  QPE_pricing_strategy<Rep>          Base;
    typedef  QPE__partial_base<Rep>             Partial_base;
    typedef  QPE__filtered_base<Rep>            Filtered_base;
    typedef  QPE_partial_filtered_pricing<Rep>  Self;

    // types from the base class
    typedef  typename Base::ET                            ET;
    typedef  typename Partial_base::Index_iterator        Index_iterator;
    typedef  typename Partial_base::Index_const_iterator  Index_const_iterator;
  public:

    // number type
    typedef  NT_                        NT;
    typedef  ET2NT_                     ET2NT;

    // creation
    QPE_partial_filtered_pricing( bool     randomize = false,
				  Random&  random    = default_random,
				  ET2NT    et2nt     = ET2NT());

    // operations
    int  pricing( );

    void  init( );
    void  transition( );
};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// construction
template < class Rep_, class NT_, class ET2NT_ >  inline
QPE_partial_filtered_pricing<Rep_,NT_,ET2NT_>::
QPE_partial_filtered_pricing( bool randomize, Random& random, ET2NT et2nt)
    : Base( "partial filtered"),
      Partial_base( randomize, random),
      Filtered_base( et2nt)
{ }
    
// operations
template < class Rep_, class NT_, class ET2NT_ >  inline
void  QPE_partial_filtered_pricing<Rep_,NT_,ET2NT_>::
init( )
{
     Partial_base::init();
    Filtered_base::init();
}

template < class Rep_, class NT_, class ET2NT_ >  inline
void  QPE_partial_filtered_pricing<Rep_,NT_,ET2NT_>::
transition( )
{
     Partial_base::transition();
    Filtered_base::transition();
}

template < class Rep_, class NT_, class ET2NT_ >
int  QPE_partial_filtered_pricing<Rep_,NT_,ET2NT_>::
pricing( )
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

	// compute mu_j
	mu = mu_j_NT( *it);

	CGAL_qpe_debug {
	    this->vout() << "  mu_" << *it << " [NT]: " << mu << std::endl;
	}

	// new minimum?
	if ( mu < min_mu) { min_it = it; min_mu = mu; }
    }

    // exact check of entering variable
    if ( min_mu < this->nt0) {
	if ( mu_j( *min_it) >= this->et0) {

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

	    // compute mu_j
	    mu = mu_j_NT( *it);

	    CGAL_qpe_debug {
		this->vout() << "  mu_" << *it << " [NT]: " << mu << std::endl;
	    }

	    // candidate for entering?
	    if ( mu < this->nt0) {

		// make variable active
		active_it = it;
		activating( active_it);

		// new minimum?
		if ( mu < min_mu) { min_it = active_it; min_mu = mu; }
	    }
	}

	// exact check of entering variable
	if ( min_mu < this->nt0) {
	    if ( mu_j( *min_it) >= this->et0) {

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

	    // certify 'mu_j >= 0'
	    if ( ! certify_mu_j_NT( *it)) {

		// entering variable missed by inexact arithmetic
		min_it = it;
		break;
	    }
	}
    }
    this->vout() << std::endl;

    // return index of entering variable, if any
    if ( min_mu < this->nt0) {
	int  j = *min_it;
	entering_basis( min_it);
	return j;
    }

    // no entering variable found
    return -1;
}

CGAL_END_NAMESPACE

#endif // CGAL_QPE_PARTIAL_FILTERED_PRICING_H

// ===== EOF ==================================================================
