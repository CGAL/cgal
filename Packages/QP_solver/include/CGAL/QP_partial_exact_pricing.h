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
// file          : include/CGAL/QP_engine/QPE_partial_exact_pricing.h
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Pricing Strategy with partial exact pricing
// ============================================================================

#ifndef CGAL_QPE_PARTIAL_EXACT_PRICING_H
#define CGAL_QPE_PARTIAL_EXACT_PRICING_H

// includes
#include <CGAL/QP_engine/QPE__partial_base.h>

CGAL_BEGIN_NAMESPACE

// =================
// class declaration
// =================
template < class Rep_ >
class QPE_partial_exact_pricing;

// ===============
// class interface
// ===============
template < class Rep_ >
class QPE_partial_exact_pricing : public QPE__partial_base<Rep_> {

    // self
    typedef  Rep_                            Rep;
    typedef  QPE_pricing_strategy<Rep>       Base;
    typedef  QPE__partial_base<Rep>          Partial_base;
    typedef  QPE_partial_exact_pricing<Rep>  Self;

    // types from the pricing base class
    typedef  typename Base::ET                            ET;
    typedef  typename Partial_base::Index_iterator        Index_iterator;
    typedef  typename Partial_base::Index_const_iterator  Index_const_iterator;

  public:

    // creation
    QPE_partial_exact_pricing( bool     randomize = false,
			       Random&  random    = default_random);

    // operations
    int  pricing( );
    
    // creation
    ~QPE_partial_exact_pricing(){ };
};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// construction
template < class Rep_ >  inline
QPE_partial_exact_pricing<Rep_>::
QPE_partial_exact_pricing( bool  randomize, Random&  random)
    : Base( "partial exact"),
      Partial_base( randomize, random)
{ }
    
// operations
template < class Rep_ >
int  QPE_partial_exact_pricing<Rep_>::
pricing( )
{
    Index_const_iterator  it, min_it;
    ET                            mu, min_mu =  0;

    // loop over all active non-basic variables
    CGAL_qpe_debug {
	this->vout() << "active variables:" << std::endl;
    }
    for ( it = this->active_set_begin(); it != this->active_set_end(); ++it) {

        // don't price artificial variables
	if (this->solver().is_artificial( *it)) continue;

	// compute mu_j
	mu = mu_j( *it);

	CGAL_qpe_debug {
	    this->vout() << "  mu_" << *it << ": " << mu << std::endl;
	}

	// new minimum?
	if ( mu < min_mu) { min_it = it; min_mu = mu; }
    }

    // no entering variable found so far?
    if ( ( min_mu == this->et0) && ( this->inactive_set_begin() <
                                     this->inactive_set_end())) {

	// loop over all inactive non-basic variables
	CGAL_qpe_debug {
	    this->vout() << "inactive variables:" << std::endl;
	}
	Index_const_iterator  active_it;
	for ( it = this->inactive_set_begin(); it != this->inactive_set_end(); ++it) {

	    // don't price artificial variables
	    if (this->solver().is_artificial( *it)) continue;
	    
	    // compute mu_j
	    mu = mu_j( *it);

	    CGAL_qpe_debug {
		this->vout() << "  mu_" << *it << ": " << mu << std::endl;
	    }

	    // candidate for entering?
	    if ( mu < this->et0) {

		// make variable active
		active_it = it;
		activating( active_it);

		// new minimum?
		if ( mu < min_mu) { min_it = active_it; min_mu = mu; }
	    }
	}
    }
    this->vout() << std::endl;

    // return index of entering variable, if any
    if ( min_mu < this->et0) {
	int  j = *min_it;
	entering_basis( min_it);
	return j;
    }

    // no entering variable found
    return -1;
}

CGAL_END_NAMESPACE

#endif // CGAL_QPE_PARTIAL_EXACT_PRICING_H

// ===== EOF ==================================================================
