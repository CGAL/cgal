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
// file          : include/CGAL/QP_solver/QP_full_filtered_pricing.h
// package       : $CGAL_Package: QP_solver $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Pricing Strategy with full filtered pricing
// ============================================================================

#ifndef CGAL_QP_FULL_FILTERED_PRICING_H
#define CGAL_QP_FULL_FILTERED_PRICING_H

// includes
#include <CGAL/QP_solver/QP__filtered_base.h>

CGAL_BEGIN_NAMESPACE

// =================
// class declaration
// =================
template < class Rep_, class NT_ = double, class ET2NT_ =
    To_double<typename Rep_::ET> >
class QP_full_filtered_pricing;

// ===============
// class interface
// ===============
template < class Rep_, class NT_, class ET2NT_ >
class QP_full_filtered_pricing : public QP__filtered_base<Rep_,NT_,ET2NT_> {

    // self
    typedef  Rep_                            Rep;
    typedef  QP_pricing_strategy<Rep>       Base;
    typedef  QP__filtered_base<Rep, NT_, ET2NT_>         Filtered_base;
    typedef  QP_full_filtered_pricing<Rep, NT_, ET2NT_>  Self;

    // types from the base class
    typedef  typename Base::ET          ET;

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

};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// construction
template < class Rep_, class NT_, class ET2NT_ >  inline
QP_full_filtered_pricing<Rep_,NT_,ET2NT_>::
QP_full_filtered_pricing( ET2NT et2nt)
    : Base( "full filtered"),
      Filtered_base( et2nt)
{ }
    
// operations
template < class Rep_, class NT_, class ET2NT_ >
int  QP_full_filtered_pricing<Rep_,NT_,ET2NT_>::
pricing(int& direction )
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
    this->vout() << std::endl;

    // return index of entering variable
    return min_j;
}

CGAL_END_NAMESPACE

#endif // CGAL_QP_FULL_FILTERED_PRICING_H

// ===== EOF ==================================================================
