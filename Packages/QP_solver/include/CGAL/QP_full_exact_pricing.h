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
// file          : include/CGAL/QP_full_exact_pricing.h
// package       : $CGAL_Package: QP_solver $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Pricing Strategy with full exact pricing
// ============================================================================

#ifndef CGAL_QP_FULL_EXACT_PRICING_H
#define CGAL_QP_FULL_EXACT_PRICING_H

// includes
#include <CGAL/QP_pricing_strategy.h>

CGAL_BEGIN_NAMESPACE

// =================
// class declaration
// =================
template < class Rep_ >
class QP_full_exact_pricing;

// ===============
// class interface
// ===============
template < class Rep_ >
class QP_full_exact_pricing : public QP_pricing_strategy<Rep_> {

    // self
    typedef  Rep_                         Rep;
    typedef  QP_pricing_strategy<Rep>    Base;
    typedef  QP_full_exact_pricing<Rep>  Self;

    // types from the base class
    typedef  typename Base::ET                      ET;
    typedef  typename Base::Is_in_standard_form     Is_in_standard_form;
    typedef  typename CGAL::QP_solver<Rep>          QP_solver;

  public:

    // creation
    QP_full_exact_pricing( );

    // operations
    int  pricing(int& direction );
    
    // cleanup
    ~QP_full_exact_pricing() { };
    
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
QP_full_exact_pricing<Rep_>::
QP_full_exact_pricing( )
    : QP_pricing_strategy<Rep_>( "full exact")
{ }
    
// operations
template < class Rep_ >
int  QP_full_exact_pricing<Rep_>::
pricing(int& direction )
{
    return (pricing_helper(direction, Is_in_standard_form()));
}

template < class Rep_ >
int  QP_full_exact_pricing<Rep_>::
pricing_helper(int& direction, Tag_true is_in_standard_form)
{
    // get properties of quadratic program
    int  w = this->solver().number_of_working_variables();

    // loop over all non-basic variables
    int  j,  min_j  = -1;
    ET   mu, min_mu = this->et0;
    for ( j = 0; j < w; ++j) {

        // variable non-basic?
        if ( ! this->solver().is_basic( j)) {
	
            // don't price artificial variables
            if (this->solver().is_artificial( j)) continue;

            // compute mu_j
            mu = this->mu_j( j);

            CGAL_qpe_debug {
                this->vout() << "mu_" << j << ": " << mu << std::endl;
            }

            // new minimum?
            if ( mu < min_mu) { min_j = j; min_mu = mu; }
        }
    }
    this->vout() << std::endl;

    // return index of entering variable
    return min_j;
    
}

template < class Rep_ >
int  QP_full_exact_pricing<Rep_>::
pricing_helper(int& direction, Tag_false is_in_standard_form)
{
    typedef typename QP_solver::Bound_index Bound_index;
    
    // get properties of quadratic program
    int  w = this->solver().number_of_working_variables();

    // loop over all non-basic variables
    int  j,  min_j  = -1;
    // Note: for mu_j > 0 we will compare -mu_j and min_mu,
    ET   mu, min_mu = this->et0;
    for ( j = 0; j < w; ++j) {

        // variable non-basic?
        if ( ! this->solver().is_basic( j)) {
	
            // don't price artificial variables
            if (this->solver().is_artificial( j)) continue;
            
            // original variable
            if (this->solver().is_original(j)) {
                Bound_index bnd_ind =
                    this->solver().nonbasic_original_variable_bound_index(j);
                switch (bnd_ind) {
                    case QP_solver::LOWER:
                        // compute mu_j
                        mu = this->mu_j(j);
                        
                        CGAL_qpe_debug { 
                            this->vout() << "mu_" << j << ": " << mu
                                << " LOWER" << std::endl;
                        }
                        
                        if (mu < this->et0) {
                            // new minimum?
                            if ( mu < min_mu) {
                                min_j = j; min_mu = mu;
                                direction = 1;
                            }
                        }
                        break;
                    case QP_solver::ZERO:
                        // compute mu_j
                        mu = this->mu_j(j);
                        
                        CGAL_qpe_debug { 
                            this->vout() << "mu_" << j << ": " << mu
                                << " ZERO" << std::endl;
                        }

                        if (mu < this->et0) {
                            // new minimum?
                            if ( mu < min_mu) {
                                min_j = j; min_mu = mu;
                                direction = 1;
                            }                            
                        } else if (mu > this->et0) {
                            // new minimum?
                            if ( -mu < min_mu) {
                                min_j = j; min_mu = -mu;
                                direction = -1;
                            }                                                    
                        }
                        break;
                    case QP_solver::UPPER:
                        // compute mu_j
                        mu = this->mu_j(j);
                        
                        CGAL_qpe_debug { 
                            this->vout() << "mu_" << j << ": " << mu
                                << " UPPER" << std::endl;
                        }
                        
                        if (mu > this->et0) {
                            // new minimum?
                            if ( -mu < min_mu) {
                                min_j = j; min_mu = -mu;
                                direction = -1;
                            }
                        }                    
                        break;
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
                mu = this->mu_j( j);

                CGAL_qpe_debug {
                    this->vout() << "mu_" << j << ": " << mu << std::endl;
                }

                // new minimum?
                if ( mu < min_mu) {
                    min_j = j; min_mu = mu;
                    direction = 1;
                }

            }
        }
    }
    this->vout() << std::endl;

    // return index of entering variable
    return min_j;
    
}


CGAL_END_NAMESPACE

#endif // CGAL_QP_FULL_EXACT_PRICING_H

// ===== EOF ==================================================================
