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

#ifndef CGAL_QP_PRICING_STRATEGY_H
#define CGAL_QP_PRICING_STRATEGY_H

#include <CGAL/license/QP_solver.h>


// includes
#include <CGAL/QP_solver/QP_solver.h>
#include <CGAL/IO/Verbose_ostream.h>

#include <string>

namespace CGAL {

// ==================
// class declarations
// ==================
template < typename Q, typename ET, typename Tags >
class QP_pricing_strategy;

template < typename Q, typename ET, typename Tags >
class QP_solver;

// ===============
// class interface
// ===============
template < typename Q, typename ET, typename Tags >
class QP_pricing_strategy {

  public:

    // self
    typedef  QP_pricing_strategy<Q, ET, Tags>  Self;

    // types
    typedef  CGAL::QP_solver<Q, ET, Tags>      QP_solver;
    typedef  CGAL::Verbose_ostream      Verbose_ostream;
    typedef  typename Tags::Is_nonnegative
                                        Is_nonnegative;
    typedef  typename Tags::Is_linear    Is_linear;
    typedef typename QP_solver::Bound_index Bound_index;

  public:

    // initialization
    void  set ( const QP_solver& solver, Verbose_ostream& vout);
    void  init( int dummy);

    // operations
    virtual  int   pricing(int& /*direction*/ ) = 0;

    virtual  void  leaving_basis( int /*i*/) { }
    virtual  void  transition( ) { }
    
  protected:
    
    // construction & destruction
    QP_pricing_strategy( const std::string& strategy_name);
public:
    virtual ~QP_pricing_strategy( ) { }
protected:
    QP_pricing_strategy( );            // detects error in virtual inheritance
        
    // initialization (of derived classes)
    virtual  void  set ( ) { }
    virtual  void  init( ) { }
    
    // operations
    ET  mu_j( int j) const;

    // access
    const QP_solver&  solver( ) const { return *solverP; }
    Verbose_ostream&  vout  ( ) const { return *voutP;   }

    // constants (protected)
    const ET  et0;

    // used during pricing loop; finds out whether j is a
    // new best candidate for the entering variable w.r.t.
    // Dantzig's pivot rule (reduced cost pricing); the
    // return value is true iff j is the new best candidate
    template <typename NT>
    bool price_dantzig (int j, const NT& mu, const NT& nt0,
		     int& min_j, NT& min_mu, int& direction);

    // returns whether j satisfying mu_j = mu is a candidate
    // for the entering variable
    template <typename NT>
    bool is_improving (int j, const NT& mu, const NT& nt0) const;

  private:

    // data members
    const QP_solver*  solverP;          // the ambient QP solver
    Verbose_ostream*  voutP;            // used for verbose output
    std::string       name;             // derived strategy's name
};

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// construction
template < typename Q, typename ET, typename Tags >  inline
QP_pricing_strategy<Q, ET, Tags>::
QP_pricing_strategy( const std::string& strategy_name)
    : et0( 0), name( strategy_name)
{ }

// detects error in virtual inheritance
template < typename Q, typename ET, typename Tags >  inline
QP_pricing_strategy<Q, ET, Tags>::
QP_pricing_strategy( )
  : et0(0)
{
    CGAL_qpe_assertion_msg
      ( false, "call to 'QP_pricing_strategy<Q,ET,Tags>::\n'" \
	"QP_pricing_strategy( const std::string&  strategy_name)'\n" \
	"is missing in most derived pricing class!");
}

// initialization
template < typename Q, typename ET, typename Tags >  inline
void  QP_pricing_strategy<Q, ET, Tags>::
set( const QP_solver&  solver, Verbose_ostream&  vout)
{
    solverP = &solver;
    voutP   = &vout;
    set();
}

template < typename Q, typename ET, typename Tags >  inline
void  QP_pricing_strategy<Q, ET, Tags>::
init( int)
{
    CGAL_qpe_debug {
	vout() << "pricing: " << name << std::endl;
    }
    init();
}

// operations
template < typename Q, typename ET, typename Tags >  inline
ET  QP_pricing_strategy<Q, ET, Tags>::
mu_j( int j) const
{
  return this->solver().mu_j(j);
}

template < typename Q, typename ET, typename Tags >
template <typename NT> 
bool  QP_pricing_strategy<Q, ET, Tags>::
is_improving (int j, const NT& mu, const NT& nt0 ) const
{  
  CGAL_qpe_assertion(!this->solver().is_basic(j));
  CGAL_qpe_assertion(!this->solver().is_artificial(j));
  if (this->solver().is_original(j)) {
    const Bound_index bnd_ind =
      this->solver().nonbasic_original_variable_bound_index(j);
    switch (bnd_ind) {
    case QP_solver::LOWER:
      return mu < nt0;
    case QP_solver::UPPER:
      return mu > nt0;
    case QP_solver::ZERO:
      {
	const int where =
	  this->solver().state_of_zero_nonbasic_variable(j);
	return ((where >= 0 && mu > nt0) || (where <= 0 && mu < nt0));
      }
    case QP_solver::FIXED:
      return false;
    case QP_solver::BASIC:
    default:
      CGAL_qpe_assertion(false);
      return false;
    }
  } else {
    // artficial variable
    return mu < nt0;
  }
}

template < typename Q, typename ET, typename Tags >
template <typename NT> 
bool QP_pricing_strategy<Q, ET, Tags>::
price_dantzig (int j, const NT& mu, const NT& nt0,
	 int& min_j, NT& min_mu, int& direction) {
  CGAL_qpe_assertion(!this->solver().is_basic(j));
  CGAL_qpe_assertion(!this->solver().is_artificial(j));
  if (this->solver().is_original(j)) {
    // original variable
    const Bound_index bnd_ind =
      this->solver().nonbasic_original_variable_bound_index(j);
    switch (bnd_ind) {
    case QP_solver::LOWER:
      {
	CGAL_qpe_debug { 
	  this->vout() << "mu_" << j << ": " << mu
		       << " LOWER" << std::endl;
	}
	    
	if (mu < nt0) {
	  // new minimum?
	  if (mu < min_mu) {
	    min_j = j; min_mu = mu;
	    direction = 1;
	  }
	}
	break;
      }
    case QP_solver::ZERO:
      {
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
	    mu > nt0) {
	  // new minimum?
	  if (-mu < min_mu) {
	    min_j = j; min_mu = -mu;
	    direction = -1;
	  }                                                    
	}
	if (where <= 0 &&       // middle or on lower bound?
	    mu < nt0) {
	  // new minimum?
	  if (mu < min_mu) {
	    min_j = j; min_mu = mu;
	    direction = 1;
	  }                            
	}
	break;
      }
    case QP_solver::UPPER:
      {
	CGAL_qpe_debug { 
	  this->vout() << "mu_" << j << ": " << mu
		       << " UPPER" << std::endl;
	}
	    
	if (mu > nt0) {
	  // new minimum?
	  if (-mu < min_mu) {
	    min_j = j; min_mu = -mu;
	    direction = -1;
	  }
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
  } else {                                    
    // slack variable
    CGAL_qpe_debug {
      this->vout() << "mu_" << j << ": " << mu 
		   << " LOWER (slack)" << std::endl;
    }

    // new minimum?
    if (mu < min_mu) {
      min_j = j; min_mu = mu;
      direction = 1;
    }
  }
  return (min_j == j);
}



} //namespace CGAL

#endif // CGAL_QP_PRICING_STRATEGY_H

// ===== EOF ==================================================================
