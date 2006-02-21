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
// $URL$
// $Id$
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp <fransw@inf.ethz.ch>
//                 Kaspar Fischer <fischerk@inf.ethz.ch>

#ifndef CGAL_QP_PRICING_STRATEGY_H
#define CGAL_QP_PRICING_STRATEGY_H

// includes
#include <CGAL/QP_solver.h>
#include <CGAL/IO/Verbose_ostream.h>

#include <string>

CGAL_BEGIN_NAMESPACE

// ==================
// class declarations
// ==================
template < class Rep_ >
class QP_pricing_strategy;

template < class Rep_ >
class QP_solver;

// ===============
// class interface
// ===============
template < class Rep_ >
class QP_pricing_strategy {

  public:

    // self
    typedef  Rep_                       Rep;
    typedef  QP_pricing_strategy<Rep>  Self;

    // types
    typedef  typename Rep::ET           ET;
    typedef  CGAL::QP_solver<Rep>      QP_solver;
    typedef  CGAL::Verbose_ostream      Verbose_ostream;
    typedef  typename Rep::Is_in_standard_form
                                        Is_in_standard_form;
    typedef  typename Rep::Is_linear    Is_linear;


  public:

    // initialization
    void  set ( const QP_solver& solver, Verbose_ostream& vout);
    void  init( int dummy);

    // operations
    virtual  int   pricing(int& direction ) = 0;

    virtual  void  leaving_basis( int i) { }
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
template < class Rep_ >  inline
QP_pricing_strategy<Rep_>::
QP_pricing_strategy( const std::string& strategy_name)
    : et0( 0), name( strategy_name)
{ }

// detects error in virtual inheritance
template < class Rep_ >  inline
QP_pricing_strategy<Rep_>::
QP_pricing_strategy( )
  : et0(0)
{
    CGAL_qpe_assertion_msg( false, "call to 'QP_pricing_strategy<Rep>::\n'" \
	"QP_pricing_strategy( const std::string&  strategy_name)'\n" \
	"is missing in most derived pricing class!");
}

// initialization
template < class Rep_ >  inline
void  QP_pricing_strategy<Rep_>::
set( const QP_solver&  solver, Verbose_ostream&  vout)
{
    solverP = &solver;
    voutP   = &vout;
    set();
}

template < class Rep_ >  inline
void  QP_pricing_strategy<Rep_>::
init( int)
{
    CGAL_qpe_debug {
	vout() << "pricing: " << name << std::endl;
    }
    init();
}

// operations
template < class Rep_ >  inline
typename Rep_::ET  QP_pricing_strategy<Rep_>::
mu_j( int j) const
{
  return this->solver().mu_j(j);
}

CGAL_END_NAMESPACE

#endif // CGAL_QP_PRICING_STRATEGY_H

// ===== EOF ==================================================================
