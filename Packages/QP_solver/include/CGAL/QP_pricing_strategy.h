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
// file          : include/CGAL/QP_pricing_strategy.h
// package       : $CGAL_Package: QP_solver $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Base Class for Pricing Strategies of the QP Solver
// ============================================================================

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

  public:

    // initialization
    void  set ( const QP_solver& solver, Verbose_ostream& vout);
    void  init( int dummy);

    // operations
    virtual  int   pricing( ) = 0;

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
    return solver().mu_j( j,
			  solver().lambda_numerator_begin(),
			  solver().basic_original_variables_numerator_begin(),
			  solver().variables_common_denominator());
}

CGAL_END_NAMESPACE

#endif // CGAL_QP_PRICING_STRATEGY_H

// ===== EOF ==================================================================
