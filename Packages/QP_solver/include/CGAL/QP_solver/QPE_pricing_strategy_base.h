// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
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
// file          : include/CGAL/QP_engine/QPE_pricing_strategy_base.h
// package       : $CGAL_Package: QP_engine $
// chapter       : Generalized Linear Programming
//
// revision      : 3.0alpha
// revision_date : 2002/05/23
//
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Base Class for Pricing Strategies for the QPE Solver
// ============================================================================

#ifndef CGAL_QPE_PRICING_STRATEGY_BASE_H
#define CGAL_QPE_PRICING_STRATEGY_BASE_H

// includes
#include <CGAL/QPE_solver.h>
#include <CGAL/IO/Verbose_ostream.h>

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class Rep_ >
class QPE_pricing_strategy_base;

template < class Rep_ >
class QPE_solver;

// Class interface
// ===============
struct To_double {
    template < class NT >
    double  operator () ( const NT& x) { return CGAL::to_double( x); }
};

template < class Rep_ >
class QPE_pricing_strategy_base {
  public:
    // self
    typedef  Rep_                            Rep;
    typedef  QPE_pricing_strategy_base<Rep>  Self;

    // the ambient QP solver
    typedef  CGAL::QPE_solver<Rep>      QP_solver;

    // types from the ambient QP solver
    typedef  typename QP_solver::ET     ET;

    typedef  typename QP_solver::A_iterator
                                        A_iterator;
    typedef  typename QP_solver::B_iterator
                                        B_iterator;
    typedef  typename QP_solver::C_iterator
                                        C_iterator;
    typedef  typename QP_solver::D_iterator
                                        D_iterator;

    typedef  typename QP_solver::A_artificial_iterator
                                        A_artificial_iterator;
    typedef  typename QP_solver::C_auxiliary_iterator
                                        C_auxiliary_iterator;

    typedef  typename QP_solver::Basic_variable_index_iterator
                                        Basic_variable_index_iterator;
    typedef  typename QP_solver::Basic_variable_value_iterator
                                        Basic_variable_value_iterator;
    typedef  typename QP_solver::Basic_variable_numerator_iterator
                                        Basic_variable_numerator_iterator;

    typedef  typename QP_solver::Lambda_value_iterator
                                        Lambda_value_iterator;
    typedef  typename QP_solver::Lambda_numerator_iterator
                                        Lambda_numerator_iterator;

    // types from the representation class
    typedef  typename Rep::Is_linear    Is_linear;
    typedef  typename Rep::Is_symmetric Is_symmetric;
    typedef  typename Rep::Has_no_inequalities
                                        Has_no_inequalities;

  protected:
    // protected types
    typedef  CGAL::Tag_true             Tag_true;
    typedef  CGAL::Tag_false            Tag_false;

  public:
    // creation
    QPE_pricing_strategy_base( )
        : et0( 0), et1( 1), et2( 2)
    { }
    
    // destruction
    virtual ~QPE_pricing_strategy_base( ) { }
    
    // initialization
    void  set( const QP_solver&        qp_solver,
               CGAL::Verbose_ostream&  verbose_out)
        {
            solverP = &qp_solver;
            voutP   = &verbose_out;
            set();
        }
    
    virtual  void  set( ) { }
    
    virtual  void  init( ) { }
    
    // access
    const QP_solver&        solver( ) const { return *solverP; }
    CGAL::Verbose_ostream&  vout  ( ) const { return *voutP; }
    
    // operations
    virtual  int   pricing( ) = 0;
    
    virtual  void  leaving_basis( int) { }
    
    virtual  void  transition( ) { }
    
    ET  mu_j( int j) const
    {
        // move this function from solver to here
        return et0;
    }
    

  private:
    // some constants
    const ET  et0, et1, et2;

    // data members
    const QP_solver*         solverP;   // the ambient QP solver
    CGAL::Verbose_ostream*   voutP;     // used for verbose output
};

CGAL_END_NAMESPACE

#endif // CGAL_QPE_PRICING_STRATEGY_BASE_H

// ===== EOF ==================================================================
