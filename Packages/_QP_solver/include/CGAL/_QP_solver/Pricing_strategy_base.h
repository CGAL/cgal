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
// file          : include/CGAL/_QP_solver/Pricing_strategy_base.h
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.4
// revision_date : 2000/08/17
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Base Class for Pricing Strategies for the QP Solver
// ============================================================================
                                                                               
#ifndef CGAL_PRICING_STRATEGY_BASE_H
#define CGAL_PRICING_STRATEGY_BASE_H

// includes
#include <CGAL/_QP_solver/QP_solver.h>
#include <CGAL/IO/Verbose_ostream.h>

CGAL_BEGIN_NAMESPACE
                    

// Class declaration
// =================
template < class _Rep >
class Pricing_strategy_base;

template < class _Rep >
class QP_solver;
                

// Class interface
// ===============
struct To_double {
    template < class NT >
    double  operator () ( const NT& x) { return CGAL::to_double( x); }
};

template < class _Rep >
class Pricing_strategy_base {
  public:
    // self
    typedef  _Rep                       Rep;
    typedef  Pricing_strategy_base<Rep> Self;

    // the ambient QP solver
    typedef  CGAL::QP_solver<Rep>       Solver;

    // types from the ambient QP solver
    typedef  typename Solver::NT        NT;
    typedef  typename Solver::ET        ET;

    typedef  typename Solver::A_iterator
                                        A_iterator;
    typedef  typename Solver::B_iterator
                                        B_iterator;
    typedef  typename Solver::C_iterator
                                        C_iterator;
    typedef  typename Solver::D_iterator
                                        D_iterator;

    typedef  typename Solver::A_artificial_iterator
                                        A_artificial_iterator;
    typedef  typename Solver::C_auxiliary_iterator
                                        C_auxiliary_iterator;

    typedef  typename Solver::Basic_variable_index_iterator
                                        Basic_variable_index_iterator;
    typedef  typename Solver::Basic_variable_value_iterator
                                        Basic_variable_value_iterator;
    typedef  typename Solver::Basic_variable_numerator_iterator
                                        Basic_variable_numerator_iterator;

    typedef  typename Solver::Lambda_value_iterator
                                        Lambda_value_iterator;
    typedef  typename Solver::Lambda_numerator_iterator
                                        Lambda_numerator_iterator;

    // types from the representation class
    typedef  typename Rep::Is_lp        Is_lp;

  protected:
    // protected types
    typedef  CGAL::Tag_true             Tag_true;
    typedef  CGAL::Tag_false            Tag_false;

  public:
    
    // creation
    Pricing_strategy_base( ) { }
    
    // destruction
    virtual ~Pricing_strategy_base( ) { }
    
    
    // initialization
    void  set( const Solver&           solver,
               CGAL::Verbose_ostream&  verbose_out)
        {
            solverP = &solver;
            voutP   = &verbose_out;
            set();
        }
    
    
    virtual  void  set( ) { }
    
    
    virtual  void  init( ) { }
    
    
    // access
    const Solver&        solver( ) const { return *solverP; }
    CGAL::Verbose_ostream&  vout  ( ) const { return *voutP; }
    
    
    // operations
    virtual  int   pricing( ) = 0;
    
    
    virtual  void  leaving_basis( int) { }
    
    virtual  void  transition( ) { }
    
    

  private:
    // data members
    const Solver*            solverP;   // the ambient QP solver
    CGAL::Verbose_ostream*   voutP;     // used for verbose output
};
  

CGAL_END_NAMESPACE
                  

#endif // CGAL_PRICING_STRATEGY_BASE_H

// ===== EOF ==================================================================
