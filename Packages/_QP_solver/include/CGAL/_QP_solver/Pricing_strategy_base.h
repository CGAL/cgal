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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>
                                                                               
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
