// ============================================================================
//
// Copyright (c) 1997-2003 The CGAL Consortium
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
// file          : include/CGAL/QP_engine/QPE__partial_base.h
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2003/07
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Base Class for Partial Pricing of the QPE Solver
// ============================================================================

#ifndef CGAL_QPE__PARTIAL_BASE_H
#define CGAL_QPE__PARTIAL_BASE_H

// includes
#include <CGAL/QP_engine/QPE_pricing_strategy.h>
#include <CGAL/Random.h>
#include <algorithm>
#include <vector>
#include <cmath>

CGAL_BEGIN_NAMESPACE

// ==================
// class declarations
// ==================
template < class Rep_ >
class QPE__partial_base;

// ===============
// class interface
// ===============
template < class Rep_ >
class QPE__partial_base : virtual public QPE_pricing_strategy<Rep_> {

    // self
    typedef  Rep_                       Rep;
    typedef  QPE__partial_base<Rep>     Self;
    typedef  QPE_pricing_strategy<Rep>  Base;

  protected:

    // types
    typedef  typename Base::QP_solver::Indices         Indices;
    typedef  typename Base::QP_solver::Index_iterator  Index_iterator;
    typedef  typename Base::QP_solver::Index_const_iterator  Index_const_iterator;

    // construction
    QPE__partial_base( bool  randomize, Random&  random);

    // initialization
    virtual  void  init( );

    // access
    Index_const_iterator    active_set_begin( ) const { return N.begin();   }
    Index_const_iterator    active_set_end  ( ) const { return N.begin()+s; }

    Index_const_iterator  inactive_set_begin( ) const { return N.begin()+s; }
    Index_const_iterator  inactive_set_end  ( ) const { return N.end  ();   }

    // operations
    void  entering_basis( Index_const_iterator  it);
    void  activating    ( Index_const_iterator& it);

    virtual  void  leaving_basis( int i);
    virtual  void  transition( );

  private:

    // data members
    Indices                  N;         // non-basis;
    int                      s;         // size of active set

    bool                     permute;   // flag: permute initial non-basis
    Random&                  rand_src;  // random source
};

// ----------------------------------------------------------------------------

// ===============================
// class implementation (template)
// ===============================

// construction
template < class Rep_ >  inline
QPE__partial_base<Rep_>::
QPE__partial_base( bool  randomize, Random&  random)
    : permute( randomize), rand_src( random)
{ }

// initialization
template < class Rep_ >
void
QPE__partial_base<Rep_>::
init( )
{
    // initialize indices of non-basic variables
    int  w = solver().number_of_working_variables();
    int  b = solver().number_of_basic_variables();

    if ( ! N.empty()) N.clear();
    N.reserve( w-b);                                        // use 'w' ???
    for ( int i = 0; i < w; ++i) {
	if ( ! solver().is_basic( i)) N.push_back( i);
    }
    if ( permute) std::random_shuffle( N.begin(), N.end(), rand_src);

    // initialize size of active set
    int  n = solver().number_of_variables();
    int  m = solver().number_of_constraints();

    s = static_cast< int>( m*std::sqrt( n/2.0));
}

// operations
template < class Rep_ >  inline
void
QPE__partial_base<Rep_>::
entering_basis( Index_const_iterator it)
{
    CGAL_qpe_precondition( it >= active_set_begin() && it < active_set_end());

    // remove from active set
    --s;
    const_cast< typename Indices::value_type&>( *it) = N[ s];
    N[ s] = N.back();
    N.pop_back();
}

template < class Rep_ >  inline
void
QPE__partial_base<Rep_>::
activating( Index_const_iterator& it)
{
    CGAL_qpe_precondition(
	it >= inactive_set_begin() && it < inactive_set_end());

    // 'append' to active set
    std::swap( const_cast< typename Indices::value_type&>( *it), N[ s]);
    it = N.begin()+s;
    ++s;
}

template < class Rep_ >
void
QPE__partial_base<Rep_>::
leaving_basis( int i)
{
    // all non-basic variables active?
    if ( s == static_cast< int>( N.size())) {

	// append at the end of all non-basic variables
	N.push_back( i);

    } else {

	// insert at the end of the current active subset
	N.push_back( N[ s]);
	N[ s] = i;
    }
    ++s;
}


template < class Rep_ >
void
QPE__partial_base<Rep_>::
transition( )
{
    // remove artificial variables from non-basis
    int  w = solver().number_of_working_variables();
    N.erase( std::partition( N.begin(), N.end(),
			     std::bind2nd( std::less<int>(), w)),
	     N.end());

    // initialize size of active set
    int  n = solver().number_of_variables();
    int  m = solver().number_of_constraints();

    s = std::min( static_cast< unsigned int>( m*std::sqrt( n/2.0)), N.size());
}

CGAL_END_NAMESPACE

#endif // CGAL_QPE__PARTIAL_BASE_H

// ===== EOF ==================================================================
