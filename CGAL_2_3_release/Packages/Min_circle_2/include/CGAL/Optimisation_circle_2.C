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
// file          : include/CGAL/Optimisation_circle_2.C
// package       : $CGAL_Package: Min_circle_2 $
// chapter       : Geometric Optimisation
//
// source        : web/Min_circle_2.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>, Bernd Gärtner
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: 2D Optimisation Circle
// ============================================================================

// includes
#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#  include <CGAL/Optimisation/assertions.h>
#endif

CGAL_BEGIN_NAMESPACE

// Class implementation (continued)
// ================================

// I/O
// ---
template < class K_ >
std::ostream&
operator << ( std::ostream& os, const CGAL::Optimisation_circle_2<K_>& c)
{
    switch ( CGAL::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << "CGAL::Optimisation_circle_2( "
           << c.center() << ", "
           << c.squared_radius() << ')';
        break;

      case CGAL::IO::ASCII:
        os << c.center() << ' ' << c.squared_radius();
        break;

      case CGAL::IO::BINARY:
        os << c.center();
        CGAL::write( os, c.squared_radius());
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::get_mode( os) invalid!");
        break; }

    return( os);
}

template < class K_ >
std::istream&
operator >> ( std::istream& is, CGAL::Optimisation_circle_2<K_>& c)
{
    typedef  CGAL::Optimisation_circle_2<K_>::Point     Point;
    typedef  CGAL::Optimisation_circle_2<K_>::Distance  Distance;

    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
        cerr << std::endl;
        cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;

      case CGAL::IO::ASCII: {
        Point     center;
        Distance  squared_radius;
        is >> center >> squared_radius;
        c.set( center, squared_radius); }
        break;

      case CGAL::IO::BINARY: {
        Point     center;
        Distance  squared_radius;
        is >> center;
        CGAL::read( is, squared_radius);
        c.set( center, squared_radius); }
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::get_mode( is) invalid!");
        break; }

    return( is);
}

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
