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
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Min_circle_2.C
// package       : Min_circle_2 3.10.5 (19 Mar 2001)
// chapter       : Geometric Optimisation
//
// source        : web/Min_circle_2.aw
// revision      : 5.30
// revision_date : 2001/03/19 18:06:10
//
// author(s)     : Sven Schönherr, Bernd Gärtner
// maintainer    : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: 2D Smallest Enclosing Circle
// ============================================================================

CGAL_BEGIN_NAMESPACE

// Class implementation (continued)
// ================================
// I/O
// ---
template < class _Traits >
std::ostream&
operator << ( std::ostream& os,
              const Min_circle_2<_Traits>& min_circle)
{
    CGAL_USING_NAMESPACE_STD

    typedef  Min_circle_2<_Traits>::Point  Point;
    typedef  ostream_iterator<Point>       Os_it;

    switch ( CGAL::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << endl;
        os << "CGAL::Min_circle_2( |P| = " << min_circle.number_of_points()
           << ", |S| = " << min_circle.number_of_support_points() << endl;
        os << "  P = {" << endl;
        os << "    ";
        copy( min_circle.points_begin(), min_circle.points_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  S = {" << endl;
        os << "    ";
        copy( min_circle.support_points_begin(),
              min_circle.support_points_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  circle = " << min_circle.circle() << endl;
        os << ")" << endl;
        break;

      case CGAL::IO::ASCII:
        copy( min_circle.points_begin(), min_circle.points_end(),
              Os_it( os, "\n"));
        break;

      case CGAL::IO::BINARY:
        copy( min_circle.points_begin(), min_circle.points_end(),
              Os_it( os));
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::get_mode( os) invalid!");
        break; }

    return( os);
}

template < class _Traits >
std::istream&
operator >> ( std::istream& is, CGAL::Min_circle_2<_Traits>& min_circle)
{
    CGAL_USING_NAMESPACE_STD

    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
        cerr << endl;
        cerr << "Stream must be in ascii or binary mode" << endl;
        break;

      case CGAL::IO::ASCII:
      case CGAL::IO::BINARY:
        typedef  CGAL::Min_circle_2<_Traits>::Point  Point;
        typedef  istream_iterator<Point>            Is_it;
        min_circle.clear();
        min_circle.insert( Is_it( is), Is_it());
        break;

      default:
        CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
        break; }

    return( is);
}

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
