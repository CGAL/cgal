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
// file          : test/Min_annulus_d/test_Min_annulus_d.h
// package       : $CGAL_Package: Min_annulus_d $
// chapter       : Geometric Optimisation
//
// source        : web/Min_annulus_d.aw
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test function for smallest enclosing annulus
// ============================================================================

#ifndef CGAL_TEST_MIN_ANNULUS_D_H
#define CGAL_TEST_MIN_ANNULUS_D_H

// includes
#  include <CGAL/IO/Verbose_ostream.h>
#include <cassert>

namespace CGAL {

#define COVER(text,code) \
            verr0.out().width( 32); verr0 << text << "..." << flush; \
            verrX.out().width(  0); verrX << "==> " << text << endl \
              << "----------------------------------------" << endl; \
            { code } verr0 << "ok."; verr << endl;

template < class ForwardIterator, class Traits >
void
test_Min_annulus_d( ForwardIterator first, ForwardIterator last,
                    const Traits& traits, int verbose)
{
    using namespace std;

    typedef  CGAL::Min_annulus_d< Traits >  Min_annulus;
    typedef  typename Traits::Point_d      Point;

    CGAL::Verbose_ostream verr ( verbose >= 0);
    CGAL::Verbose_ostream verr0( verbose == 0);
    CGAL::Verbose_ostream verrX( verbose >  0);
    CGAL::set_pretty_mode( verr.out());

    bool  is_valid_verbose = ( verbose > 0);

    // constructors
    COVER( "default constructor",
        Min_annulus  ms( traits);
        assert( ms.is_valid( is_valid_verbose));
        assert( ms.is_empty());
    )

    COVER( "point set constructor",
        Min_annulus  ms( first, last, traits);
        assert( ms.is_valid( is_valid_verbose));
    )

    Min_annulus  min_annulus( first, last);
    COVER( "ambient dimension",
        Min_annulus  ms;
        assert( ms.ambient_dimension() == -1);
        verrX << min_annulus.ambient_dimension() << endl;
    )

    COVER( "(number of) points",
        verrX << min_annulus.number_of_points() << endl;
        typename Min_annulus::Point_iterator
            point_it = min_annulus.points_begin();
        for ( ; point_it != min_annulus.points_end(); ++point_it) {
            verrX << *point_it << endl;
        }
        assert( ( min_annulus.points_end() - min_annulus.points_begin())
                == min_annulus.number_of_points());
    )

    COVER( "(number of) support points",
        verrX << min_annulus.number_of_support_points() << endl;
        typename Min_annulus::Support_point_iterator
            point_it = min_annulus.support_points_begin();
        for ( ; point_it != min_annulus.support_points_end(); ++point_it) {
            verrX << *point_it << endl;
        }
        assert( ( min_annulus.support_points_end()
                  - min_annulus.support_points_begin())
                == min_annulus.number_of_support_points());
    )

    COVER( "(number of) inner support points",
        verrX << min_annulus.number_of_inner_support_points() << endl;
        typename Min_annulus::Inner_support_point_iterator
            point_it = min_annulus.inner_support_points_begin();
	typename Min_annulus::Inner_support_point_index_iterator
            point_index = min_annulus.inner_support_points_indices_begin();
        for ( ; point_it != min_annulus.inner_support_points_end();
              ++point_it, ++point_index) {
            verrX << *point_it << endl;
	    assert(*point_it == first[*point_index]);
        }
        assert( ( min_annulus.inner_support_points_end()
                  - min_annulus.inner_support_points_begin())
                == min_annulus.number_of_inner_support_points());
    )

    COVER( "(number of) outer support points",
        verrX << min_annulus.number_of_outer_support_points() << endl;
        typename Min_annulus::Outer_support_point_iterator
            point_it = min_annulus.outer_support_points_begin();
        typename Min_annulus::Outer_support_point_index_iterator
            point_index = min_annulus.outer_support_points_indices_begin();
        for ( ; point_it != min_annulus.outer_support_points_end();
              ++point_it, ++point_index) {
            verrX << *point_it << endl;
	    assert(*point_it == first[*point_index]);
        }
        assert( ( min_annulus.outer_support_points_end()
                  - min_annulus.outer_support_points_begin())
                == min_annulus.number_of_outer_support_points());
    )

    COVER( "center and squared radii",
        verrX << "center:";
        typename Min_annulus::Coordinate_iterator  coord_it;
        for ( coord_it  = min_annulus.center_coordinates_begin();
              coord_it != min_annulus.center_coordinates_end();
              ++coord_it) {
            verrX << ' ' << *coord_it;
        }
        verrX << endl << "squared inner radius: "
              << min_annulus.squared_inner_radius_numerator()   << " / "
              << min_annulus.squared_radii_denominator() << endl;
        verrX << endl << "squared outer radius: "
              << min_annulus.squared_outer_radius_numerator()   << " / "
              << min_annulus.squared_radii_denominator() << endl;
    )

    COVER( "predicates",
        CGAL::Bounded_side  bounded_side;
        bool                has_on_bounded_side;
        bool                has_on_boundary;
        bool                has_on_unbounded_side;
        Point               p;
        typename Min_annulus::Point_iterator
            point_it = min_annulus.points_begin();
        for ( ; point_it != min_annulus.points_end(); ++point_it) {
            p = *point_it;
            bounded_side          = min_annulus.bounded_side( p);
            has_on_bounded_side   = min_annulus.has_on_bounded_side( p);
            has_on_boundary       = min_annulus.has_on_boundary( p);
            has_on_unbounded_side = min_annulus.has_on_unbounded_side( p);
            verrX.out().width( 2);
            verrX << bounded_side          << "  "
                  << has_on_bounded_side   << ' '
                  << has_on_boundary       << ' '
                  << has_on_unbounded_side << endl;
            assert( bounded_side != CGAL::ON_UNBOUNDED_SIDE);
            assert( has_on_bounded_side || has_on_boundary);
            assert( ! has_on_unbounded_side);
        }
    )

    COVER( "clear",
        min_annulus.clear();
        verrX << "min_annulus is" << ( min_annulus.is_empty() ? "" : " not")
              << " empty." << endl;
        assert( min_annulus.is_empty());
    )

    COVER( "insert (single point)",
        min_annulus.insert( *first);
        assert( min_annulus.is_valid( is_valid_verbose));
        assert( min_annulus.is_degenerate());
    )

    COVER( "insert (point set)",
        min_annulus.insert( first, last);
        assert( min_annulus.is_valid( is_valid_verbose));
    )

    COVER( "traits class access",
        min_annulus.traits();
    )

    COVER( "I/O",
	   verrX << min_annulus  << endl;
    )
    verr0 << endl;
}

} //namespace CGAL

#endif // CGAL_TEST_MIN_ANNULUS_D_H

// ===== EOF ==================================================================
