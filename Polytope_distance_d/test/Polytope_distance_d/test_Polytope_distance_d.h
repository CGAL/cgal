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
// file          : test/Polytope_distance_d/test_Polytope_distance_d.h
// package       : $CGAL_Package: Polytope_distance_d $
// chapter       : Geometric Optimization
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test function for polytope distance
// ============================================================================

#ifndef CGAL_TEST_POLYTOPE_DISTANCE_D_H
#define CGAL_TEST_POLYTOPE_DISTANCE_D_H

// includes
#  include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/use.h>
#include <cassert>

namespace CGAL {

#define COVER(text,code) \
            verr0.out().width( 26); verr0 << text << "..." << flush; \
            verrX.out().width(  0); verrX << "==> " << text << endl \
              << "----------------------------------------" << endl; \
            { code } verr0 << "ok."; verr0 << endl;

template < class ForwardIterator, class Traits >
void
test_Polytope_distance_d( ForwardIterator p_first, ForwardIterator p_last,
                          ForwardIterator q_first, ForwardIterator q_last,
                          const Traits& traits, int verbose)
{
    using namespace std;

    typedef  CGAL::Polytope_distance_d< Traits >  Poly_dist;
    CGAL_USE_TYPE(typename Traits::Point_d);

    CGAL::Verbose_ostream verr ( verbose >= 0);
    CGAL::Verbose_ostream verr0( verbose == 0);
    CGAL::Verbose_ostream verr1( verbose == 1);
    CGAL::Verbose_ostream verrX( verbose == 2);
    CGAL::IO::set_pretty_mode( verr.out());

    bool  is_valid_verbose = ( verbose > 1);

    // constructors
    COVER( "default constructor",
        Poly_dist  pd( traits);
        assert( pd.is_valid( is_valid_verbose));
        assert( ! pd.is_finite());
    )

    COVER( "point set constructor",
        Poly_dist  pd( p_first, p_last, q_first, q_last,
                       traits);
        assert( pd.is_valid( is_valid_verbose));
    )

    Poly_dist  poly_dist( p_first, p_last, q_first, q_last);
    COVER( "ambient dimension",
        Poly_dist  pd;
        assert( pd.ambient_dimension() == -1);
        verrX << poly_dist.ambient_dimension() << endl;
    )

    COVER( "(number of) points",
        verrX << poly_dist.number_of_points() << endl;
        verrX << endl << poly_dist.number_of_points_p() << endl;
        typename Poly_dist::Point_iterator
            point_p_it = poly_dist.points_p_begin();
        for ( ; point_p_it != poly_dist.points_p_end(); ++point_p_it) {
            verrX << *point_p_it << endl;
        }
        verrX << endl << poly_dist.number_of_points_q() << endl;
        typename Poly_dist::Point_iterator
            point_q_it = poly_dist.points_q_begin();
        for ( ; point_q_it != poly_dist.points_q_end(); ++point_q_it) {
            verrX << *point_q_it << endl;
        }
        assert( ( poly_dist.number_of_points_p()
                  + poly_dist.number_of_points_q())
                == poly_dist.number_of_points());
        assert( ( poly_dist.points_p_end() - poly_dist.points_p_begin())
                == poly_dist.number_of_points_p());
        assert( ( poly_dist.points_q_end() - poly_dist.points_q_begin())
                == poly_dist.number_of_points_q());
    )

    COVER( "(number of) support points",
        verrX << poly_dist.number_of_support_points() << endl;
        verrX << endl << poly_dist.number_of_support_points_p() << endl;
        typename Poly_dist::Support_point_index_iterator
            point_p_index = poly_dist.support_points_p_indices_begin();
        typename Poly_dist::Support_point_iterator
            point_p_it = poly_dist.support_points_p_begin();
        for ( ; point_p_it != poly_dist.support_points_p_end();
              ++point_p_it, ++point_p_index) {
            verrX << *point_p_it << endl;
            assert(*point_p_it == p_first [*point_p_index]);
        }
        verrX << endl << poly_dist.number_of_support_points_q() << endl;
        typename Poly_dist::Support_point_index_iterator
            point_q_index = poly_dist.support_points_q_indices_begin();
        typename Poly_dist::Support_point_iterator
            point_q_it = poly_dist.support_points_q_begin();
        for ( ; point_q_it != poly_dist.support_points_q_end();
              ++point_q_it,++point_q_index ) {
            verrX << *point_q_it << endl;
            assert(*point_q_it == q_first [*point_q_index]);
        }
        assert( ( poly_dist.number_of_support_points_p()
                  + poly_dist.number_of_support_points_q())
                == poly_dist.number_of_support_points());
        assert( ( poly_dist.support_points_p_end()
                  - poly_dist.support_points_p_begin())
                == poly_dist.number_of_support_points_p());
        assert( ( poly_dist.support_points_q_end()
                  - poly_dist.support_points_q_begin())
                == poly_dist.number_of_support_points_q());

    )

    COVER( "realizing points",
        typename Poly_dist::Coordinate_iterator  coord_it;
        verrX << "p:";
        for ( coord_it  = poly_dist.realizing_point_p_coordinates_begin();
              coord_it != poly_dist.realizing_point_p_coordinates_end();
              ++coord_it) {
            verrX << ' ' << *coord_it;
        }
        verrX << endl;
        verrX << "q:";
        for ( coord_it  = poly_dist.realizing_point_q_coordinates_begin();
              coord_it != poly_dist.realizing_point_q_coordinates_end();
              ++coord_it) {
            verrX << ' ' << *coord_it;
        }
        verrX << endl;
    )

    COVER( "squared distance",
        verrX << poly_dist.squared_distance_numerator()   << " / "
              << poly_dist.squared_distance_denominator() << endl;
    )

    COVER( "clear",
        poly_dist.clear();
        verrX << "poly_dist is" << ( poly_dist.is_finite() ? "" : " not")
              << " finite." << endl;
        assert( ! poly_dist.is_finite());
    )

    COVER( "insert (single point)",
        poly_dist.insert_p( *p_first);
        poly_dist.insert_q( *q_first);
        assert( poly_dist.is_valid( is_valid_verbose));
    )

    COVER( "insert (point set)",
        poly_dist.insert( q_first, q_last, p_first, p_last);
        assert( poly_dist.is_valid( is_valid_verbose));
        poly_dist.clear();
        poly_dist.insert_p( p_first, p_last);
        poly_dist.insert_q( q_first, q_last);
        assert( poly_dist.is_valid( is_valid_verbose));
    )

    COVER( "traits class access",
        poly_dist.traits();
    )

    COVER( "I/O",
        verrX << poly_dist;
    )
}

} //namespace CGAL

#endif // CGAL_TEST_POLYTOPE_DISTANCE_D_H

// ===== EOF ==================================================================
