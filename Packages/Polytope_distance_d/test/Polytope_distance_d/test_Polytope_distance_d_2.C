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
// file          : test/Polytope_distance_d/test_Polytope_distance_d_2.C
// package       : Polytope_distance_d 1.0.2 (20 Mar 2001)
// chapter       : Geometric Optimisation
//
// source        : web/Polytope_distance_d.aw
// revision      : 1.8
// revision_date : 2000/10/04 14:04:28
//
// author(s)     : Sven Schönherr
// maintainer    : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for polytope distance (2D traits class)
// ============================================================================

// includes and typedefs
// ---------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Optimisation_d_traits_2.h>

// test variant 1 (needs LEDA)
#ifdef CGAL_USE_LEDA
# include <CGAL/leda_integer.h>
  typedef  CGAL::Cartesian<leda_integer>       R_1;
  typedef  CGAL::Optimisation_d_traits_2<R_1>  Traits_1;
# define TEST_VARIANT_1 \
    "Optimisation_d_traits_2< Cartesian<leda_integer> >"
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( leda_integer)
#endif

// test variant 2 (needs GMP)
#ifdef CGAL_USE_GMP
# include <CGAL/_QP_solver/Double.h>
  typedef  CGAL::Cartesian< int >                                 R_2;
  typedef  CGAL::Optimisation_d_traits_2<R_2,GMP::Double,double>  Traits_2;
# define TEST_VARIANT_2 \
    "Optimisation_d_traits_2< Cartesian<int>, GMP::Double, double >"
#endif

#include <CGAL/Random.h>
#include <vector>

#include "test_Polytope_distance_d.h"


// main
// ----
int
main( int argc, char* argv[])
{
    CGAL_USING_NAMESPACE_STD
    
    int verbose = -1;
    if ( argc > 1) verbose = atoi( argv[ 1]);
    CGAL::Verbose_ostream  verr( verbose >= 0);
    
    // test variant 1
    // --------------
    #ifdef TEST_VARIANT_1
    
        verr << endl
             << "Testing `Polytope_distance_d' with traits class model" << endl
             << "==> " << TEST_VARIANT_1 << endl
             << "==================================="
             << "===================================" << endl
             << endl;
    
        // generate point set
        std::vector<R_1::Point_2>  p_points_1, q_points_1;
        p_points_1.reserve( 50);
        q_points_1.reserve( 50);
        {
            int  i;
            for ( i = 0; i < 50; ++i) {
                p_points_1.push_back( R_1::Point_2(
                    CGAL::default_random( 0x100000),
                    CGAL::default_random( 0x100000)));
            }
            for ( i = 0; i < 50; ++i) {
                q_points_1.push_back( R_1::Point_2(
                    -CGAL::default_random( 0x100000),
                    -CGAL::default_random( 0x100000)));
            }
        }
    
        // call test function
        // call test function
        CGAL::test_Polytope_distance_d( p_points_1.begin(), p_points_1.end(),
                                        q_points_1.begin(), q_points_1.end(),
                                        Traits_1(), verbose);
    
    #endif
    
    // test variant 2
    // --------------
    #ifdef TEST_VARIANT_2
    
        verr << endl
             << "Testing `Polytope_distance_d' with traits class model" << endl
             << "==> " << TEST_VARIANT_2 << endl
             << "==================================="
             << "===================================" << endl
             << endl;
    
        // generate point set
        std::vector<R_2::Point_2>  p_points_2, q_points_2;
        p_points_2.reserve( 50);
        q_points_2.reserve( 50);
        {
            int  i;
            for ( i = 0; i < 50; ++i) {
                p_points_2.push_back( R_2::Point_2(
                    CGAL::default_random( 0x100000),
                    CGAL::default_random( 0x100000)));
            }
            for ( i = 0; i < 50; ++i) {
                q_points_2.push_back( R_2::Point_2(
                    -CGAL::default_random( 0x100000),
                    -CGAL::default_random( 0x100000)));
            }
        }
    
        // call test function
        // call test function
        CGAL::test_Polytope_distance_d( p_points_2.begin(), p_points_2.end(),
                                        q_points_2.begin(), q_points_2.end(),
                                        Traits_2(), verbose);
    
    #endif
    
    return 0;
}

// ===== EOF ==================================================================
