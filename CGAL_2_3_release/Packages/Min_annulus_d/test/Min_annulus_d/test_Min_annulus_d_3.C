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
// file          : test/Min_annulus_d/test_Min_annulus_d_3.C
// package       : $CGAL_Package: Min_annulus_d $
// chapter       : Geometric Optimisation
//
// source        : web/Min_annulus_d.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for smallest enclosing annulus (3D traits class)
// ============================================================================

// includes and typedefs
// ---------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Min_annulus_d.h>
#include <CGAL/Optimisation_d_traits_3.h>

// test variant 1 (needs LEDA)
#ifdef CGAL_USE_LEDA
# include <CGAL/leda_integer.h>
  typedef  CGAL::Cartesian<leda_integer>       K_1;
  typedef  CGAL::Optimisation_d_traits_3<K_1>  Traits_1;
# define TEST_VARIANT_1 \
    "Optimisation_d_traits_3< Cartesian<leda_integer> >"
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( leda_integer)
#endif

// test variant 2 (needs GMP)
#ifdef CGAL_USE_GMP
# include <CGAL/_QP_solver/Double.h>
  typedef  CGAL::Cartesian< int >                                 K_2;
  typedef  CGAL::Optimisation_d_traits_3<K_2,GMP::Double,double>  Traits_2;
# define TEST_VARIANT_2 \
    "Optimisation_d_traits_3< Cartesian<int>, GMP::Double, double >"
#endif

#include <CGAL/Random.h>
#include <vector>

#include "test_Min_annulus_d.h"

#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <iterator>

// main
// ----
int
main( int argc, char* argv[])
{
    CGAL_USING_NAMESPACE_STD
    
    int verbose = -1;
    if ( argc > 1) verbose = atoi( argv[ 1]);
    CGAL::Verbose_ostream  verr ( verbose >= 0); verr  << "";
    
    // test variant 1
    // --------------
    #ifdef TEST_VARIANT_1
    
        verr << endl
             << "==================================="
             << "===================================" << endl
             << "Testing `Min_annulus_d' with traits class model" << endl
             << "==> " << TEST_VARIANT_1 << endl
             << "==================================="
             << "===================================" << endl
             << endl;
    
        // generate point set
        std::vector<K_1::Point_3>  points_1;
        points_1.reserve( 100);
        CGAL::copy_n( CGAL::Random_points_on_sphere_3<K_1::Point_3>( 0x100000),
                      100, std::back_inserter( points_1));
    
        // call test function
        CGAL::test_Min_annulus_d( points_1.begin(), points_1.end(),
                                  Traits_1(), verbose);
    
    #endif
    
    // test variant 2
    // --------------
    #ifdef TEST_VARIANT_2
    
        verr << endl
             << "==================================="
             << "===================================" << endl
             << "Testing `Min_annulus_d' with traits class model" << endl
             << "==> " << TEST_VARIANT_2 << endl
             << "==================================="
             << "===================================" << endl
             << endl;
    
        // generate point set
        std::vector<K_2::Point_3>  points_2;
        points_2.reserve( 100);
        CGAL::copy_n( CGAL::Random_points_on_sphere_3<K_2::Point_3>( 0x100000),
                      100, std::back_inserter( points_2));
    
        // call test function
        CGAL::test_Min_annulus_d( points_2.begin(), points_2.end(),
                                  Traits_2(), verbose);
    
    #endif
    
    return 0;
}

// ===== EOF ==================================================================
