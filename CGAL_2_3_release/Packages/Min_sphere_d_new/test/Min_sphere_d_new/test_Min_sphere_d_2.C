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
// file          : test/Min_sphere_d_new/test_Min_sphere_d_2.C
// package       : $CGAL_Package: Min_sphere_d_new $
// chapter       : Geometric Optimisation
//
// source        : web/Min_sphere_d.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Bernd Gärtner, Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for smallest enclosing sphere (2D traits class)
// ============================================================================

// includes
// --------
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Min_sphere_d_new.h>
#include <CGAL/Optimisation_d_traits_2.h>

#include <CGAL/Random.h>
#include <vector>

#include "test_Min_sphere_d.h"

#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <iterator>

// typedefs
// --------
// test variant 1 (needs LEDA)
#ifdef CGAL_USE_LEDA
# include <CGAL/leda_integer.h>
  typedef  CGAL::Cartesian<leda_integer>       K_1;
  typedef  CGAL::Optimisation_d_traits_2<K_1>  Traits_1;
# define TEST_VARIANT_1 \
    "Optimisation_d_traits_2< Cartesian<leda_integer> >"
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( leda_integer)
#endif

// test variant 2 (needs GMP)
#ifdef CGAL_USE_GMP
# include <CGAL/_QP_solver/Double.h>
  typedef  CGAL::Cartesian< int >                                 K_2;
  typedef  CGAL::Optimisation_d_traits_2<K_2,GMP::Double,double>  Traits_2;
# define TEST_VARIANT_2 \
    "Optimisation_d_traits_2< Cartesian<int>, GMP::Double, double >"
#endif


// comparing (needs LEDA)
#ifdef CGAL_USE_LEDA
  typedef  CGAL::Homogeneous<leda_integer>   K_3;
  typedef  CGAL::Min_circle_2_traits_2<K_3>  Traits_3;
  typedef  CGAL::Min_sphere_d<Traits_1>      Min_sphere_d;
  typedef  CGAL::Min_circle_2<Traits_3>      O_Min_sphere_d;
#endif

// main
// ----
int
main( int argc, char* argv[])
{
    CGAL_USING_NAMESPACE_STD

    // command line arguments
    int verbose = -1;
    if ( argc > 1) verbose = atoi( argv[ 1]);
    CGAL::Verbose_ostream  verr ( verbose >= 0); verr  << "";
    CGAL::Verbose_ostream  verr0( verbose == 0); verr0 << "";
    CGAL::Verbose_ostream  verrX( verbose >  0); verrX << "";

    // code coverage
    // -------------
    #ifdef TEST_VARIANT_1
        verr << endl
             << "==================================="
             << "===================================" << endl
             << "Testing `Min_sphere_d' with traits class model" << endl
             << "==> " << TEST_VARIANT_1 << endl
             << "==================================="
             << "===================================" << endl
             << endl;
    
        // generate point set
        std::vector<K_1::Point_2>  points_1;
        points_1.reserve( 100);
        CGAL::copy_n( CGAL::Random_points_on_circle_2<K_1::Point_2>( 0x100000),
                      100, std::back_inserter( points_1));
    
        // call test function
        CGAL::test_Min_sphere_d( points_1.begin(), points_1.end(),
                                 Traits_1(), verbose);
    #endif

    #ifdef TEST_VARIANT_2
        verr << endl
             << "==================================="
             << "===================================" << endl
             << "Testing `Min_sphere_d' with traits class model" << endl
             << "==> " << TEST_VARIANT_2 << endl
             << "==================================="
             << "===================================" << endl
             << endl;
    
        // generate point set
        std::vector<K_2::Point_2>  points_2;
        points_2.reserve( 100);
        CGAL::copy_n( CGAL::Random_points_on_circle_2<K_2::Point_2>( 0x100000),
                      100, std::back_inserter( points_2));
    
        // call test function
        CGAL::test_Min_sphere_d( points_2.begin(), points_2.end(),
                                 Traits_2(), verbose);
    #endif

    // additional tests
    // ----------------
    #ifdef CGAL_USE_LEDA
    
        verr << endl
             << "==================================="
             << "===================================" << endl
             << "Comparing `Min_sphere_d' with `Min_circle_2'" << endl
             << "==================================="
             << "===================================" << endl
             << endl;
    
        // convert point set
        std::vector<K_3::Point_2>  points_3;
        points_3.reserve( points_1.size());
        {
            unsigned int i;
            for ( i = 0; i < points_1.size(); ++i) {
                points_3.push_back( K_3::Point_2( points_1[ i][ 0],
                                                  points_1[ i][ 1]));
            }
        }
    
        // compute smallest enclosing spheres
        Min_sphere_d  ms( points_1.begin(), points_1.end(),
                          Traits_1(), verbose);
        verrX << endl << ms << endl;
        assert( ms.is_valid( verbose > 0));
        
        O_Min_sphere_d  o_ms( points_3.begin(), points_3.end(), false);
        verrX << endl << o_ms << endl;
        assert( o_ms.is_valid( verbose > 0));
        verrX << endl;
    
        // check center and squared radius
        COVER( "center",
            O_Min_sphere_d::Point  o_ms_center = o_ms.circle().center();
        
            verrX << "center (as point): " << ms.center()
                  << "  [NOTE: coordinates are truncated!]" << endl;
        
            int           d     = points_1[ 0].dimension();
            leda_integer  den   = ms.center_coordinates_begin()[ d];
            leda_integer  o_den = o_ms_center.homogeneous( d);
            for ( int j = 0; j < d; ++j) {
                assert( ms.center_coordinates_begin()[ j]*o_den
                        == o_ms_center.homogeneous( j)*den);
            }
            verrX << "centers are equal." << endl;
        );
        
        COVER( "squared radius",
            verrX << "squared radius: " << ms.squared_radius()
                  << "  [NOTE: value is truncated!]" << endl;
        
            assert( CGAL::Quotient<leda_integer>(
                        ms.squared_radius_numerator(),
                        ms.squared_radius_denominator())
                    == o_ms.circle().squared_radius());
            verrX << "squared radii are equal." << endl;
        );
    
    #endif

    return 0;
}

// ===== EOF ==================================================================
