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
// file          : test/Polytope_distance_d/test_Polytope_distance_d_2.C
// package       : $CGAL_Package: Polytope_distance_d $
// chapter       : Geometric Optimisation
//
// source        : web/Polytope_distance_d.aw
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for polytope distance (2D traits class)
// ============================================================================

#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf RT;
#else
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
typedef CGAL::MP_Float RT;
#endif
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

// fast and exact
typedef  CGAL::Cartesian< double >                      CK1;
typedef  CGAL::Polytope_distance_d_traits_2<CK1, RT, double> CTraits1;
typedef  CGAL::Homogeneous< double >                    HK1;
typedef  CGAL::Polytope_distance_d_traits_2<HK1, RT, double> HTraits1;

#include <CGAL/Random.h>
#include <vector>

#include "test_Polytope_distance_d.h"

template <class K, class Traits>
void process ()
{
  // generate point set
  std::vector<typename K::Point_2>  p_points, q_points;
  p_points.reserve( 50);
  q_points.reserve( 50);
  {
    int  i;
    double hom = 2.0;
    for ( i = 0; i < 50; ++i) {
      p_points.push_back
	(typename K::Point_2(
			     CGAL::get_default_random()( 0x100000),
			     CGAL::get_default_random()( 0x100000),
			     hom));
    }
    hom = 3.0;
    for ( i = 0; i < 50; ++i) {
      q_points.push_back
	(typename K::Point_2(
			     -CGAL::get_default_random()( 0x100000),
			     -CGAL::get_default_random()( 0x100000),
			     hom));
    }
  }

  // call test function
  CGAL::test_Polytope_distance_d( p_points.begin(), p_points.end(),
				  q_points.begin(), q_points.end(),
				  Traits(), 1);
}

// main
// ----
int
main()
{
    process<CK1, CTraits1>();
    process<HK1, HTraits1>();
}

// ===== EOF ==================================================================
