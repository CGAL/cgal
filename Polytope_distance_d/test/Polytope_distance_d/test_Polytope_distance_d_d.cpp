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
// file          : test/Polytope_distance_d/test_Polytope_distance_d_3.C
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
// implementation: test program for polytope distance (dD traits class)
// ============================================================================

#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous_d.h>
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
typedef  CGAL::Cartesian_d< double >                      CK1;
typedef  CGAL::Polytope_distance_d_traits_d<CK1, RT, double>   CTraits1;
typedef  CGAL::Homogeneous_d< double >                    HK1;
typedef  CGAL::Polytope_distance_d_traits_d<HK1, RT, double>   HTraits1;

#include <CGAL/Random.h>
#include <vector>

#include "test_Polytope_distance_d.h"

template <class K, class Traits>
void process () 
{
  // generate point set
  std::vector<typename K::Point_d>  p_points, q_points;
  p_points.reserve( 50);
  q_points.reserve( 50);
  {
    int d = 10;
    std::vector<double>  coords( d+1);
    int  i, j;
    double hom = 2.0;
    for ( i = 0; i < 50; ++i) {
      for (j=0; j<d; ++j) 
	coords[ j] = CGAL::get_default_random()( 0x100000);
      coords[d] = hom;
      p_points.push_back
	(typename K::Point_d(d, coords.begin(), coords.end()));
    }
    hom = 3.0;
    for ( i = 0; i < 50; ++i) {
      for (j=0; j<d; ++j) 
	coords[ j] = -CGAL::get_default_random()( 0x100000);
      coords[d] = hom;
      q_points.push_back
	(typename K::Point_d(d, coords.begin(), coords.end()));
    }
    
  // call test function
  CGAL::test_Polytope_distance_d( p_points.begin(), p_points.end(),
				  q_points.begin(), q_points.end(),
				  Traits(), 2);
  }
}
    
// main
// ----
int main()
{
  std::cerr << "Testing with " << typeid(CTraits1()).name() << std::endl;
    process<CK1, CTraits1>();
    process<HK1, HTraits1>();
}
 
// ===== EOF ==================================================================
