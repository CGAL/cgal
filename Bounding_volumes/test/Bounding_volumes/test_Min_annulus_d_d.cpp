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
// file          : test/Min_annulus_d/test_Min_annulus_d_d.cpp
// package       : $CGAL_Package: Min_annulus_d $
// chapter       : Geometric Optimisation
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for Min_annulus (dD traits class)
// ============================================================================

#include <CGAL/Min_annulus_d.h>
#include <CGAL/Min_sphere_annulus_d_traits_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous_d.h>
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpzf RT;
typedef CGAL::Gmpq FT;
#else
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
typedef CGAL::MP_Float RT;
typedef CGAL::Quotient<CGAL::MP_Float> FT;
#endif
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

// exact kernels
typedef  CGAL::Cartesian_d<FT>                     CK1;
typedef  CGAL::Min_sphere_annulus_d_traits_d<CK1>        CTraits1;
typedef  CGAL::Homogeneous_d<RT>                   HK1;
typedef  CGAL::Min_sphere_annulus_d_traits_d<HK1>        HTraits1;

// inexact kernels
typedef  CGAL::Cartesian_d<double>                      CK2;
typedef  CGAL::Min_sphere_annulus_d_traits_d<CK2, RT, double> CTraits2;
typedef  CGAL::Homogeneous_d<double>                    HK2;
typedef  CGAL::Min_sphere_annulus_d_traits_d<HK2, RT, double> HTraits2;

#include <CGAL/Random.h>
#include <vector>

#include "test_Min_annulus_d.h"

template <class K, class Traits>
void process () 
{
  // generate point set
  std::vector<typename K::Point_d>  points;
  points.reserve( 100);
  {
    int d = 10;
    std::vector<double>  coords( d+1);
    int  i, j;
    double hom = 2.0;
    for ( i = 0; i < 100; ++i) {
      for (j=0; j<d; ++j) 
	coords[ j] = CGAL::get_default_random()( 0x100000);
      coords[d] = hom;
      points.push_back
	(typename K::Point_d(d, coords.begin(), coords.end()));
    }
    
  // call test function
  CGAL::test_Min_annulus_d(points.begin(), points.end(), Traits(), 0);
  }
}
    
// main
// ----
int main()
{
  // the following takes forever under Quotient<MP_Float>
#ifdef CGAL_USE_GMP 
    process<CK1, CTraits1>();
#endif
    process<HK1, HTraits1>(); 
    process<CK2, CTraits2>();
    process<HK2, HTraits2>();
}
 
// ===== EOF ==================================================================
