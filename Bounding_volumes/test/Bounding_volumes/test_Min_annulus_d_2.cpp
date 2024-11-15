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
// file          : test/Min_annulus_d/test_Min_annulus_d_2.cpp
// package       : $CGAL_Package: Min_annulus_d $
// chapter       : Geometric Optimization
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for Min_annulus (2D traits class)
// ============================================================================

#include <CGAL/Min_annulus_d.h>
#include <CGAL/Min_sphere_annulus_d_traits_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
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
typedef  CGAL::Cartesian<FT>                     CK1;
typedef  CGAL::Min_sphere_annulus_d_traits_2<CK1>      CTraits1;
typedef  CGAL::Homogeneous<RT>                   HK1;
typedef  CGAL::Min_sphere_annulus_d_traits_2<HK1>      HTraits1;

// inexact kernels
typedef  CGAL::Cartesian<double>                        CK2;
typedef  CGAL::Min_sphere_annulus_d_traits_2<CK2, RT, double> CTraits2;
typedef  CGAL::Homogeneous<double>                      HK2;
typedef  CGAL::Min_sphere_annulus_d_traits_2<HK2, RT, double> HTraits2;

#include <CGAL/Random.h>
#include <vector>

#include "test_Min_annulus_d.h"

template <class K, class Traits>
void process ()
{
  // generate point set
  std::vector<typename K::Point_2>  points;
  points.reserve( 100);
  {
    double hom = 2.0;
    for ( int i = 0; i < 100; ++i) {
      points.push_back
        (typename K::Point_2
         (CGAL::get_default_random()( 0x100000),
          CGAL::get_default_random()( 0x100000),
          hom));
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
