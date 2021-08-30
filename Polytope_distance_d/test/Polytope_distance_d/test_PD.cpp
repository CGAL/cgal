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
// file          : test/Polytope_distance_d/test_PD.C
// package       : $CGAL_Package: Polytope_distance_d $
// chapter       : Geometric Optimisation
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for polytope distance (dD traits class
// ============================================================================

// includes and typedefs
// ---------------------
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Optimisation_d_traits_d.h>
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

// fast and exact
typedef  CGAL::Cartesian_d< double >                    CK1;
typedef  CGAL::Optimisation_d_traits_d<CK1, RT, double> CTraits1;
typedef  CGAL::Homogeneous_d< double >                  HK1;
typedef  CGAL::Optimisation_d_traits_d<HK1, RT, double> HTraits1;

// slow and exact
typedef  CGAL::Cartesian_d< FT >                        CK2;
typedef  CGAL::Optimisation_d_traits_d<CK2>             CTraits2;
typedef  CGAL::Homogeneous_d< RT >                      HK2;
typedef  CGAL::Optimisation_d_traits_d<HK2>             HTraits2;

// overloaded function to create point from infile

CK1::Point_d create_point (int d, std::ifstream& infile, const CK1&)
{
  std::vector<double>  coords( d);
  for (int j = 0; j < d; ++j)
      infile >> coords[ j];
  return
    CK1::Point_d (d, coords.begin(), coords.end());
}

HK1::Point_d create_point (int d, std::ifstream& infile, const HK1&)
{
  // use some nontrivial homogeneous coordinates
  static double hom = 1.0;
  std::vector<double>  coords( d+1);
  for (int j = 0; j < d; ++j) {
    infile >> coords[ j];
    coords [ j] *= hom;
  }
  coords[d] = hom++;
  return
    HK1::Point_d (d, coords.begin(), coords.end());
}

CK2::Point_d create_point (int d, std::ifstream& infile, const CK2&)
{
  std::vector<FT>  coords( d);
  for (int j = 0; j < d; ++j)
      infile >> coords[ j];
  return
    CK2::Point_d (d, coords.begin(), coords.end());
}

HK2::Point_d create_point (int d, std::ifstream& infile, const HK2&)
{
  // use some nontrivial homogeneous coordinates
  static RT hom (1);
  std::vector<RT>  coords( d+1);
  for (int j = 0; j < d; ++j) {
    infile >> coords[ j];
    coords [ j] *= hom;
  }
  coords[d] = hom;
  hom += RT(1);
  return
    HK2::Point_d (d, coords.begin(), coords.end());
}

#include "test_Polytope_distance_d.h"

// processing of file
template <class K, class Traits>
void process (const std::string& filename)
{

  std::cout << "Reading file: " << filename << "\n";

  std::ifstream infile(filename.c_str());

  // read dimension
  int d;
  infile >> d;

  // read n1, n2
  int n1, n2;
  infile >> n1 >> n2;

  std::vector<typename K::Point_d>  p_points, q_points;

  // read points in P
  for (int i = 0; i < n1; ++i) {
    p_points.push_back(create_point (d, infile, K()));
  }

  // read points in Q
  for (int i = 0; i < n2; ++i) {
    q_points.push_back(create_point (d, infile, K()));
  }

  std::cout << "Testing...\n";
  // call test function
  CGAL::test_Polytope_distance_d( p_points.begin(), p_points.end(),
                                  q_points.begin(), q_points.end(),
                                  Traits(), 1);
  std::cout << std::endl;

}

int main(/* const int ac,const char **av */) {
  // we assume that av contains a list of files in the following
  // format:
  //     d n1 n2 p_1 p_2 ... p_n1 q_1 q_2 ... q_n2
  // where d is the dimension, n1 is the number of points in P,
  // n2 the number of points in Q, and p_1 p_2 ... p_n1, q_1 q_2 ... q_n2
  // are the actual points, stored coordinate-wise

  // read from standard input
  std::istream in(std::cin.rdbuf());

  // process input file(s):
  std::string filename;
  while (in >> filename) {
    process<CK1, CTraits1>(filename);
    process<HK1, HTraits1>(filename);
    // the following takes "forever" due to Quotient<MP_Float>
#ifdef CGAL_USE_GMP
    process<CK2, CTraits2>(filename);
#endif
    process<HK2, HTraits2>(filename);
  }

  return 0;
}
// ===== EOF ==================================================================
