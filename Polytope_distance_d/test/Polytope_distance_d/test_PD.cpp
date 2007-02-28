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
#include <CGAL/QP_solver/gmp_double.h> // temporarily, can be removed
                                       // once gmp_double is in CGAL
#include <CGAL/Cartesian_d.h>
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Optimisation_d_traits_d.h>
#include <CGAL/QP_solver/Double.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
 
typedef  CGAL::Cartesian_d< double >                                K;
typedef  CGAL::Optimisation_d_traits_d<K,CGAL::Double,double>  Traits;

#include "test_Polytope_distance_d.h"

// processing of file
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

  std::vector<K::Point_d>  p_points, q_points;
  
  // read points in P 
  for (int i = 0; i < n1; ++i) {
    std::vector<double>  coords( d);
    for (int j = 0; j < d; ++j)
      infile >> coords[ j];
    p_points.push_back( K::Point_d( d, coords.begin(), coords.end()));
  }

  // read points in Q
  for (int i = 0; i < n2; ++i) {
    std::vector<double>  coords( d);
    for (int j = 0; j < d; ++j)
      infile >> coords[ j];
    q_points.push_back( K::Point_d( d, coords.begin(), coords.end()));
  }
   
  std::cout << "Testing...\n";
  // call test function
  CGAL::test_Polytope_distance_d( p_points.begin(), p_points.end(),
                                  q_points.begin(), q_points.end(),
                                  Traits(), 1);
  std::cout << std::endl;
    
}

int main(const int ac,const char **av) {
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
  while (in >> filename) process(filename);

  return 0;
}
// ===== EOF ==================================================================
