// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb, Mael Rouxel-Labbé
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>

#include <CGAL/Timer.h>

template <typename K>
struct Tester
{
  typedef CGAL::Robust_weighted_circumcenter_filtered_traits_3<K>   Gt;
  typedef typename Gt::FT                                           FT;
  typedef typename Gt::Point_3                                      Point_3;
  typedef typename Gt::Weighted_point_3                             Weighted_point_3;

  void operator()() const
  {
    //-------------------------------------------------------
    // Data generation : get 4 nearly coplanar points
    //-------------------------------------------------------
    FT little(1e-10);
    FT tiny(1e-25);

    Weighted_point_3 p1(little, 1, tiny);
    Weighted_point_3 p2(1, little, 0);
    Weighted_point_3 p3(-1*little, 1, 0);
    Weighted_point_3 p4(1, -1*little, 0);
    Weighted_point_3 p5(0, 0, 1);
    Weighted_point_3 p6(0, 1, 0);

    std::cerr << "Using points: p1[" << p1 << "]\tp2[" << p2
              << "]\tp3[" << p3 << "]\tp4[" << p4 << "]\tp5[" << p5
              << "]\tp6[" << p6 << "]\n";

    //-------------------------------------------------------
    // Test correctness
    //-------------------------------------------------------
    typename Gt::Construct_weighted_circumcenter_3 circumcenter =
        Gt().construct_weighted_circumcenter_3_object();

    Point_3 center = circumcenter(p1,p2);
    std::cerr << "\tcircumcenter(p1,p2)=[" << center << "]\n";

    center = circumcenter(p1,p3,p6);
    std::cerr << "\tcircumcenter(p1,p3,p6)=[" << center << "]\n";

    center = circumcenter(p1,p2,p5);
    std::cerr << "\tcircumcenter(p1,p3,p5)=[" << center << "]\n";

    // Use positive orientation
    center = circumcenter(p2,p3,p4,p1);
    std::cerr << "\tcircumcenter(p2,p3,p4,p1)=[" << center << "]\n";

    center = circumcenter(p1,p3,p2,p5);
    std::cerr << "\tcircumcenter(p1,p3,p2,p5)=[" << center << "]\n";

    //-------------------------------------------------------
    // Test speed
    //-------------------------------------------------------
    std::cerr << "Test speed: compute loops of: 999*c(p1,p3,p2,p5) "
              << "and 1*c(p2,p3,p4,p1)\n";

    CGAL::Timer timer;
    timer.start();
    int nb_loop = 0;
    while ( timer.time() < 0.5 )
    {
      // Compute 1000 fast queries
      for ( int i = 0 ; i < 999 ; ++i)
        circumcenter(p1,p3,p2,p5);

      // Compute 1 exact query
      circumcenter(p2,p3,p4,p1);
      ++nb_loop;
    }
    timer.stop();
    std::cerr << "\t" << nb_loop*1000/timer.time()
              << " circumcenter computation / second\n";


    std::cerr << "Test speed: compute loops of: 999*c(p2,p3,p4,p1) "
              << "and 1*c(p1,p3,p2,p5)\n";

    timer.reset();
    timer.start();
    nb_loop = 0;
    while ( timer.time() < 0.5 )
    {
      // Compute 1 exact queries
      for ( int i = 0 ; i < 999 ; ++i)
        circumcenter(p2,p3,p4,p1);

      // Compute 1 fast query
      circumcenter(p1,p3,p2,p5);
      ++nb_loop;
    }
    timer.stop();
    std::cerr << "\t" << nb_loop*1000/timer.time()
              << " circumcenter computation / second\n";
  }
};

int main()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
  typedef CGAL::Exact_predicates_exact_constructions_kernel   Epeck;

  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n";
  Tester<Epick> test_epic;
  test_epic();

  std::cerr << "\nTESTING WITH Exact_predicates_exact_constructions_kernel...\n";
  Tester<Epeck> test_epec;
  test_epec();

//  std::cerr << "\nTESTING WITH Filtered_kernel<Simple_cartesian<float> > kernel...\n";
//  Tester<Filtered_kernel<CGAL::Simple_cartesian<float> > > test_scf;
//  test_scf();

  return EXIT_SUCCESS;
}


