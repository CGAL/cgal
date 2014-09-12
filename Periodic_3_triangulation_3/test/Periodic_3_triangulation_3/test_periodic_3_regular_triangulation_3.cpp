// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
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
// $Date$
// 
//
// Author(s)     : Aymeric PELLE <aymeric.pelle@sophia.inria.fr>

#include <CGAL/Periodic_3_Regular_triangulation_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/Gmpz.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Periodic_3_Regular_triangulation_3.h>

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random.hpp>

#include <cassert>
#include <iostream>
#include <fstream>


typedef CGAL::Epick K;
typedef CGAL::Epick::FT FT;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Regular_traits;
typedef CGAL::Periodic_3_Regular_triangulation_traits_3<Regular_traits> Traits;
typedef typename Traits::Weighted_point Weighted_point;
typedef typename Traits::Bare_point Bare_point;
typedef typename Traits::Iso_cuboid_3 Iso_cuboid;

template class CGAL::Periodic_3_Regular_triangulation_3<Traits>;

typedef CGAL::Periodic_3_Regular_triangulation_3<Traits> P3RT3;


Weighted_point read_wpoint (std::istream& os)
{
  FT x = 0., y = 0., z = 0., w = 0.;
  os >> x;
  os >> y;
  os >> z;
  os >> w;
  return Weighted_point(Bare_point(x, y, z), w);
}


void test_construction ()
{
  P3RT3 p3rt3;
  assert(p3rt3.is_valid(true));
}

void test_insert_1 ()
{
  P3RT3 p3rt3;

  Weighted_point p(Bare_point(0,0,0), 0.1);
  p3rt3.insert(p);

  assert(p3rt3.is_valid(true));
  assert(p3rt3.number_of_vertices() == 1);
  assert(p3rt3.number_of_stored_vertices() == 27);
}

void test_insert_rnd_100 ()
{
  P3RT3 p3rt3;

  srand(time(NULL));

//  boost::random::random_device rd;
//  boost::mt19937 gen(rd);
//  boost::random::uniform_real_distribution<> c_dis(0., 1.);
//  auto gen_coord = boost::bind(c_dis, gen);
//  boost::random::uniform_real_distribution<> w_dis(0., 0.1);
//  auto gen_weight = boost::bind(w_dis, gen);
  auto gen_coord = [](){ return static_cast<double>(rand() % 1000) / 1000.; };
  auto gen_weight = [](){ return static_cast<double>(rand() % 100) / 1000.; };

  std::ofstream stream("out");
  assert(stream);

  for (unsigned cnt = 100; cnt--; )
  {
    Weighted_point p(Bare_point(gen_coord(), gen_coord(), gen_coord()), gen_weight());
    std::cout << p << std::endl;
    stream << p << std::endl;
    p3rt3.insert(p);
    assert(p3rt3.is_valid());
  }

  assert(p3rt3.number_of_vertices() == 100);
  assert(p3rt3.number_of_stored_vertices() == 2700);
  assert(p3rt3.is_valid(true));
}

void test_insert_from_file (const char* filename)
{
  P3RT3 p3rt3;

  srand(time(NULL));

  std::ifstream stream(filename);
  assert(stream);

  while (stream && !(stream.eof()))
  {
    Weighted_point p = read_wpoint(stream);
    std::cout << p << std::endl;
    p3rt3.insert(p);
    assert(p3rt3.is_valid());
  }
  assert(p3rt3.is_valid(true));
}

int main (int argc, char** argv)
{
  std::cout << "TESTING ..." << std::endl;

  test_construction();
  test_insert_1();
  if (argc > 1 && strlen(argv[1]) > 0)
    test_insert_from_file(argv[1]);

  std::cout << "EXIT SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
