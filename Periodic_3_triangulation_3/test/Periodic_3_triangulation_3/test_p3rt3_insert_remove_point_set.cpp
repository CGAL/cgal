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

#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/Gmpz.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <set>


typedef CGAL::Regular_triangulation_euclidean_traits_3<CGAL::Epeck> RT_Epeck;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<RT_Epeck> Traits_Epeck;
/* Explicit instantiation.
 * It allows us to test if the template class, instantiated with given template parameters, compiles.
 */
template class CGAL::Periodic_3_regular_triangulation_3<Traits_Epeck>;


typedef CGAL::Regular_triangulation_euclidean_traits_3<CGAL::Epick> RT_Epick;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<RT_Epick> Traits_Epick;
/* Explicit instantiation.
 * It allows us to test if the template class, instantiated with given template parameters, compiles.
 */
template class CGAL::Periodic_3_regular_triangulation_3<Traits_Epick>;


template <class Kernel>
class Tests
{
public:
  typedef Kernel K;
  typedef typename K::FT FT;
  typedef CGAL::Regular_triangulation_euclidean_traits_3<K> RT;
  typedef CGAL::Periodic_3_regular_triangulation_traits_3<RT> Traits;

  typedef CGAL::Periodic_3_regular_triangulation_3<Traits> P3RT3;
  typedef typename P3RT3::Vertex_handle Vertex_handle;
  typedef typename P3RT3::Cell_handle Cell_handle;
  typedef typename P3RT3::Facet Facet;
  typedef typename P3RT3::Cell_iterator Cell_iterator;
  typedef typename P3RT3::Vertex_iterator Vertex_iterator;
  typedef typename P3RT3::Segment Segment;
  typedef typename P3RT3::Triangle Triangle;
  typedef typename P3RT3::Tetrahedron Tetrahedron;
  typedef typename P3RT3::Periodic_tetrahedron Periodic_tetrahedron;
  typedef typename P3RT3::Locate_type Locate_type;
  typedef typename P3RT3::Cell Cell;
  typedef typename P3RT3::Offset Offset;

  typedef typename Traits::Weighted_point Weighted_point;
  typedef typename Traits::Bare_point Bare_point;
  typedef typename Traits::Iso_cuboid_3 Iso_cuboid;


  static void test_insert_rnd_then_remove_all (unsigned pt_count, unsigned seed)
  {
    std::cout << "--- test_insert_rnd (" << pt_count << ", " << seed << ')' << std::endl;

    CGAL::Random random(seed);
    typedef CGAL::Creator_uniform_3<double,Bare_point>  Creator;
    CGAL::Random_points_in_cube_3<Bare_point, Creator> in_cube(0.5, random);

    Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
    P3RT3 p3rt3(iso_cuboid);

    std::ofstream stream("out_p3rt3_test");
    assert(stream);

    std::vector<Weighted_point> insert_set;
    insert_set.reserve(pt_count);
    std::vector<Weighted_point> remove_set;
    remove_set.reserve(pt_count);

    std::cout << "-- insert" << std::endl;
    for (unsigned cnt = 1; cnt <= pt_count; ++cnt)
    {
      Weighted_point p(*in_cube++, random.get_double(0., 0.015625));
      //    std::cout << cnt << " : " << p << std::endl;
      assert(p.weight() < 0.015625);
      stream << p << std::endl;

      std::ptrdiff_t hidden_point_count = 0;
      for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
        hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());

      Vertex_handle vh = p3rt3.insert(p);

      std::ptrdiff_t hidden_point_count_2 = 0;
      for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
        hidden_point_count_2 += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
      assert(hidden_point_count <= hidden_point_count_2);
      assert(hidden_point_count_2 + p3rt3.number_of_vertices() == cnt);

      std::cout << cnt << " - p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;
      if (vh == Vertex_handle())
      {
        if (find(insert_set.begin(), insert_set.end(), p) == insert_set.end())
          insert_set.push_back(p);
      }
      else
        insert_set.push_back(p);
    }

    stream.close();

    assert(p3rt3.is_valid());

    std::cout << "-- remove" << std::endl;
    unsigned cnt = 1;
    for (; p3rt3.number_of_vertices() != 0; ++cnt)
    {
      Vertex_iterator iter = p3rt3.vertices_begin();
      for (int j = random.get_int(0, static_cast<int>(p3rt3.number_of_vertices())); j; --j)
        ++iter;

      std::cout << cnt << " : " << iter->point() << std::endl;

      remove_set.push_back(iter->point());
      p3rt3.remove(iter);

      std::cout << "    p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;

      std::ptrdiff_t hidden_point_count = 0;
      for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
        hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());

      assert(hidden_point_count + cnt + p3rt3.number_of_vertices() == insert_set.size());
    }

    std::sort(insert_set.begin(), insert_set.end());
    std::sort(remove_set.begin(), remove_set.end());
    assert(insert_set == remove_set);

    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_stored_vertices() == 0);
  }

  static void test ()
  {
    //////    Iso_cuboid unitaire ->  0 <= weight < 0.015625
    test_insert_rnd_then_remove_all(800, 7);
    test_insert_rnd_then_remove_all(800, 12);
  }
};

int main (int argc, char** argv)
{
  std::cout << "TESTING ..." << std::endl;

  CGAL::force_ieee_double_precision();

  std::cout << "Epeck ..." << std::endl;
  Tests<CGAL::Epeck>::test();
  std::cout << "Epick ..." << std::endl;
  Tests<CGAL::Epick>::test();

  std::cout << "EXIT SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
