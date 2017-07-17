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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>

typedef CGAL::Periodic_3_regular_triangulation_traits_3<CGAL::Epeck> Traits_Epeck;
/* Explicit instantiation.
 * It allows us to test if the template class, instantiated with given template parameters, compiles.
 */
template class CGAL::Periodic_3_regular_triangulation_3<Traits_Epeck>;


typedef CGAL::Periodic_3_regular_triangulation_traits_3<CGAL::Epick> Traits_Epick;
/* Explicit instantiation.
 * It allows us to test if the template class, instantiated with given template parameters, compiles.
 */
template class CGAL::Periodic_3_regular_triangulation_3<Traits_Epick>;


template <class Kernel>
class Tests
{
public:
  typedef Kernel                                              K;
  typedef typename K::FT                                      FT;
  typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>  Traits;

  typedef CGAL::Periodic_3_regular_triangulation_3<Traits>    P3RT3;
  typedef typename P3RT3::Vertex_handle                       Vertex_handle;
  typedef typename P3RT3::Cell_handle                         Cell_handle;
  typedef typename P3RT3::Facet                               Facet;
  typedef typename P3RT3::Cell                                Cell;
  typedef typename P3RT3::Vertex_iterator                     Vertex_iterator;
  typedef typename P3RT3::Cell_iterator                       Cell_iterator;
  typedef typename P3RT3::Segment                             Segment;
  typedef typename P3RT3::Triangle                            Triangle;
  typedef typename P3RT3::Tetrahedron                         Tetrahedron;
  typedef typename P3RT3::Periodic_tetrahedron                Periodic_tetrahedron;
  typedef typename P3RT3::Locate_type                         Locate_type;
  typedef typename P3RT3::Offset                              Offset;

  typedef typename Traits::Weighted_point_3                   Weighted_point_3;
  typedef typename Traits::Point_3                            Point_3;
  typedef typename Traits::Iso_cuboid_3                       Iso_cuboid;


  static Weighted_point_3 read_wpoint (std::istream& stream)
  {
    FT x = 0., y = 0., z = 0., w = 0.;
    stream >> x;
    assert(stream && !stream.eof());
    stream >> y;
    assert(stream && !stream.eof());
    stream >> z;
    assert(stream && !stream.eof());
    stream >> w;
    assert(stream);
    return Weighted_point_3(Point_3(x, y, z), w);
  }

  static void test_construction ()
  {
    std::cout << "--- test_construction" << std::endl;

    P3RT3 p3rt3;
    assert(p3rt3.is_valid());
  }

  static void test_insert_1 ()
  {
    std::cout << "--- test_insert_1" << std::endl;

    P3RT3 p3rt3;

    Weighted_point_3 p(Point_3(0,0,0), 0.01);
    p3rt3.insert(p);

    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 1);
    assert(p3rt3.number_of_stored_vertices() == 27);

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 0);
  }

  static void test_insert_point ()
  {
    std::cout << "--- test_insert_point" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0,0,0, 1,1,1));

    Vertex_handle vh;
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.9,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.9,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.9),0.01));
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    Weighted_point_3 p(Point_3(0.4, 0.4, 0.4), 0.001);
    vh = p3rt3.insert(p);
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 5);
    assert(p3rt3.number_of_stored_vertices() == 135);

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 0);
  }

  static void test_insert_hidden_point ()
  {
    std::cout << "--- test_insert_hidden_point" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0,0,0, 1,1,1));

    Vertex_handle vh;
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.9,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.9,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.9),0.01));
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 0);

    Weighted_point_3 hidden_point(Point_3(0.101, 0.101, 0.101), 0.001);
    vh = p3rt3.insert(hidden_point);
    assert(vh == Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    unsigned hidden_found_count = 0;
    hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    {
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
      hidden_found_count += static_cast<unsigned>(std::find(iter->hidden_points_begin(), iter->hidden_points_end(), hidden_point) != iter->hidden_points_end());
    }
    assert(hidden_point_count == 1);
    assert(hidden_found_count == 1);
  }

  static void test_insert_hiding_point ()
  {
    std::cout << "--- test_insert_hiding_point" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0,0,0, 1,1,1));

    Vertex_handle vh;
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.9,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.9,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.9),0.01));
    assert(vh != Vertex_handle());
    Weighted_point_3 hidden_point(Point_3(0.101, 0.101, 0.101), 0.001);
    vh = p3rt3.insert(hidden_point);
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 0);

    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    unsigned hidden_found_count = 0;
    hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    {
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
      hidden_found_count += static_cast<unsigned>(std::find(iter->hidden_points_begin(), iter->hidden_points_end(), hidden_point) != iter->hidden_points_end());
    }
    assert(hidden_point_count == 1);
    assert(hidden_found_count == 1);
  }

  static void test_insert_a_point_twice ()
  {
    std::cout << "--- test_insert_a_point_twice" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0,0,0, 1,1,1));

    Vertex_handle vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 1);
    assert(p3rt3.number_of_stored_vertices() == 27);

    Vertex_handle vh2 = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1),0.01));
    assert(vh2 == vh);
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 1);
    assert(p3rt3.number_of_stored_vertices() == 27);

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 0);
  }

  static void test_insert_two_points_with_the_same_position ()
  {
    std::cout << "--- test_insert_two_points_with_the_same_position" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0,0,0, 1,1,1));

    Weighted_point_3 hidden_point(Point_3(0.1,0.1,0.1),0.01);
    Vertex_handle vh = p3rt3.insert(hidden_point);
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 1);
    assert(p3rt3.number_of_stored_vertices() == 27);

    Vertex_handle vh2 = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1),0.015));
    assert(vh2 != Vertex_handle());
    assert(vh2 != vh);
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 1);
    assert(p3rt3.number_of_stored_vertices() == 27);

    unsigned hidden_found_count = 0;
    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    {
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
      hidden_found_count += static_cast<unsigned>(std::find(iter->hidden_points_begin(), iter->hidden_points_end(), hidden_point) != iter->hidden_points_end());
    }
    assert(hidden_point_count == 1);
    assert(hidden_found_count == 1);
  }

  static void test_remove ()
  {
    std::cout << "--- test_remove" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0,0,0, 1,1,1));

    Vertex_handle vh;
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.9,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.9,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.9),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    Weighted_point_3 hidden_point(Point_3(0.101, 0.101, 0.101), 0.001);
    Vertex_handle vhh = p3rt3.insert(hidden_point);
    assert(vhh == Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    p3rt3.remove(vh);
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 0);

    unsigned point_found_count = 0;
    for (Vertex_iterator iter = p3rt3.vertices_begin(), end_iter = p3rt3.vertices_end(); iter != end_iter; ++iter)
      point_found_count += (iter->point() == hidden_point);
    assert(point_found_count == 27);
  }

  static void test_insert_rnd_as_delaunay (unsigned pt_count, double weight)
  {
    std::cout << "--- test_insert_rnd_as_delaunay (" << pt_count << ',' << weight << ')' << std::endl;

    CGAL::Random random(7);
    typedef CGAL::Creator_uniform_3<double, Point_3>  Creator;
    CGAL::Random_points_in_cube_3<Point_3, Creator> in_cube(0.5, random);

    Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.6, 0.6, 0.6);
    P3RT3 p3rt3(iso_cuboid);

    for (unsigned cnt = 1; cnt <= pt_count; ++cnt)
    {
      Point_3 p(*in_cube++);
      Weighted_point_3 wp(p, weight);
      assert(iso_cuboid.has_on_bounded_side(p));
      assert(weight < 0.015625);
      p3rt3.insert(wp);
    }

    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == pt_count);
    assert(p3rt3.number_of_sheets() == CGAL::make_array(3,3,3) ?
        p3rt3.number_of_stored_vertices() == 27 * pt_count
        : p3rt3.number_of_stored_vertices() == pt_count);

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 0);
  }

  static void test_insert_from_file (const char* filename)
  {
    std::cout << "--- test_insert_from_file" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5));

    std::ifstream stream(filename);
    assert(stream);

    unsigned cnt = 1;
    while (stream && !(stream.eof()))
    {
      Weighted_point_3 p = read_wpoint(stream);
      std::cout << cnt << " : " << p << std::endl;
      assert(p.weight() < 0.015625);
      std::cout << "p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;
      p3rt3.insert(p);
      std::cout << "p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;
      ++cnt;
    }

    assert(p3rt3.is_valid());
  }

  static void test_insert_rt3_pointset ()
  {
    P3RT3 p3rt3(Iso_cuboid(Point_3(-100,-100,-100), Point_3(100,100,100)));

    for (int a=0;a!=10;a++)
      for (int b=0;b!=5;b++)
        for (int d=0;d!=5;d++)
        {
          Weighted_point_3 p( Point_3(a*b-d*a + (a-b)*10 +a , a-b+d +5*b, a*a-d*d+b), a*b-a*d );
          std::cout << p << std::endl;
          p3rt3.insert(p);
        }
    assert(p3rt3.is_valid());
  }

  static void test_27_to_1_sheeted_covering ()
  {
    std::cout << "--- test_27_to_1_sheeted_covering" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0, 0, 0, 1, 1, 1));

    unsigned count = 1;
    for (unsigned i = 0; i < 6; ++i)
      for (unsigned j = 0; j < 6; ++j)
        for (unsigned k = 0; k < 8; ++k)
        {
          FT x = FT(i) / FT(6);
          if (k % 2)
            x += FT(1) / FT(12);
          FT y = FT(j) / FT(6);
          if (k % 2)
            y += FT(1) / FT(12);
          FT z = FT(k) / FT(8);
          std::cout << count++ << " - " << i << " " << j << " " << k << std::endl;
          Weighted_point_3 point(Point_3(x, y, z), 0);
          p3rt3.insert(point);
          if (CGAL::make_array(i,j,k) != CGAL::make_array<unsigned>(5,5,7))
          {
            assert(p3rt3.number_of_sheets() == CGAL::make_array(3,3,3));
          }
        }

    assert(p3rt3.number_of_sheets() == CGAL::make_array(1,1,1));
    assert(p3rt3.number_of_vertices() == 6*6*8);
    assert(p3rt3.is_valid());

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());

    assert(hidden_point_count == 0);
  }

  static void test_dummy_points ()
  {
    std::cout << "--- test_dummy_points" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0, 0, 0, 1, 1, 1));

    std::vector<Vertex_handle> vertices;
    vertices.reserve(p3rt3.number_of_vertices());

    unsigned count = 1;
    for (unsigned i = 0; i < 6; ++i)
      for (unsigned j = 0; j < 6; ++j)
        for (unsigned k = 0; k < 8; ++k)
        {
          FT x = FT(i) / FT(6);
          if (k % 2)
            x += FT(1) / FT(12);
          FT y = FT(j) / FT(6);
          if (k % 2)
            y += FT(1) / FT(12);
          FT z = FT(k) / FT(8);
          std::cout << count++ << " - " << i << " " << j << " " << k << std::endl;
          Weighted_point_3 point(Point_3(x, y, z), 0);
          vertices.push_back(p3rt3.insert(point));
          if (CGAL::make_array(i,j,k) != CGAL::make_array<unsigned>(5,5,7))
          {
            assert(p3rt3.number_of_sheets() == CGAL::make_array(3,3,3));
          }
        }

    assert(p3rt3.number_of_sheets() == CGAL::make_array(1,1,1));
    assert(p3rt3.number_of_vertices() == 6*6*8);
    assert(p3rt3.is_valid());

    P3RT3 p3rt3_b(Iso_cuboid(0, 0, 0, 1, 1, 1));
    p3rt3_b.insert_dummy_points();

    assert(p3rt3_b.number_of_sheets() == CGAL::make_array(1,1,1));
    assert(p3rt3_b.number_of_vertices() == 6*6*8);
    assert(p3rt3_b.is_valid());

    assert(p3rt3 == p3rt3_b);
  }

  static void test_insert_range (unsigned pt_count, unsigned /* seed */)
  {
    std::cout << "--- test_insert_range" << std::endl;

    // the expect values for number of vertices / hidden points are hardcoded,
    // thus if the code below is used, the assertions will also need to be changed
//    CGAL::Random random(seed);
//    typedef CGAL::Creator_uniform_3<double, Point_3>  Creator;
//    CGAL::Random_points_in_cube_3<Point_3, Creator> in_cube(0.5, random);

    Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
    P3RT3 p3rt3(iso_cuboid);

    std::vector<Weighted_point_3> points;
    points.reserve(pt_count);

    std::ifstream input_stream("data/p3rt3_point_set__s7_n800");

    while (points.size() != pt_count)
    {
//      Weighted_point_3 p(*in_cube++, random.get_double(0., 0.015625));
      Weighted_point_3 p;
      input_stream >> p;
      points.push_back(p);
    }
    p3rt3.insert(points.begin(), points.end(), true);
    std::cout << "--- done " << std::endl;

    for (Vertex_iterator iter = p3rt3.vertices_begin(), end_iter = p3rt3.vertices_end(); iter != end_iter; ++iter)
    {
      typename std::vector<Weighted_point_3>::iterator it = std::find(points.begin(), points.end(), iter->point());
      assert(it != points.end());
    }
    unsigned hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    {
      for (typename Cell::Point_iterator it = iter->hidden_points_begin(), end_it = iter->hidden_points_end(); it != end_it; ++it)
      {
        assert(std::find(points.begin(), points.end(), *it) != points.end());
        ++hidden_point_count;
      }
    }

    assert(p3rt3.number_of_vertices() == 659);
    assert(p3rt3.number_of_vertices() + hidden_point_count == 800);

    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_sheets() == CGAL::make_array(1,1,1));
  }

  static void test_construction_and_insert_range (unsigned pt_count, unsigned /* seed */)
  {
    std::cout << "--- test_construction_and_insert_range" << std::endl;

    // the expect values for number of vertices / hidden points are hardcoded,
    // thus if the code below is used, the assertions will also need to be changed
//    CGAL::Random random(seed);
//    typedef CGAL::Creator_uniform_3<double, Point_3>  Creator;
//    CGAL::Random_points_in_cube_3<Point_3, Creator> in_cube(0.5, random);

    Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);

    std::vector<Weighted_point_3> points;
    points.reserve(pt_count);

    std::ifstream input_stream("data/p3rt3_point_set__s7_n800");

    while (points.size() != pt_count)
    {
//      Weighted_point_3 p(*in_cube++, random.get_double(0., 0.015625));
      Weighted_point_3 p;
      input_stream >> p;
      points.push_back(p);
    }

    P3RT3 p3rt3(points.begin(), points.end(), iso_cuboid);

    for (Vertex_iterator iter = p3rt3.vertices_begin(), end_iter = p3rt3.vertices_end(); iter != end_iter; ++iter)
    {
      typename std::vector<Weighted_point_3>::iterator it = std::find(points.begin(), points.end(), iter->point());
      assert(it != points.end());
    }
    unsigned hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
    {
      for (typename Cell::Point_iterator it = iter->hidden_points_begin(), end_it = iter->hidden_points_end(); it != end_it; ++it)
      {
        assert(std::find(points.begin(), points.end(), *it) != points.end());
        ++hidden_point_count;
      }
    }
    assert(p3rt3.number_of_vertices() == 659);
    assert(p3rt3.number_of_vertices() + hidden_point_count == 800);

    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_sheets() == CGAL::make_array(1,1,1));
  }

  static void test_locate_geometry ()
  {
    std::cout << "--- test_locate_geometry" << std::endl;

    unsigned pt_count = 600;

    CGAL::Random random(7);
    typedef CGAL::Creator_uniform_3<double, Point_3>  Creator;
    CGAL::Random_points_in_cube_3<Point_3, Creator> in_cube(0.5, random);

    Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
    P3RT3 p3rt3(iso_cuboid);

    std::vector<Weighted_point_3> points;
    points.reserve(pt_count);

    while (points.size() != pt_count)
    {
      Weighted_point_3 p(*in_cube++, random.get_double(0., 0.015625));
      points.push_back(p);
    }

    p3rt3.insert(points.begin(), points.end(), true);

    Point_3 point(-0.49, -0.49, -0.49);
    Weighted_point_3 wpoint(point);
    Offset lo;
    Cell_handle ch = p3rt3.locate(wpoint, lo);
    Periodic_tetrahedron p_tetra = p3rt3.periodic_tetrahedron(ch, lo);
    Tetrahedron tetra = p3rt3.construct_tetrahedron(p_tetra);

    assert(p3rt3.orientation(point, tetra[1], tetra[2], tetra[3]) == CGAL::POSITIVE);
    assert(p3rt3.orientation(tetra[0], point, tetra[2], tetra[3]) == CGAL::POSITIVE);
    assert(p3rt3.orientation(tetra[0], tetra[1], point, tetra[3]) == CGAL::POSITIVE);
    assert(p3rt3.orientation(tetra[0], tetra[1], tetra[2], point) == CGAL::POSITIVE);

    CGAL::Vector_3<K> v(tetra[0], tetra[1]);
    v = v * 0.5;
    point = tetra[0] + v;
    wpoint = Weighted_point_3(point);
    Locate_type lt;
    int li, lj;
    ch = p3rt3.locate(wpoint, lo, lt, li, lj);
    assert(lt == P3RT3::EDGE);
    Segment segment = p3rt3.construct_segment(p3rt3.periodic_segment(ch, lo, li, lj));
    assert(CGAL::collinear(segment[0], segment[1], point));

    v = CGAL::Vector_3<K>(point, tetra[2]);
    v = v * 0.5;
    point = point + v;
    wpoint = Weighted_point_3(point);
    ch = p3rt3.locate(wpoint, lo, lt, li, lj);
    assert(lt == P3RT3::FACET);
    Triangle triangle = p3rt3.construct_triangle(p3rt3.periodic_triangle(ch, lo, li));
    assert(p3rt3.coplanar(triangle[0], triangle[1], triangle[2], point));
  }

  static void test_number_of_hidden_points ()
  {
    std::cout << "--- test_number_of_hidden_points" << std::endl;

    P3RT3 p3rt3(Iso_cuboid(0,0,0, 1,1,1));

    Vertex_handle vh;
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.9,0.1,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.9,0.1),0.01));
    assert(vh != Vertex_handle());
    vh = p3rt3.insert(Weighted_point_3(Point_3(0.1,0.1,0.9),0.01));
    assert(vh != Vertex_handle());
    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_vertices() == 4);
    assert(p3rt3.number_of_stored_vertices() == 108);

    std::size_t hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 0);
    assert(p3rt3.number_of_hidden_points() == 0);

    vh = p3rt3.insert(Weighted_point_3(Point_3(0.101, 0.101, 0.101), 0.001));
    assert(vh == Vertex_handle());
    hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 1);
    assert(p3rt3.number_of_hidden_points() == 1);

    vh = p3rt3.insert(Weighted_point_3(Point_3(0.101, 0.101, 0.102), 0.001));
    assert(vh == Vertex_handle());
    hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 2);
    assert(p3rt3.number_of_hidden_points() == 2);

    vh = p3rt3.insert(Weighted_point_3(Point_3(0.101, 0.102, 0.101), 0.001));
    assert(vh == Vertex_handle());
    hidden_point_count = 0;
    for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
      hidden_point_count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    assert(hidden_point_count == 3);
    assert(p3rt3.number_of_hidden_points() == 3);
  }

  static void test_find_conflicts ()
  {
    std::cout << "--- test_find_conflicts" << std::endl;

    CGAL::Random random(7);
    typedef CGAL::Creator_uniform_3<double, Point_3>  Creator;
    CGAL::Random_points_in_cube_3<Point_3, Creator> in_cube(0.5, random);

    Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.6, 0.6, 0.6);
    P3RT3 p3rt3(iso_cuboid);

    std::vector<Weighted_point_3> points;
    points.reserve(800);

    for(int i=0; i<800; ++i)
    {
      Weighted_point_3 p(*in_cube++, random.get_double(0., 0.015625));
      points.push_back(p);
    }

    p3rt3.insert(points.begin(), points.end(), false);

    std::vector<Facet> bd_facets;
    std::vector<Cell_handle> conflict_cells;
    std::vector<Facet> int_facets;
    Point_3 bp(-0.5,-0.5,0.5);
    Weighted_point_3 p(bp);
    Cell_handle ch = p3rt3.locate(p);
    p3rt3.find_conflicts(p, ch,
                         std::back_inserter(bd_facets),
                         std::back_inserter(conflict_cells),
                         std::back_inserter(int_facets));
    for (unsigned int i = 0; i < bd_facets.size(); i++)
    {
      assert(
          (p3rt3.side_of_power_sphere(bd_facets[i].first, p) == CGAL::ON_BOUNDED_SIDE)
          ^ (p3rt3.side_of_power_sphere(bd_facets[i].first->neighbor(bd_facets[i].second), p)
              == CGAL::ON_BOUNDED_SIDE));
    }
    for (unsigned int i = 0; i < conflict_cells.size(); i++)
    {
      assert(p3rt3.side_of_power_sphere(conflict_cells[i], p) == CGAL::ON_BOUNDED_SIDE);
    }
    for (unsigned int i = 0; i < int_facets.size(); i++)
    {
      assert((p3rt3.side_of_power_sphere(int_facets[i].first, p) == CGAL::ON_BOUNDED_SIDE));
      assert(
          (p3rt3.side_of_power_sphere(int_facets[i].first->neighbor(int_facets[i].second), p)
              == CGAL::ON_BOUNDED_SIDE));
    }
  }

  static void test ()
  {
    test_find_conflicts();
    test_insert_range(800, 7);
    test_construction_and_insert_range(800, 7);
    test_number_of_hidden_points();
    test_locate_geometry();
    test_dummy_points();
    test_construction();
    test_insert_1();
    test_insert_point();
    test_insert_hidden_point();
    test_insert_hiding_point();
    test_insert_a_point_twice();
    test_insert_two_points_with_the_same_position();
    test_remove();
    test_27_to_1_sheeted_covering();
    //////    Iso_cuboid unitaire ->  0 <= weight < 0.015625
    test_insert_rnd_as_delaunay(100, 0.);
    test_insert_rnd_as_delaunay(100, 0.01);
  }
};

int main (int, char**)
{
  std::cout << "TESTING ..." << std::endl;

  CGAL::Set_ieee_double_precision pfr;

  std::cout << "Epeck ..." << std::endl;
  Tests<CGAL::Epeck>::test();
  std::cout << "Epick ..." << std::endl;
  Tests<CGAL::Epick>::test();

  std::cout << "EXIT SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
