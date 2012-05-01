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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

//#define CGAL_MESH_3_DEBUG_CELL_CRITERIA
//#define CGAL_MESH_3_DEBUG_FACET_CRITERIA

#include "test_utilities.h"
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_3/Creator_weighted_point_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Mesh_criteria_3.h>

#include <string>
#include <iostream>
#include <sstream>


template <typename K>
struct Tester
{
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_traits;

  typedef typename CGAL::Mesh_triangulation_3<Mesh_traits>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
  typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
  typedef typename Mesh_criteria::Cell_criteria Cell_criteria;

  typedef typename C3t3::Cell_handle Cell_handle;
  typedef typename C3t3::Facet Facet;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;

  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::FT FT;
  typedef typename Gt::Point_3 Point;
  typedef CGAL::Mesh_3::Creator_weighted_point_3<FT, Point> Point_creator;

  C3t3 c3t3_;
  // Cells & facets
  int c1_,c2_,f1_,f2_;
  // Facets center
  Point p1_,p2_;

  Tester()
    : c3t3_()
    , c1_(0)
    , c2_(1)
    , f1_(0)
    , f2_(4)
  {
    Tr& tr = c3t3_.triangulation();

    Point_creator creator;
    Point p1 = creator(0,0,0);
    Point p2 = creator(1,0,0);
    Point p3 = creator(0,1,0);
    Point p4 = creator(0,0,100);
    Point p5 = creator(0,0,-1);

    p1_ = creator(0,.5,50);
    p2_ = creator(1,1,-1);

    tr.insert(p1);
    tr.insert(p2);
    tr.insert(p3);
    tr.insert(p4);
    tr.insert(p5);
  }

  void init_facet(const Facet& f, const Point& p) const
  {
    f.first->set_surface_patch_index(f.second,Surface_patch_index(0,1));
    f.first->set_facet_surface_center(f.second,p);

    for ( int i=0; i<4 ; ++i )
    {
      if ( i != f.second )
        f.first->vertex(i)->set_dimension(2);
    }
  }

  void init_facets() const
  {
    init_facet(facet_handle(f1_),p1_);
    init_facet(facet_handle(f2_),p2_);
  }

  void remove_surface_vertex() const
  {
    // Removes surface vertex on f1. f1 & f2 share one vertex, so it must be
    // sufficient to have facet_on_surface_criterion false for f1 & f2
    Facet f = facet_handle(f1_);
    for ( int i=0; i<4 ; ++i )
    {
      if ( i != f.second )
        f.first->vertex(i)->set_dimension(3);
    }
  }

  void print(const Cell_handle c, const std::string& name) const
  {
    std::cerr << name << ":[" ;

    for ( int i=0; i<4 ; ++i )
      std::cerr << "[" << c->vertex(i)->point() << "]";

    std::cerr << "] ";
  }

  void print(const Facet& f, const std::string& name) const
  {
    std::cerr << name << ":[" ;

    for ( int i=0; i<4 ; ++i )
    {
      if ( i != f.second )
        std::cerr << "[" << f.first->vertex(i)->point() << "]";
    }

    std::cerr << "] ";
  }

  void print_cells() const
  {
    Cell_handle c1 = cell_handle(c1_);
    Cell_handle c2 = cell_handle(c2_);

    std::cerr << "Testing using cells: ";
    print(c1,"c1");
    print(c2,"c2");
    std::cerr << std::endl;
  }

  void print_facets() const
  {
//    for ( int i = 0 ; i < 12 ; ++i)
//    {
//      Facet fi = facet(i);
//      std::stringstream s;
//      s << i;
//      print(fi,s.str());
//    }

    Facet f1 = facet_handle(f1_);
    Facet f2 = facet_handle(f2_);

    std::cerr << "Testing using facets: ";
    print(f1,"f1");
    print(f2,"f2");
    std::cerr << std::endl;
  }


  Cell_handle cell_handle(int i=0) const
  {
    typedef typename Tr::Finite_cells_iterator Iterator;
    Iterator it = c3t3_.triangulation().finite_cells_begin();

    int j = 0;
    while ( j<i && it != c3t3_.triangulation().finite_cells_end() )
    {
      ++it;
      ++j;
    }

    if ( it == c3t3_.triangulation().finite_cells_end() )
      --it;

    return it;
  }


  Facet facet_handle(int i=0) const
  {
    typedef typename Tr::Finite_facets_iterator Iterator;
    Iterator it = c3t3_.triangulation().finite_facets_begin();

    int j = 0;
    while ( j<i && it != c3t3_.triangulation().finite_facets_end() )
    {
      ++it;
      ++j;
    }

    if ( it == c3t3_.triangulation().finite_facets_end() )
      --it;

    return *it;
  }


  void test_cell(const FT radius,
                 const FT radius_edge,
                 const bool is_cell1_bad = false,
                 const bool is_cell2_bad = false,
                 const bool compare_cells = false) const
  {
    typedef typename Cell_criteria::Cell_badness Badness;

    Cell_handle cell1 = cell_handle(c1_);
    Cell_handle cell2 = cell_handle(c2_);

    Cell_criteria cell_criteria(radius_edge,radius);

    Badness b1 = cell_criteria(cell1);
    Badness b2 = cell_criteria(cell2);

    std::cerr << "\t[Radius bound: " << radius
              << " - Radius-edge bound: " << radius_edge
              << "]\tc1 is ";

    if ( b1 )
      std::cerr << "BAD q=<" << b1->first << ";" << b1->second << ">";
    else
      std::cerr << "GOOD ";

    std::cerr << "\t\tc2 is ";
    if ( b2 )
      std::cerr << "BAD q=<" << b2->first << ";" << b2->second << ">";
    else
      std::cerr << "GOOD";

    assert( is_cell1_bad == (bool)b1 );
    assert( is_cell2_bad == (bool)b2 );

    // b1 is badder than b2
    if ( compare_cells )
    {
      std::cerr << "\t\tq(c1) < q(c2)";
      assert(*b1<*b2);
    }

    std::cerr << std::endl;
  }

  void test_facet(const FT angle,
                  const FT radius,
                  const FT distance,
                  const bool is_facet1_bad = false,
                  const bool is_facet2_bad = false,
                  const bool compare_facets = false) const
    {
      typedef typename Facet_criteria::Facet_badness Badness;

      Facet f1 = facet_handle(f1_);
      Facet f2 = facet_handle(f2_);

      Facet_criteria criteria(angle, radius, distance);

      Badness b1 = criteria(f1);
      Badness b2 = criteria(f2);

      std::cerr << "\t[Angle bound: " << angle
                << " - Radius bound: " << radius
                << " - Distance bound: " << distance
                << "]   \tf1 is ";

      if ( b1 )
        std::cerr << "BAD q=<" << b1->first << ";" << b1->second << ">";
      else
        std::cerr << "GOOD ";

      std::cerr << "\t\tf2 is ";
      if ( b2 )
        std::cerr << "BAD q=<" << b2->first << ";" << b2->second << ">";
      else
        std::cerr << "GOOD";

      assert( is_facet1_bad == (bool)b1 );
      assert( is_facet2_bad == (bool)b2 );

      // b1 is badder than b2
      if ( compare_facets )
      {
        std::cerr << "\t\tq(f1) < q(f2)";
        assert(*b1<*b2);
      }

      std::cerr << std::endl;
    }


  void test_cell() const
  {
    print_cells();

    // Test with different values and known results
    test_cell(0,0);
    test_cell(49,0,true);
    test_cell(51,0);
    test_cell(0,51);
    test_cell(0,49,true);
    test_cell(.5,0,true,true,true);
    test_cell(0,.5,true,true,true);
    test_cell(.5,20,true,true,true);
    test_cell(20,.5,true,true,true);
  }

  void test_facet() const
  {
    init_facets();
    print_facets();

    // Test with different values and known results
    test_facet(59.5,0,0,true);
    test_facet(60.5,0,0,true,true,true);
    test_facet(0,49,0,true);
    test_facet(0,51,0);
    // distance from center(f2) to p2 is 2/sqrt(3)=1.15..
    test_facet(0,0,1.12,false,true);
    test_facet(0,0,1.16);
    // Test with various criterion
    test_facet(30,5,1,true,true,true);
    // Test facet on surface criterion
    remove_surface_vertex();
    test_facet(0,0,0,true,true);
    test_facet(0.5,5,5,true,true,true);
  }
};

int main()
{
  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n";
  Tester<K_e_i> test_epic;
  test_epic.test_cell();
  test_epic.test_facet();

  std::cerr << "\nTESTING WITH Exact_predicates_exact_constructions_kernel...\n";
  Tester<K_e_e> test_epec;
  test_epec.test_cell();
  test_epec.test_facet();

  return EXIT_SUCCESS;
}
