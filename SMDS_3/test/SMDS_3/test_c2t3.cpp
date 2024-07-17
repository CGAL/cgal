// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : Test C3T3 class.
//******************************************************************************

#include <CGAL/Bbox_3.h>

#include "test_utilities.h"

#include <CGAL/Polyhedral_mesh_domain_3.h>

// IO
#include <fstream>
#include <iostream>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/File_tetgen.h>



template <typename K>
struct Tester
{
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
  typedef CGAL::Mesh_complex_2_in_triangulation_3<Tr, int, int> C2t3;

  typedef typename Tr::Bare_point Bare_point;
  typedef typename Tr::Weighted_point Weighted_point;

  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::FT FT;

  typedef typename C2t3::Facet Facet;
  typedef typename C2t3::Facets_in_complex_iterator Facet_iterator;
  typedef typename C2t3::Vertex_handle Vertex_handle;

  typedef typename C2t3::Vertices_in_complex Vertices_in_complex;

  typedef typename C2t3::Surface_patch_index Surface_patch_index;
  typedef typename C2t3::Index Index;

  typedef typename C2t3::size_type size_type;

  void operator()() const
  {
    //-------------------------------------------------------
    // Test default constructed c3t3
    //-------------------------------------------------------
    C2t3 c2t3;
    Tr& tr = c2t3.triangulation();

    assert(c2t3.facets_in_complex_begin() == c2t3.facets_in_complex_end());
    assert(c2t3.number_of_facets_in_complex() == 0);

    //-------------------------------------------------------
    // Data generation : fill a triangulation with 4 vertices
    //-------------------------------------------------------
    Weighted_point p1(0,0,0);
    Weighted_point p2(1,0,0);
    Weighted_point p3(0,1,0);
    Weighted_point p4(0,0,1);

    tr.insert(p1);
    tr.insert(p2);
    tr.insert(p3);
    tr.insert(p4);

    Surface_patch_index surface_patch_index (0,1);
    Surface_patch_index surface_patch_index_bis (2,3);
    Index vertex_index (2);

    //-------------------------------------------------------
    // Test empty c2t3
    //-------------------------------------------------------
    std::cerr << "\tNumber of facets in c2t3: "
              << c2t3.number_of_facets_in_complex() << std::endl;

    assert(c2t3.facets_in_complex_begin() == c2t3.facets_in_complex_end());
    assert(c2t3.number_of_facets_in_complex() == 0);

    ////-------------------------------------------------------
    //// Test move construction
    ////-------------------------------------------------------
    //C2t3 c3t3_moved{std::move(c3t3)};
    //assert(c3t3_moved.is_valid());
    //assert(c3t3.is_valid());
    //assert(ch == (Cell_handle)c3t3_moved.cells_in_complex_begin());
    //assert(c3t3_moved.number_of_cells_in_complex() == 1);
    //assert(c3t3_moved.number_of_cells_in_complex() ==
    //       (size_type)std::distance(c3t3_moved.cells_in_complex_begin(),
    //                                c3t3_moved.cells_in_complex_end()));
    //assert(c3t3_moved.is_in_complex(ch));
    //assert(c3t3_moved.subdomain_index(ch) == subdomain_index);

    //assert(c3t3.number_of_cells_in_complex() == 0);
    //assert(c3t3.number_of_cells_in_complex() == (size_type)std::distance(c3t3.cells_in_complex_begin(),
    //                                                                     c3t3.cells_in_complex_end()));
    //c3t3 = std::move(c3t3_moved);
    //assert(ch == (Cell_handle)c3t3.cells_in_complex_begin());
    //assert(c3t3.number_of_cells_in_complex() == 1);
    //assert(c3t3.number_of_cells_in_complex() == (size_type)std::distance(c3t3.cells_in_complex_begin(),
    //                                                                     c3t3.cells_in_complex_end()));
    //assert(c3t3.is_in_complex(ch));
    //assert(c3t3.subdomain_index(ch) == subdomain_index);

    //assert(c3t3_moved.number_of_cells_in_complex() == 0);
    //assert(c3t3_moved.number_of_cells_in_complex() ==
    //       (size_type)std::distance(c3t3_moved.cells_in_complex_begin(),
    //                                c3t3_moved.cells_in_complex_end()));

    //-------------------------------------------------------
    // Add facet to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert one facet in c2t3" << std::endl;

    Facet f = *( ++tr.finite_facets_begin() );
    c2t3.add_to_complex(f,surface_patch_index);

    std::cerr << "\tNumber of facets in c2t3: "
              << c2t3.number_of_facets_in_complex() << std::endl;

    assert(*(c2t3.facets_in_complex_begin()) == f);
    assert(c2t3.number_of_facets_in_complex() == 1);
    assert(c2t3.number_of_facets_in_complex() == (size_type)std::distance(c2t3.facets_in_complex_begin(),
                                                                          c2t3.facets_in_complex_end()));
    assert(c2t3.is_in_complex(f));
    assert(c2t3.surface_patch_index(f) == surface_patch_index);

    //-------------------------------------------------------
    // Remove facet from c2t3 and verify
    //-------------------------------------------------------
    std::cerr << "Remove facet from c2t3" << std::endl;

    c2t3.remove_from_complex(f);

    std::cerr << "\tNumber of facets in c2t3: "
              << c2t3.number_of_facets_in_complex() << std::endl;

    assert(c2t3.facets_in_complex_begin() == c2t3.facets_in_complex_end());
    assert(c2t3.number_of_facets_in_complex() == 0);
    assert(!c2t3.is_in_complex(f));
    assert(c2t3.surface_patch_index(f) == Surface_patch_index());

    //-------------------------------------------------------
    // Add facet to c2t3 and verify (with f=(c,i))
    //-------------------------------------------------------
    c2t3.add_to_complex(f.first,f.second,surface_patch_index);

    assert(*(c2t3.facets_in_complex_begin()) == f);
    assert(c2t3.number_of_facets_in_complex() == 1);
    assert(c2t3.number_of_facets_in_complex() == (size_type)std::distance(c2t3.facets_in_complex_begin(),
                                                                          c2t3.facets_in_complex_end()));
    assert(c2t3.is_in_complex(f));
    assert(c2t3.surface_patch_index(f) == surface_patch_index);

    c2t3.remove_from_complex(f);

    //-------------------------------------------------------
    // Add 4 facets to c2t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert 4 facets in c2t3" << std::endl;

    typename Tr::Finite_facets_iterator fit = tr.finite_facets_begin();
    while ( fit != tr.finite_facets_end() )
      c2t3.add_to_complex(*(fit++), surface_patch_index);

    std::cerr << "\tNumber of facets in c2t3: "
              << c2t3.number_of_facets_in_complex() << std::endl;

    assert(c2t3.number_of_facets_in_complex() == 4);
    assert(c2t3.number_of_facets_in_complex() == (size_type)std::distance(c2t3.facets_in_complex_begin(),
                                                                          c2t3.facets_in_complex_end()));
    //-------------------------------------------------------
    // Create c2t3_bis
    //-------------------------------------------------------
    std::cout << "Insert 6 points from domain in c2t3_bis, add 1 cell to c2t3_bis\n";
    Polyhedron polyhedron;
    std::ifstream input(CGAL::data_file_path("meshes/sphere.off"));
    input >> polyhedron;
    input.close();
    Mesh_domain domain(polyhedron);

    typedef std::vector<std::pair<Bare_point, Index> > Initial_points_vector;
    Initial_points_vector initial_points;
    domain.construct_initial_points_object()(std::back_inserter(initial_points), 6);

    //-------------------------------------------------------
    // Swap c2t3 and c2t3_bis
    //-------------------------------------------------------
    std::cout << "Swap c2t3 and c2t3_bis\n";
    C2t3 c2t3_bis;
    typedef typename C2t3::size_type size_type;

    size_type c2t3_facet_nb = c2t3.number_of_facets_in_complex();
    size_type c2t3_vertex_nb = c2t3.triangulation().number_of_vertices();

    size_type c2t3_bis_facet_nb = c2t3_bis.number_of_facets_in_complex();
    size_type c2t3_bis_vertex_nb = c2t3_bis.triangulation().number_of_vertices();

    c2t3.swap(c2t3_bis);

    std::cout << "\tNumber of facets in c2t3: "
              << c2t3.number_of_facets_in_complex() << std::endl;
    std::cout << "\tNumber of vertices in c2t3 triangulation: "
              << c2t3.triangulation().number_of_vertices() << std::endl;

    std::cout << "\tNumber of facets in c2t3_bis: "
              << c2t3_bis.number_of_facets_in_complex() << std::endl;
    std::cout << "\tNumber of vertices in c2t3_bis triangulation: "
              << c2t3_bis.triangulation().number_of_vertices() << std::endl;

    assert(c2t3_facet_nb == c2t3_bis.number_of_facets_in_complex());
    assert(c2t3_vertex_nb == c2t3_bis.triangulation().number_of_vertices());

    assert(c2t3_bis_facet_nb == c2t3.number_of_facets_in_complex());
    assert(c2t3_bis_vertex_nb == c2t3.triangulation().number_of_vertices());

    // reset
    c2t3.swap(c2t3_bis);

    //-------------------------------------------------------
    // Modify indices and dimension and verify
    //-------------------------------------------------------
    std::cerr << "Play with indices\n";
    Facet f2 = *c2t3.facets_in_complex_begin();
    Vertex_handle vh = c2t3.triangulation().vertices_begin();

    c2t3.set_surface_patch_index(f2, surface_patch_index_bis);

    c2t3.set_dimension(vh, 1);
    c2t3.set_index(vh, vertex_index);

    assert(c2t3.surface_patch_index(f2) == surface_patch_index_bis);
    assert(c2t3.in_dimension(vh) == 1);
    assert(c2t3.index(vh) == vertex_index);

    c2t3.set_surface_patch_index(f2.first, f2.second, surface_patch_index);
    assert(c2t3.surface_patch_index(f2) == surface_patch_index);

    // -----------------------------------
    // Test surface patch iterators
    // -----------------------------------
    std::cerr << "Test surface patch iterators\n";
    c2t3.set_surface_patch_index(f2.first, f2.second, surface_patch_index_bis);

    typename C2t3::Facets_in_complex_iterator patch_fit =
      c2t3.facets_in_complex_begin(surface_patch_index);
    typename C2t3::Facets_in_complex_iterator patch_fit_bis =
      c2t3.facets_in_complex_begin(surface_patch_index_bis);
    typename C2t3::Facets_in_complex_iterator fend =
      c2t3.facets_in_complex_end();

    std::cout << "\tNumber of facets of index '<" << surface_patch_index.first
              << "," << surface_patch_index.second << ">': "
              << std::distance(patch_fit,fend) << std::endl;
    std::cout << "\tNumber of facets of index '<" << surface_patch_index_bis.first
              << "," << surface_patch_index.second << ">': "
              << std::distance(patch_fit_bis,fend) << std::endl;

    assert ( std::distance(patch_fit,fend) == 3 );
    assert ( std::distance(patch_fit_bis,fend) == 1 );
    assert ( c2t3.surface_patch_index(*patch_fit) == surface_patch_index );
    assert ( c2t3.surface_patch_index(*patch_fit_bis) == surface_patch_index_bis );
  }
};


int main()
{
  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n";
  Tester<K_e_i> test_epic;
  test_epic();

  std::cerr << "\nTESTING WITH Exact_predicates_exact_constructions_kernel...\n";
  Tester<K_e_e> test_epec;
  test_epec();

  return EXIT_SUCCESS;
}
