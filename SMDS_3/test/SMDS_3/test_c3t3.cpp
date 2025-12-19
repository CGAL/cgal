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
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  typedef typename Tr::Bare_point Bare_point;
  typedef typename Tr::Weighted_point Weighted_point;

  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::FT FT;

  typedef typename C3t3::Cell_handle Cell_handle;
  typedef typename C3t3::Facet Facet;
  typedef typename C3t3::Facets_in_complex_iterator Facet_iterator;
  typedef typename C3t3::Vertex_handle Vertex_handle;

  typedef typename C3t3::Vertices_in_complex Vertices_in_complex;
  typedef typename C3t3::Cells_in_complex Cells_in_complex;

  typedef typename C3t3::Subdomain_index Subdomain_index;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename C3t3::Index Index;

  typedef typename C3t3::size_type size_type;

  void operator()() const
  {
    //-------------------------------------------------------
    // Test default constructed c3t3
    //-------------------------------------------------------
    C3t3 c3t3;
    Tr& tr = c3t3.triangulation();

    assert(c3t3.cells_in_complex_begin() == c3t3.cells_in_complex_end());
    assert(c3t3.facets_in_complex_begin() == c3t3.facets_in_complex_end());
    assert(c3t3.number_of_cells_in_complex() == 0);
    assert(c3t3.number_of_facets_in_complex() == 0);

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

    Subdomain_index subdomain_index (1);
    Subdomain_index subdomain_index_bis (2);
    Surface_patch_index surface_patch_index (0,1);
    Surface_patch_index surface_patch_index_bis (2,3);
    Index vertex_index (2);

    //-------------------------------------------------------
    // Test empty c3t3
    //-------------------------------------------------------
    std::cerr << "\tNumber of cells in c3t3: "
              << c3t3.number_of_cells_in_complex() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: "
              << c3t3.number_of_facets_in_complex() << std::endl;

    assert(c3t3.cells_in_complex_begin() == c3t3.cells_in_complex_end());
    assert(c3t3.facets_in_complex_begin() == c3t3.facets_in_complex_end());
    assert(c3t3.number_of_cells_in_complex() == 0);
    assert(c3t3.number_of_facets_in_complex() == 0);

    //-------------------------------------------------------
    // Add cell to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert one cell in c3t3" << std::endl;

    Cell_handle ch = tr.finite_cells_begin();
    c3t3.add_to_complex(ch,subdomain_index);

    std::cerr << "\tNumber of cells in c3t3: "
              << c3t3.number_of_cells_in_complex() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: "
              << c3t3.number_of_facets_in_complex() << std::endl;

    assert(ch == (Cell_handle)c3t3.cells_in_complex_begin());
    assert(c3t3.number_of_cells_in_complex() == 1);
    assert(c3t3.number_of_cells_in_complex() == (size_type)std::distance(c3t3.cells_in_complex_begin(),
                                                                         c3t3.cells_in_complex_end()));
    assert(c3t3.is_in_complex(ch));
    assert(c3t3.subdomain_index(ch) == subdomain_index);

    // Test iterator range
    {
      Cell_handle ch = *c3t3.cells_in_complex().begin();
      for (auto c : c3t3.cells_in_complex()) {
        assert(c == ch);
        break;
      }
    }

    //-------------------------------------------------------
    // Test move construction
    //-------------------------------------------------------
    C3t3 c3t3_moved{std::move(c3t3)};
    assert(c3t3_moved.is_valid());
    assert(c3t3.is_valid());
    assert(ch == (Cell_handle)c3t3_moved.cells_in_complex_begin());
    assert(c3t3_moved.number_of_cells_in_complex() == 1);
    assert(c3t3_moved.number_of_cells_in_complex() ==
           (size_type)std::distance(c3t3_moved.cells_in_complex_begin(),
                                    c3t3_moved.cells_in_complex_end()));
    assert(c3t3_moved.is_in_complex(ch));
    assert(c3t3_moved.subdomain_index(ch) == subdomain_index);

    assert(c3t3.number_of_cells_in_complex() == 0);
    assert(c3t3.number_of_cells_in_complex() == (size_type)std::distance(c3t3.cells_in_complex_begin(),
                                                                         c3t3.cells_in_complex_end()));
    c3t3 = std::move(c3t3_moved);
    assert(ch == (Cell_handle)c3t3.cells_in_complex_begin());
    assert(c3t3.number_of_cells_in_complex() == 1);
    assert(c3t3.number_of_cells_in_complex() == (size_type)std::distance(c3t3.cells_in_complex_begin(),
                                                                         c3t3.cells_in_complex_end()));
    assert(c3t3.is_in_complex(ch));
    assert(c3t3.subdomain_index(ch) == subdomain_index);

    assert(c3t3_moved.number_of_cells_in_complex() == 0);
    assert(c3t3_moved.number_of_cells_in_complex() ==
           (size_type)std::distance(c3t3_moved.cells_in_complex_begin(),
                                    c3t3_moved.cells_in_complex_end()));
    // -----------------------------------
    // Test Cell_in_complex_iterator
    // The goal here is to test operators and conversion on iterator type
    // -----------------------------------
    typename C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();
    ch = cit;
    typename C3t3::Triangulation::Cell& c1 = *ch;
    typename C3t3::Triangulation::Cell& c2 = *cit;

    assert( c1.subdomain_index() == c2.subdomain_index() );
    assert ( cit->vertex(0) == ch->vertex(0) );

    //-------------------------------------------------------
    // Remove cell from c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Remove cell from c3t3" << std::endl;

    c3t3.remove_from_complex(ch);

    std::cerr << "\tNumber of cells in c3t3: "
              << c3t3.number_of_cells_in_complex() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: "
              << c3t3.number_of_facets_in_complex() << std::endl;

    assert(c3t3.number_of_cells_in_complex() == 0);
    assert(! c3t3.is_in_complex(ch));
    assert(c3t3.subdomain_index(ch) == Subdomain_index());

    //-------------------------------------------------------
    // Add facet to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert one facet in c3t3" << std::endl;

    Facet f = *( ++tr.finite_facets_begin() );
    c3t3.add_to_complex(f,surface_patch_index);

    std::cerr << "\tNumber of cells in c3t3: "
              << c3t3.number_of_cells_in_complex() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: "
              << c3t3.number_of_facets_in_complex() << std::endl;

    assert(*(c3t3.facets_in_complex_begin()) == f);
    assert(c3t3.number_of_facets_in_complex() == 1);
    assert(c3t3.number_of_facets_in_complex() == (size_type)std::distance(c3t3.facets_in_complex_begin(),
                                                                          c3t3.facets_in_complex_end()));
    assert(c3t3.is_in_complex(f));
    assert(c3t3.surface_patch_index(f) == surface_patch_index);

    //-------------------------------------------------------
    // Remove facet from c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Remove facet from c3t3" << std::endl;

    c3t3.remove_from_complex(f);

    std::cerr << "\tNumber of cells in c3t3: "
              << c3t3.number_of_cells_in_complex() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: "
              << c3t3.number_of_facets_in_complex() << std::endl;

    assert(c3t3.facets_in_complex_begin() == c3t3.facets_in_complex_end());
    assert(c3t3.number_of_facets_in_complex() == 0);
    assert(!c3t3.is_in_complex(f));
    assert(c3t3.surface_patch_index(f) == Surface_patch_index());

    //-------------------------------------------------------
    // Add facet to c3t3 and verify (with f=(c,i))
    //-------------------------------------------------------
    c3t3.add_to_complex(f.first,f.second,surface_patch_index);

    assert(*(c3t3.facets_in_complex_begin()) == f);
    assert(c3t3.number_of_facets_in_complex() == 1);
    assert(c3t3.number_of_facets_in_complex() == (size_type)std::distance(c3t3.facets_in_complex_begin(),
                                                                          c3t3.facets_in_complex_end()));
    assert(c3t3.is_in_complex(f));
    assert(c3t3.surface_patch_index(f) == surface_patch_index);

    c3t3.remove_from_complex(f);

    //-------------------------------------------------------
    // Add 4 facets to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert 4 facets in c3t3" << std::endl;

    typename Tr::Finite_facets_iterator fit = tr.finite_facets_begin();
    while ( fit != tr.finite_facets_end() )
      c3t3.add_to_complex(*(fit++), surface_patch_index);

    std::cerr << "\tNumber of cells in c3t3: "
              << c3t3.number_of_cells_in_complex() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: "
              << c3t3.number_of_facets_in_complex() << std::endl;

    assert(c3t3.number_of_facets_in_complex() == 4);
    assert(c3t3.number_of_facets_in_complex() == (size_type)std::distance(c3t3.facets_in_complex_begin(),
                                                                          c3t3.facets_in_complex_end()));
    //-------------------------------------------------------
    // Create c3t3_bis
    //-------------------------------------------------------
    std::cout << "Insert 6 points from domain in c3t3_bis, add 1 cell to c3t3_bis\n";
    Polyhedron polyhedron;
    std::ifstream input(CGAL::data_file_path("meshes/sphere.off"));
    input >> polyhedron;
    input.close();
    Mesh_domain domain(polyhedron);

    typedef std::vector<std::pair<Bare_point, Index> > Initial_points_vector;
    Initial_points_vector initial_points;
    domain.construct_initial_points_object()(std::back_inserter(initial_points), 6);

    C3t3 c3t3_bis;
    c3t3_bis.insert_surface_points(initial_points.begin(), initial_points.end());

    Cell_handle ch_bis = c3t3_bis.triangulation().finite_cells_begin();
    c3t3_bis.add_to_complex(ch_bis,subdomain_index);

    std::cout << "\tNumber of cells in c3t3_bis: "
              << c3t3_bis.number_of_cells_in_complex() << std::endl;
    std::cout << "\tNumber of facets in c3t3_bis: "
              << c3t3_bis.number_of_facets_in_complex() << std::endl;
    std::cout << "\tNumber of vertices in c3t3_bis triangulation: "
              << c3t3_bis.triangulation().number_of_vertices() << std::endl;

    //-------------------------------------------------------
    // Swap c3t3 and c3t3_bis
    //-------------------------------------------------------
    std::cout << "Swap c3t3 and c3t3_bis\n";
    typedef typename C3t3::size_type size_type;

    size_type c3t3_cell_nb = c3t3.number_of_cells_in_complex();
    size_type c3t3_facet_nb = c3t3.number_of_facets_in_complex();
    size_type c3t3_vertex_nb = c3t3.triangulation().number_of_vertices();

    size_type c3t3_bis_cell_nb = c3t3_bis.number_of_cells_in_complex();
    size_type c3t3_bis_facet_nb = c3t3_bis.number_of_facets_in_complex();
    size_type c3t3_bis_vertex_nb = c3t3_bis.triangulation().number_of_vertices();

    c3t3.swap(c3t3_bis);

    std::cout << "\tNumber of cells in c3t3: "
              << c3t3.number_of_cells_in_complex() << std::endl;
    std::cout << "\tNumber of facets in c3t3: "
              << c3t3.number_of_facets_in_complex() << std::endl;
    std::cout << "\tNumber of vertices in c3t3 triangulation: "
              << c3t3.triangulation().number_of_vertices() << std::endl;

    std::cout << "\tNumber of cells in c3t3_bis: "
              << c3t3_bis.number_of_cells_in_complex() << std::endl;
    std::cout << "\tNumber of facets in c3t3_bis: "
              << c3t3_bis.number_of_facets_in_complex() << std::endl;
    std::cout << "\tNumber of vertices in c3t3_bis triangulation: "
              << c3t3_bis.triangulation().number_of_vertices() << std::endl;

    assert(c3t3_cell_nb == c3t3_bis.number_of_cells_in_complex());
    assert(c3t3_facet_nb == c3t3_bis.number_of_facets_in_complex());
    assert(c3t3_vertex_nb == c3t3_bis.triangulation().number_of_vertices());

    assert(c3t3_bis_cell_nb == c3t3.number_of_cells_in_complex());
    assert(c3t3_bis_facet_nb == c3t3.number_of_facets_in_complex());
    assert(c3t3_bis_vertex_nb == c3t3.triangulation().number_of_vertices());

    // reset
    c3t3.swap(c3t3_bis);

    //-------------------------------------------------------
    // Modify indices and dimension and verify
    //-------------------------------------------------------
    std::cerr << "Play with indices\n";
    c3t3.add_to_complex(ch_bis, subdomain_index);
    Facet f2 = *c3t3.facets_in_complex_begin();
    Vertex_handle vh = c3t3.triangulation().vertices_begin();

    c3t3.set_subdomain_index(ch, subdomain_index_bis);
    c3t3.set_surface_patch_index(f2, surface_patch_index_bis);
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
    c3t3.set_surface_index(f2, surface_patch_index_bis);
    c3t3.set_surface_index(f2.first, f2.second, surface_patch_index_bis);
#endif
    c3t3.set_dimension(vh, 1);
    c3t3.set_index(vh, vertex_index);

    assert(c3t3.subdomain_index(ch) == subdomain_index_bis);
    assert(c3t3.surface_patch_index(f2) == surface_patch_index_bis);
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
    assert(c3t3.surface_index(f2) == surface_patch_index_bis);
    assert(c3t3.surface_index(f2.first, f2.second) == surface_patch_index_bis);
#endif
    assert(c3t3.in_dimension(vh) == 1);
    assert(c3t3.index(vh) == vertex_index);

    c3t3.set_surface_patch_index(f2.first, f2.second, surface_patch_index);
    assert(c3t3.surface_patch_index(f2) == surface_patch_index);

    // -----------------------------------
    // Test subdomain iterators
    // -----------------------------------
    std::cerr << "Test subdomain iterators\n";

    typename C3t3::Cells_in_complex_iterator subdomain_cit =
      c3t3.cells_in_complex_begin(subdomain_index);
    typename C3t3::Cells_in_complex_iterator subdomain_cit_bis =
      c3t3.cells_in_complex_begin(subdomain_index_bis);
    typename C3t3::Cells_in_complex_iterator cend =
      c3t3.cells_in_complex_end();

    std::cout << "\tNumber of cells of index '" << subdomain_index << "': "
              << std::distance(subdomain_cit,cend) << std::endl;
    std::cout << "\tNumber of cells of index '" << subdomain_index_bis << "': "
              << std::distance(subdomain_cit_bis,cend) << std::endl;

    assert ( std::distance(subdomain_cit,cend) == 0 );
    assert ( std::distance(subdomain_cit_bis,cend) == 1 );
    assert ( c3t3.subdomain_index(subdomain_cit_bis) == subdomain_index_bis );

    // -----------------------------------
    // Test surface patch iterators
    // -----------------------------------
    std::cerr << "Test surface patch iterators\n";
    c3t3.set_surface_patch_index(f2.first, f2.second, surface_patch_index_bis);

    typename C3t3::Facets_in_complex_iterator patch_fit =
      c3t3.facets_in_complex_begin(surface_patch_index);
    typename C3t3::Facets_in_complex_iterator patch_fit_bis =
      c3t3.facets_in_complex_begin(surface_patch_index_bis);
    typename C3t3::Facets_in_complex_iterator fend =
      c3t3.facets_in_complex_end();

    std::cout << "\tNumber of facets of index '<" << surface_patch_index.first
              << "," << surface_patch_index.second << ">': "
              << std::distance(patch_fit,fend) << std::endl;
    std::cout << "\tNumber of facets of index '<" << surface_patch_index_bis.first
              << "," << surface_patch_index.second << ">': "
              << std::distance(patch_fit_bis,fend) << std::endl;

    assert ( std::distance(patch_fit,fend) == 3 );
    assert ( std::distance(patch_fit_bis,fend) == 1 );
    assert ( c3t3.surface_patch_index(*patch_fit) == surface_patch_index );
    assert ( c3t3.surface_patch_index(*patch_fit_bis) == surface_patch_index_bis );

    std::ofstream out_medit("test-medit.mesh");
    CGAL::IO::write_MEDIT(out_medit, c3t3);
    out_medit.close();
    CGAL::IO::output_to_tetgen("test-tetgen", c3t3);
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
