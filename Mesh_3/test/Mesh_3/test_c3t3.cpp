// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// File Description : Test C3T3 class.
//******************************************************************************

#include <CGAL/Bbox_3.h>

#include "test_utilities.h"
#include <CGAL/Mesh_3/Creator_weighted_point_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>

// IO
#include <fstream>
#include <iostream>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_medit.h>



template <typename K>
struct Tester
{
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::FT FT;
  typedef typename Gt::Point_3 Point;
  typedef CGAL::Mesh_3::Creator_weighted_point_3<FT, Point> Point_creator;

  typedef typename C3t3::Cell_handle Cell_handle;
  typedef typename C3t3::Facet Facet;
  typedef typename C3t3::Facet_iterator Facet_iterator;
  typedef typename C3t3::Vertex_handle Vertex_handle;

  typedef typename C3t3::Subdomain_index Subdomain_index;
  typedef typename C3t3::Surface_index Surface_index;
  typedef typename C3t3::Index Index;

  typedef typename C3t3::size_type size_type;

  void operator()() const
  {
    //-------------------------------------------------------
    // Test default constructed c3t3
    //-------------------------------------------------------
    C3t3 c3t3;
    Tr& tr = c3t3.triangulation();

    assert(c3t3.cells_begin() == c3t3.cells_end());
    assert(c3t3.facets_begin() == c3t3.facets_end());
    assert(c3t3.number_of_cells() == 0);
    assert(c3t3.number_of_facets() == 0);

    //-------------------------------------------------------
    // Data generation : fill a triangulation with 4 vertices
    //-------------------------------------------------------
    Point_creator creator;
    Point p1 = creator(0,0,0);
    Point p2 = creator(1,0,0);
    Point p3 = creator(0,1,0);
    Point p4 = creator(0,0,1);

    tr.insert(p1);
    tr.insert(p2);
    tr.insert(p3);
    tr.insert(p4);

    Subdomain_index subdomain_index (1);
    Subdomain_index subdomain_index_bis (2);
    Surface_index surface_index (0,1);
    Surface_index surface_index_bis (2,3);
    Index vertex_index (2);

    //-------------------------------------------------------
    // Test empty c3t3
    //-------------------------------------------------------
    std::cerr << "\tNumber of cells in c3t3: " << c3t3.number_of_cells() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: " << c3t3.number_of_facets() << std::endl;

    assert(c3t3.cells_begin() == c3t3.cells_end());
    assert(c3t3.facets_begin() == c3t3.facets_end());
    assert(c3t3.number_of_cells() == 0);
    assert(c3t3.number_of_facets() == 0);

    //-------------------------------------------------------
    // Add cell to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert one cell in c3t3" << std::endl;

    Cell_handle ch = tr.finite_cells_begin();
    c3t3.add_to_complex(ch,subdomain_index);

    std::cerr << "\tNumber of cells in c3t3: " << c3t3.number_of_cells() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: " << c3t3.number_of_facets() << std::endl;

    assert(ch == (Cell_handle)c3t3.cells_begin());
    assert(c3t3.number_of_cells() == 1);
    assert(c3t3.number_of_cells() == (size_type)std::distance(c3t3.cells_begin(),
                                                              c3t3.cells_end()));
    assert(c3t3.is_in_complex(ch));
    assert(c3t3.subdomain_index(ch) == subdomain_index);

    //-------------------------------------------------------
    // Remove cell from c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Remove cell from c3t3" << std::endl;

    c3t3.remove_from_complex(ch);

    std::cerr << "\tNumber of cells in c3t3: " << c3t3.number_of_cells() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: " << c3t3.number_of_facets() << std::endl;

    assert(c3t3.number_of_cells() == 0);
    assert(! c3t3.is_in_complex(ch));
    assert(c3t3.subdomain_index(ch) == Subdomain_index());

    //-------------------------------------------------------
    // Add facet to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert one facet in c3t3" << std::endl;

    Facet f = *( ++tr.finite_facets_begin() );
    c3t3.add_to_complex(f,surface_index);

    std::cerr << "\tNumber of cells in c3t3: " << c3t3.number_of_cells() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: " << c3t3.number_of_facets() << std::endl;

    assert(*(c3t3.facets_begin()) == f);
    assert(c3t3.number_of_facets() == 1);
    assert(c3t3.number_of_facets() == (size_type)std::distance(c3t3.facets_begin(),
                                                               c3t3.facets_end()));
    assert(c3t3.is_in_complex(f));
    assert(c3t3.surface_index(f) == surface_index);

    //-------------------------------------------------------
    // Remove facet from c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Remove facet from c3t3" << std::endl;

    c3t3.remove_from_complex(f);

    std::cerr << "\tNumber of cells in c3t3: " << c3t3.number_of_cells() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: " << c3t3.number_of_facets() << std::endl;

    assert(c3t3.facets_begin() == c3t3.facets_end());
    assert(c3t3.number_of_facets() == 0);
    assert(!c3t3.is_in_complex(f));
    assert(c3t3.surface_index(f) == Surface_index());

    //-------------------------------------------------------
    // Add facet to c3t3 and verify (with f=(c,i))
    //-------------------------------------------------------
    c3t3.add_to_complex(f.first,f.second,surface_index);

    assert(*(c3t3.facets_begin()) == f);
    assert(c3t3.number_of_facets() == 1);
    assert(c3t3.number_of_facets() == (size_type)std::distance(c3t3.facets_begin(),
                                                               c3t3.facets_end()));
    assert(c3t3.is_in_complex(f));
    assert(c3t3.surface_index(f) == surface_index);

    c3t3.remove_from_complex(f);

    //-------------------------------------------------------
    // Add 4 facets to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert 4 facets in c3t3" << std::endl;

    typename Tr::Finite_facets_iterator fit = tr.finite_facets_begin();
    while ( fit != tr.finite_facets_end() )
      c3t3.add_to_complex(*(fit++), surface_index);

    std::cerr << "\tNumber of cells in c3t3: " << c3t3.number_of_cells() << std::endl;
    std::cerr << "\tNumber of facets in c3t3: " << c3t3.number_of_facets() << std::endl;

    assert(c3t3.number_of_facets() == 4);
    assert(c3t3.number_of_facets() == (size_type)std::distance(c3t3.facets_begin(),
                                                               c3t3.facets_end()));
    //-------------------------------------------------------
    // Test I/O operators
    //-------------------------------------------------------
    {
    std::cout << "Test I/O operators in ascii mode.\n" 
              << "Save c3t3 to \"c3t3_dump\"..." << std::endl;
    std::ofstream out_file("c3t3.dump.ascii");
    out_file << c3t3;
    assert(out_file.good());
    out_file.close();

    std::cout << "Then reload this file to c3t3_reload..." << std::endl;
    std::ifstream in_file("c3t3.dump.ascii");
    C3t3 c3t3_loaded;
    in_file >> c3t3_loaded;
    assert(in_file.good());
    in_file.close();

    std::cerr << "\tNumber of cells in c3t3_loaded: " 
              << c3t3_loaded.number_of_cells() << std::endl;
    std::cerr << "\tNumber of facets in c3t3_loaded: " 
              << c3t3_loaded.number_of_facets() << std::endl;
    assert( c3t3_loaded.number_of_cells() ==  c3t3.number_of_cells() );
    assert( c3t3_loaded.number_of_facets() ==  c3t3.number_of_facets() );
    }
    if(!boost::is_same<K, K_e_e>::value) {
    std::cout << "Test I/O operators in binary mode.\n" 
              << "Save c3t3 to \"c3t3_dump\"..." << std::endl;
    std::ofstream out_file("c3t3.dump",
                           std::ios_base::out|std::ios_base::binary);
    CGAL::set_binary_mode(out_file);
    out_file << c3t3;
    assert(out_file.good());
    out_file.close();

    std::cout << "Then reload this file to c3t3_reload..." << std::endl;
    std::ifstream in_file("c3t3.dump",
                          std::ios_base::in|std::ios_base::binary);
    CGAL::set_binary_mode(in_file);
    C3t3 c3t3_loaded;
    in_file >> c3t3_loaded;
    assert(in_file.good());
    in_file.close();

    std::cerr << "\tNumber of cells in c3t3_loaded: " 
              << c3t3_loaded.number_of_cells() << std::endl;
    std::cerr << "\tNumber of facets in c3t3_loaded: " 
              << c3t3_loaded.number_of_facets() << std::endl;
    assert( c3t3_loaded.number_of_cells() ==  c3t3.number_of_cells() );
    assert( c3t3_loaded.number_of_facets() ==  c3t3.number_of_facets() );
    }
    //-------------------------------------------------------
    // Create c3t3_bis
    //-------------------------------------------------------
    std::cout << "Insert 6 points from domain in c3t3_bis, add 1 cell to c3t3_bis\n";
    Polyhedron polyhedron;
    std::ifstream input("data/sphere.off");
    input >> polyhedron;
    input.close();
    Mesh_domain domain(polyhedron);
        
    typedef std::vector<std::pair<Point, Index> > Initial_points_vector;
    Initial_points_vector initial_points;
    domain.construct_initial_points_object()(std::back_inserter(initial_points), 6);
    
    C3t3 c3t3_bis;
    c3t3_bis.insert_surface_points(initial_points.begin(), initial_points.end());
    
    Cell_handle ch_bis = c3t3_bis.triangulation().finite_cells_begin();
    c3t3_bis.add_to_complex(ch_bis,subdomain_index);
    
    std::cout << "\tNumber of cells in c3t3_bis: " << c3t3_bis.number_of_cells() << std::endl;
    std::cout << "\tNumber of facets in c3t3_bis: " << c3t3_bis.number_of_facets() << std::endl;
    std::cout << "\tNumber of vertices in c3t3_bis: " << c3t3_bis.number_of_vertices() << std::endl;
    
    //-------------------------------------------------------
    // Swap c3t3 and c3t3_bis
    //-------------------------------------------------------
    std::cout << "Swap c3t3 and c3t3_bis\n";
    typedef typename C3t3::size_type size_type;

    size_type c3t3_cell_nb = c3t3.number_of_cells(); 
    size_type c3t3_facet_nb = c3t3.number_of_facets();
    size_type c3t3_vertex_nb = c3t3.number_of_vertices();

    size_type c3t3_bis_cell_nb = c3t3_bis.number_of_cells(); 
    size_type c3t3_bis_facet_nb = c3t3_bis.number_of_facets();
    size_type c3t3_bis_vertex_nb = c3t3_bis.number_of_vertices();

    c3t3.swap(c3t3_bis);
    
    std::cout << "\tNumber of cells in c3t3: " << c3t3.number_of_cells() << std::endl;
    std::cout << "\tNumber of facets in c3t3: " << c3t3.number_of_facets() << std::endl;
    std::cout << "\tNumber of vertices in c3t3: " << c3t3.number_of_vertices() << std::endl;
  
    std::cout << "\tNumber of cells in c3t3_bis: " << c3t3_bis.number_of_cells() << std::endl;
    std::cout << "\tNumber of facets in c3t3_bis: " << c3t3_bis.number_of_facets() << std::endl;
    std::cout << "\tNumber of vertices in c3t3_bis: " << c3t3_bis.number_of_vertices() << std::endl;

    assert(c3t3_cell_nb == c3t3_bis.number_of_cells());
    assert(c3t3_facet_nb == c3t3_bis.number_of_facets());
    assert(c3t3_vertex_nb == c3t3_bis.number_of_vertices());

    assert(c3t3_bis_cell_nb == c3t3.number_of_cells());
    assert(c3t3_bis_facet_nb == c3t3.number_of_facets());
    assert(c3t3_bis_vertex_nb == c3t3.number_of_vertices());
  
    //-------------------------------------------------------
    // Modify indices and dimension and verify
    //-------------------------------------------------------
    std::cerr << "Play with indices\n";
    c3t3.add_to_complex(ch, subdomain_index);
    Facet f2 = *c3t3.facets_begin();
    Vertex_handle vh = c3t3.triangulation().vertices_begin();

    c3t3.set_subdomain_index(ch, subdomain_index_bis);
    c3t3.set_surface_index(f2, surface_index_bis);
    c3t3.set_dimension(vh, 1);
    c3t3.set_index(vh, vertex_index);

    assert(c3t3.subdomain_index(ch) == subdomain_index_bis);
    assert(c3t3.surface_index(f2) == surface_index_bis);
    assert(c3t3.in_dimension(vh) == 1);
    assert(c3t3.index(vh) == vertex_index);

    c3t3.set_surface_index(f2.first, f2.second, surface_index);
    assert(c3t3.surface_index(f2) == surface_index);
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
