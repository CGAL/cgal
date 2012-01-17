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
// File Description : Test C3T3_with_features class.
//******************************************************************************

#include <CGAL/Bbox_3.h>

#include "test_utilities.h"
#include <CGAL/Mesh_3/Creator_weighted_point_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>

// IO
#include <fstream>
#include <iostream>
#include <CGAL/IO/File_medit.h>



template <typename K>
struct Tester
{
  typedef typename K::FT (Function)(const typename K::Point_3&);
  
  typedef CGAL::Implicit_mesh_domain_3<Function, K>                   Base_domain;
  typedef CGAL::Mesh_domain_with_polyline_features_3<Base_domain>     Md;
  typedef typename CGAL::Mesh_triangulation_3<Md>::type               Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<
    Tr, typename Md::Corner_index, typename Md::Curve_segment_index>  C3t3;

  typedef typename Tr::Geom_traits  Gt;
  typedef typename Gt::FT           FT;
  typedef typename Gt::Point_3      Point;
  typedef CGAL::Mesh_3::Creator_weighted_point_3<FT, Point> Point_creator;

  typedef typename C3t3::Cell_handle    Cell_handle;
  typedef typename C3t3::Facet          Facet;
  typedef typename C3t3::Edge           Edge;
  typedef typename C3t3::Vertex_handle  Vertex_handle;
  
  typedef typename C3t3::Edges_in_complex_iterator      Edge_iterator;
  typedef typename C3t3::Vertices_in_complex_iterator   Vertices_iterator;
  
  typedef typename C3t3::Curve_segment_index  Curve_segment_index;
  typedef typename C3t3::Corner_index         Corner_index;
  typedef typename C3t3::Index                Index;

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
    assert(c3t3.edges_in_complex_begin() == c3t3.edges_in_complex_end());
    assert(c3t3.vertices_in_complex_begin() == c3t3.vertices_in_complex_end());
    assert(c3t3.number_of_cells_in_complex() == 0);
    assert(c3t3.number_of_facets_in_complex() == 0);
    assert(c3t3.number_of_edges_in_complex() == 0);
    assert(c3t3.number_of_vertices_in_complex() == 0);
    
    //-------------------------------------------------------
    // Data generation : fill a triangulation with 4 vertices
    //-------------------------------------------------------
    Point_creator creator;
    Point p1 = creator(0,0,0);
    Point p2 = creator(1,0,0);
    Point p3 = creator(0,1,0);
    Point p4 = creator(0,0,1);

    Vertex_handle vp1 = tr.insert(p1);
    Vertex_handle vp2 = tr.insert(p2);
    Vertex_handle vp3 = tr.insert(p3);
    Vertex_handle vp4 = tr.insert(p4);

    Corner_index corner_index (1);
    Corner_index corner_index_bis (2);
    Curve_segment_index curve_segment_index (1);
    Curve_segment_index curve_segment_index_bis (2);
    Index vertex_index (curve_segment_index);

    //-------------------------------------------------------
    // Add edge to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "\tNumber of edges in c3t3: "
              << c3t3.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3: "
              << c3t3.number_of_vertices_in_complex() << std::endl;
    std::cerr << "Insert one edge in c3t3" << std::endl;

    Edge e = *(tr.finite_edges_begin());
    const Vertex_handle& ev1 = e.first->vertex(e.second);
    const Vertex_handle& ev2 = e.first->vertex(e.third);
    
    c3t3.add_to_complex(e,curve_segment_index);

    std::cerr << "\tNumber of edges in c3t3: "
              << c3t3.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3: "
              << c3t3.number_of_vertices_in_complex() << std::endl;

    assert(e == *(c3t3.edges_in_complex_begin()));
    assert(c3t3.number_of_edges_in_complex() == 1);
    assert(c3t3.number_of_edges_in_complex() == size_type(std::distance(c3t3.edges_in_complex_begin(),
                                                                        c3t3.edges_in_complex_end())));
    assert(c3t3.is_in_complex(e));
    assert(c3t3.is_in_complex(ev1, ev2));
    assert(c3t3.curve_segment_index(e) == curve_segment_index);

    //-------------------------------------------------------
    // Remove cell from c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Remove edge from c3t3" << std::endl;
    
    c3t3.remove_from_complex(ev1, ev2);

    std::cerr << "\tNumber of edges in c3t3: "
              << c3t3.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3: "
              << c3t3.number_of_vertices_in_complex() << std::endl;

    assert(c3t3.number_of_edges_in_complex() == 0);
    assert(! c3t3.is_in_complex(e));
    assert(! c3t3.is_in_complex(ev1, ev2));
    assert(c3t3.curve_segment_index(e) == Curve_segment_index());
    assert(c3t3.curve_segment_index(ev1, ev2) == Curve_segment_index());
    
    //-------------------------------------------------------
    // Add corner to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert one corner in c3t3" << std::endl;

    Vertex_handle v = ++tr.finite_vertices_begin();
    c3t3.add_to_complex(v,corner_index);

    std::cerr << "\tNumber of edges in c3t3: "
              << c3t3.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3: "
              << c3t3.number_of_vertices_in_complex() << std::endl;

    assert(Vertex_handle(c3t3.vertices_in_complex_begin()) == v);
    assert(c3t3.number_of_vertices_in_complex() == 1);
    assert(c3t3.number_of_vertices_in_complex() == size_type(std::distance(c3t3.vertices_in_complex_begin(),
                                                                          c3t3.vertices_in_complex_end())));
    assert(c3t3.is_in_complex(v));
    assert(c3t3.corner_index(v) == corner_index);
    
    //-------------------------------------------------------
    // Remove corner from c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Remove corner from c3t3" << std::endl;

    c3t3.remove_from_complex(v);

    std::cerr << "\tNumber of edges in c3t3: "
              << c3t3.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3: "
              << c3t3.number_of_vertices_in_complex() << std::endl;

    assert(c3t3.vertices_in_complex_begin() == c3t3.vertices_in_complex_begin());
    assert(c3t3.number_of_vertices_in_complex() == 0);
    assert(!c3t3.is_in_complex(v));
    assert(c3t3.corner_index(v) == Corner_index());

    //-------------------------------------------------------
    // Add 1 curve segment (3 edges + 2 corners) to c3t3 and verify
    //-------------------------------------------------------
    std::cerr << "Insert 1 curve segment (3 edges + 2 corners) in c3t3" << std::endl;

    c3t3.add_to_complex(vp1,vp2,curve_segment_index);
    c3t3.add_to_complex(vp2,vp3,curve_segment_index);
    c3t3.add_to_complex(vp3,vp4,curve_segment_index);
    c3t3.add_to_complex(vp1,corner_index);
    c3t3.add_to_complex(vp4,corner_index);
    c3t3.set_dimension(vp1,0);
    c3t3.set_dimension(vp2,1);
    c3t3.set_dimension(vp3,1);
    c3t3.set_dimension(vp4,0);
    
    std::cerr << "\tNumber of edges in c3t3: "
              << c3t3.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3: "
              << c3t3.number_of_vertices_in_complex() << std::endl;

    assert(c3t3.number_of_edges_in_complex() == 3);
    assert(c3t3.number_of_edges_in_complex() == size_type(std::distance(c3t3.edges_in_complex_begin(),
                                                                        c3t3.edges_in_complex_end())));
    assert(c3t3.number_of_vertices_in_complex() == 2);
    assert(c3t3.number_of_vertices_in_complex() == size_type(std::distance(c3t3.vertices_in_complex_begin(),
                                                                          c3t3.vertices_in_complex_end())));

    // -----------------------------------
    // Test iterators
    // The goal here is to test operators and conversion on iterator type
    // -----------------------------------
    typename C3t3::Vertices_in_complex_iterator vit = c3t3.vertices_in_complex_begin();
    v = vit;
    typename C3t3::Triangulation::Vertex& tv1 = *v;
    typename C3t3::Triangulation::Vertex& tv2 = *vit;
    
    assert(   ( v == vp1 && vit->point() == p1 )
           || ( v == vp4 && vit->point() == p4 ) );
    
    assert ( tv1.in_dimension() == tv2.in_dimension() );
    //-------------------------------------------------------
    // Check adjacencies
    //-------------------------------------------------------
    std::vector<std::pair<Vertex_handle,Curve_segment_index> > incident_vertices;
    c3t3.adjacent_vertices_in_complex(vp1,std::back_inserter(incident_vertices));
    
    assert(incident_vertices.size() == 1);
    assert(incident_vertices.front().first == vp2);
    
    incident_vertices.clear();
    c3t3.adjacent_vertices_in_complex(vp3,std::back_inserter(incident_vertices));
    
    assert(incident_vertices.size() == 2);
    assert(   (incident_vertices.front().first == vp2 && incident_vertices.back().first == vp4)
           || (incident_vertices.front().first == vp4 && incident_vertices.back().first == vp2));
    
    //-------------------------------------------------------
    // Create c3t3_bis
    //-------------------------------------------------------
    std::cout << "Insert 6 points in c3t3_bis, add 1 corner and 1 edge to c3t3_bis\n";

    std::vector<Point> points;
    points.push_back(creator(10,11,12));
    points.push_back(creator(11,13,10));
    points.push_back(creator(7,4,6));
    points.push_back(creator(5,2,14));
    points.push_back(creator(1,2,3));
    points.push_back(creator(3,9,13));
    
    C3t3 c3t3_bis;
    c3t3_bis.triangulation().insert(points.begin(),points.end());
    
    Edge e_bis = *(c3t3_bis.triangulation().finite_edges_begin());
    c3t3_bis.add_to_complex(e_bis,curve_segment_index_bis);
    Vertex_handle v_bis = ++c3t3_bis.triangulation().finite_vertices_begin();
    c3t3_bis.add_to_complex(v_bis,corner_index_bis);
    
    std::cerr << "\tNumber of edges in c3t3_bis: "
              << c3t3_bis.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3_bis: "
              << c3t3_bis.number_of_vertices_in_complex() << std::endl;
    std::cout << "\tNumber of vertices in c3t3_bis triangulation: "
              << c3t3_bis.triangulation().number_of_vertices() << std::endl;
    
    //-------------------------------------------------------
    // Swap c3t3 and c3t3_bis
    //-------------------------------------------------------
    std::cout << "Swap c3t3 and c3t3_bis\n";
    typedef typename C3t3::size_type size_type;

    size_type c3t3_edge_nb = c3t3.number_of_edges_in_complex(); 
    size_type c3t3_corner_nb = c3t3.number_of_vertices_in_complex();
    size_type c3t3_vertex_nb = c3t3.triangulation().number_of_vertices();

    size_type c3t3_bis_edge_nb = c3t3_bis.number_of_edges_in_complex(); 
    size_type c3t3_bis_corner_nb = c3t3_bis.number_of_vertices_in_complex();
    size_type c3t3_bis_vertex_nb = c3t3_bis.triangulation().number_of_vertices();

    c3t3.swap(c3t3_bis);
    
    std::cerr << "\tNumber of edges in c3t3: "
              << c3t3.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3: "
              << c3t3.number_of_vertices_in_complex() << std::endl;
    std::cout << "\tNumber of vertices in c3t3: "
              << c3t3.triangulation().number_of_vertices() << std::endl;
  
    std::cerr << "\tNumber of edges in c3t3_bis: "
              << c3t3_bis.number_of_edges_in_complex() << std::endl;
    std::cerr << "\tNumber of corners in c3t3_bis: "
              << c3t3_bis.number_of_vertices_in_complex() << std::endl;
    std::cout << "\tNumber of vertices in c3t3_bis: "
              << c3t3_bis.triangulation().number_of_vertices() << std::endl;

    assert(c3t3_edge_nb == c3t3_bis.number_of_edges_in_complex());
    assert(c3t3_corner_nb == c3t3_bis.number_of_vertices_in_complex());
    assert(c3t3_vertex_nb == c3t3_bis.triangulation().number_of_vertices());

    assert(c3t3_bis_edge_nb == c3t3.number_of_edges_in_complex());
    assert(c3t3_bis_corner_nb == c3t3.number_of_vertices_in_complex());
    assert(c3t3_bis_vertex_nb == c3t3.triangulation().number_of_vertices());
    
    // reset
    c3t3.swap(c3t3_bis);
    
    //-------------------------------------------------------
    // Test edge iterators
    //-------------------------------------------------------
    std::cout << "Test edge iterators\n";
    const Edge& edge_to_modify = *(c3t3.edges_in_complex_begin());
    c3t3.remove_from_complex(edge_to_modify);
    c3t3.add_to_complex(edge_to_modify,curve_segment_index_bis);
    
    typename C3t3::Edges_in_complex_iterator curve_eit =
      c3t3.edges_in_complex_begin(curve_segment_index);
    typename C3t3::Edges_in_complex_iterator curve_eit_bis =
      c3t3.edges_in_complex_begin(curve_segment_index_bis);
    typename C3t3::Edges_in_complex_iterator eend =
      c3t3.edges_in_complex_end();
    
    std::cout << "\tNumber of edges of index '" << curve_segment_index << "': "
              << std::distance(curve_eit,eend) << std::endl;
    std::cout << "\tNumber of edges of index '" << curve_segment_index_bis << "': "
              << std::distance(curve_eit_bis,eend) << std::endl;
    
    assert ( std::distance(curve_eit,eend) == 2 );
    assert ( std::distance(curve_eit_bis,eend) == 1 );
    assert ( c3t3.curve_segment_index(*curve_eit) == curve_segment_index );
    assert ( c3t3.curve_segment_index(*curve_eit_bis) == curve_segment_index_bis );
    
    //-------------------------------------------------------
    // Test vertex iterators
    //-------------------------------------------------------
    std::cout << "Test vertex iterators\n";
    const Vertex_handle& vertex_to_modify = c3t3.vertices_in_complex_begin();
    c3t3.remove_from_complex(vertex_to_modify);
    c3t3.add_to_complex(vertex_to_modify,corner_index_bis);
    
    typename C3t3::Vertices_in_complex_iterator corner_vit =
      c3t3.vertices_in_complex_begin(corner_index);
    typename C3t3::Vertices_in_complex_iterator corner_vit_bis =
      c3t3.vertices_in_complex_begin(corner_index_bis);
    typename C3t3::Vertices_in_complex_iterator vend =
      c3t3.vertices_in_complex_end();
    
    std::cout << "\tNumber of vertices of index '" << corner_index << "': "
              << std::distance(corner_vit,vend) << std::endl;
    std::cout << "\tNumber of vertices of index '" << corner_index_bis << "': "
              << std::distance(corner_vit_bis,vend) << std::endl;
    
    assert ( std::distance(corner_vit,vend) == 1 );
    assert ( std::distance(corner_vit_bis,vend) == 1 );
    assert ( c3t3.corner_index(corner_vit) == corner_index );
    assert ( c3t3.corner_index(corner_vit_bis) == corner_index_bis );
  }
};


int main()
{
  Tester<K_e_i> test_epic;
  test_epic();

  return EXIT_SUCCESS;
}
