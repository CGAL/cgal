// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

// file: nearest_neighbor_graph.cpp

#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <iterator>

#include <CGAL/Simple_cartesian.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

#if 0
// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

struct Gt
  : public CGAL::Segment_Delaunay_graph_filtered_traits_2<Rep> {};

typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt>  SDG2;


typedef Gt::Point_2           Point_2;
typedef SDG2::Vertex_handle   Vertex_handle;

typedef CGAL::Directed_graph<Vertex_handle>   Nearest_neighbor_graph;
typedef Nearest_neighbor_graph::Node_handle   Node_handle;
typedef Nearest_neighbor_graph                NNG;
typedef std::map<Vertex_handle,Node_handle>   V2N;
#endif

#include <CGAL/Filtered_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include "Nearest_neighbor_graph_2.h"
#include "Delaunay_triangulation_nearest_neighbor_graph_traits_2.h"

typedef CGAL::Filtered_kernel<Rep>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt>  DG;

typedef DG::Geom_traits::Point_2 Point_2;

typedef CGAL::Delaunay_triangulation_nearest_neighbor_graph_traits_2<DG>
NNG_Traits;

typedef CGAL::Nearest_neighbor_graph_2<DG,NNG_Traits>  NNG;


int main()
{
  DG   dg;

#if 0
#if 0
  Point_2 p1(1.0,1.0);
  Point_2 p2(5.0,7.0);
  Point_2 p3(4.0,9.0);
  Point_2 p4(2.0,4.0);
#else
  Point_2 p1(0.0,0.0);
  Point_2 p2(0.0,1.0);
  Point_2 p3(2.0,0.0);
  Point_2 p4(2.0,1.0);
#endif
#else
#if 0
  Point_2 p1(0.0,0.0);
  Point_2 p2(1.0,0.0);
  Point_2 p3(0.0,1.0);
  Point_2 p4(1.0,1.0);
#else
  Point_2 p1(0.0,0.0);
  Point_2 p2(1.0,0.0);
  Point_2 p3(2.0,0.0);
  Point_2 p4(4.0,0.0);
#endif
#endif

  dg.insert(p1);
  dg.insert(p2);
  dg.insert(p3);
  dg.insert(p4);

  // validate the segment Delaunay graph
  assert( dg.is_valid() );

  NNG  nng(dg);

  NNG::Traits::Access_site_2 accessor = nng.traits().access_site_2_object();

  // print out nearest neighbor graph
  for (NNG::Node_iterator nit = nng.nodes_begin();
       nit != nng.nodes_end(); ++nit) {
    std::cout << "Node: " << accessor( nng[nit] ) << std::endl;
    std::cout << "\tIn-neighbors (in-degree: " << nng.in_degree(nit)
              << "): " << std::endl;
    for (NNG::Node_iterator nnit = nng.in_neighbors_begin(nit);
         nnit != nng.in_neighbors_end(nit); ++nnit) {
      std::cout << "\t\t" << accessor( nng[nnit] ) << std::endl;
    }
    std::cout << "\tOut-neighbors (out-degree: " << nng.out_degree(nit)
              << "):" << std::endl;
    for (NNG::Node_iterator nnit = nng.out_neighbors_begin(nit);
         nnit != nng.out_neighbors_end(nit); ++nnit) {
      std::cout << "\t\t" << accessor( nng[nnit] ) << std::endl;
    }
  }


  return 0;
}
