// Copyright (c) 2015-2021  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno
//                 Mariette Yvinec

#ifndef CGAL_BGL_IO_TDS_2_OFF_H
#define CGAL_BGL_IO_TDS_2_OFF_H
#include <CGAL/IO/OFF.h>
#include <CGAL/Triangulation_data_structure_2.h>


namespace CGAL {
namespace IO {

template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
off_file_input( std::istream& is, Triangulation_data_structure_2<Vb,Fb>& tds, bool verbose = false)
{
  typedef typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle Vertex_handle;
  typedef typename Triangulation_data_structure_2<Vb,Fb>::Face_iterator Face_handle;
  typedef std::pair<Vertex_handle,Vertex_handle> Vh_pair;
  typedef std::pair<Face_handle, int>                Edge;
  // input from an OFF file
  // assume a dimension 2 triangulation
  // create an infinite-vertex and  infinite faces with the
  // boundary edges if any.
  // return the infinite vertex if created
  Vertex_handle vinf;
  File_scanner_OFF scanner(is, verbose);
  if (! is) {
    if (scanner.verbose()) {
         std::cerr << " " << std::endl;
         std::cerr << "TDS::off_file_input" << std::endl;
         std::cerr << " input error: file format is not OFF." << std::endl;
    }
    return vinf;
  }

  if(tds.number_of_vertices() != 0)    tds.clear();
  int dim = 2;
  tds.set_dimension(dim);

  std::vector<Vertex_handle > vvh(scanner.size_of_vertices());
  std::map<Vh_pair, Edge> edge_map;
  typedef typename Vb::Point   Point;

  // read vertices
  std::size_t i;
  for ( i = 0; i < scanner.size_of_vertices(); i++) {
    Point p;
    file_scan_vertex( scanner, p);
    vvh[i] = tds.create_vertex();
    vvh[i]->set_point(p);
    scanner.skip_to_next_vertex( i);
  }
  if ( ! is ) {
    is.clear( std::ios::badbit);
    return vinf;
  }
  //vinf = vvh[0];

  // create the facets
  for ( i = 0; i < scanner.size_of_facets(); i++) {
    Face_handle fh = tds.create_face();
    std::size_t no;
    scanner.scan_facet( no, i);
    if( ! is || no != 3) {
      if ( scanner.verbose()) {
        std::cerr << " " << std::endl;
        std::cerr << "TDS::off_file_input" << std::endl;
        std::cerr << "facet " << i << "does not have  3 vertices."
                  << std::endl;
      }
      is.clear( std::ios::badbit);
      return vinf;
    }

    for ( std::size_t j = 0; j < no; ++j) {
      std::size_t index;
      scanner.scan_facet_vertex_index( index, j+1, i);
      fh->set_vertex(j, vvh[index]);
      vvh[index]->set_face(fh);
    }

    for (std::size_t ih  = 0; ih < no; ++ih) {
        tds.set_adjacency(fh, ih, edge_map);
    }
  }

  // deal with  boundaries
  if ( !edge_map.empty()) {
    vinf = tds.create_vertex();
    std::map<Vh_pair, Edge> inf_edge_map;
   while (!edge_map.empty()) {
     Face_handle fh = edge_map.begin()->second.first;
     int ih = edge_map.begin()->second.second;
     Face_handle fn = tds.create_face( vinf,
                                   fh->vertex(tds.cw(ih)),
                                   fh->vertex(tds.ccw(ih)));
     vinf->set_face(fn);
     tds.set_adjacency(fn, 0, fh, ih);
     tds.set_adjacency(fn, 1, inf_edge_map);
     tds.set_adjacency(fn, 2, inf_edge_map);
     edge_map.erase(edge_map.begin());
   }
   CGAL_triangulation_assertion(inf_edge_map.empty());
  }


  // coherent orientation
  tds.reorient_faces();
  return vinf;
}

}
}
#endif // CGAL_BGL_IO_TDS_2_OFF_H
