// Copyright (c) 2020  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_TEST_FACE
#define CGAL_BOOST_GRAPH_TEST_FACE

#include <boost/graph/graph_traits.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

namespace CGAL {



  template <typename VertexRange,typename PMesh>
bool test_face(const VertexRange& vrange, const PMesh& sm)
{
  typedef boost::graph_traits<PMesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<PMesh>::halfedge_descriptor halfedge_descriptor;

  std::vector<typename boost::graph_traits<PMesh>::vertex_descriptor> face(vrange.begin(), vrange.end())

  int N = face.size();
  std::vector<vertex_descriptor> f2(face);
  std::sort(f2.begin(), f2.end());

  std::vector<vertex_descriptor>::iterator it = std::unique(f2.begin(),f2.end());

  if((N > 0) && (it != f2.end())){
    std::cerr << "twice the same vertex" << std::endl;
    return false;
  }

  if(N < 3){
    return false;
  }

  face.push_back(face.front());

  for(int i=0; i < N; i++){
    halfedge_descriptor hd;
    bool found;
    boost::tie(hd,found) = halfedge(face[i],face[i+1],sm);
    if(found && (! is_border(hd,sm))){
      std::cerr << "non-manifold edge" << std::endl;
      return false;
    }
  }

  for(int i=0; i < N; i++){
    if(halfedge(face[i],sm) == boost::graph_traits<PMesh>::null_halfedge()){
      continue;
    }

    if(! is_border(face[i],sm)){
      std::cerr << "non-manifold vertex" << std::endl;
      return false;
    }
  }
  //Test if all halfedges of the new face
    //are possibly consecutive border halfedges in the HDS.
    //Possibly because it may be not directly encoded in the HDS
    //(using next() function ). This situation can occur when one or
    //more facets share only a vertex: For example, the new facet we try to add
    //would make the vertex indices[i] a manifold but this should be forbidden
    //if a facet only incident to that vertex has already been inserted.
    //We check this for each vertex of the sequence.
  for(int i = 0; i < N; ++i) {
      std::size_t prev_index= (i-1+N)%N;
      std::size_t next_index= (i+1)%N;
      vertex_descriptor   previous_vertex = face[ prev_index ];
      vertex_descriptor   next_vertex     = face[ next_index ];

      halfedge_descriptor v = halfedge(face[i],sm);

      if ( v == boost::graph_traits<PMesh>::null_halfedge() ||
           halfedge(previous_vertex,sm) == boost::graph_traits<PMesh>::null_halfedge()||
           halfedge(next_vertex,sm) == boost::graph_traits<PMesh>::null_halfedge()
         ) continue;

      halfedge_descriptor start=v;
      //halfedges pointing to/running out from vertex indices[i]
      //and that need to be possibly consecutive
      halfedge_descriptor prev_hd= boost::graph_traits<PMesh>::null_halfedge(),next_hd= boost::graph_traits<PMesh>::null_halfedge();

      v = opposite(next(v,sm),sm);
      //look for a halfedge incident to vertex indices[i]
      //and which opposite is incident to previous_vertex
      do{
        if(target(opposite(v,sm),sm)==previous_vertex){
          prev_hd=v;
          CGAL_precondition(is_border(prev_hd,sm));
          break;
        }
        v = opposite(next(v,sm),sm);
      }
      while (v!=start);

      if (prev_hd != boost::graph_traits<PMesh>::null_halfedge()){
        v = opposite(next(v,sm),sm);
        //prev_hd and next are already consecutive in the HDS
        if (target(opposite(v,sm),sm)==next_vertex) continue;

        //look for a border halfedge which opposite is
        //incident to next_vertex: set next halfedge
        do
        {
          if (target(opposite(v,sm),sm)==next_vertex){
            next_hd = opposite(v,sm);
            break;
          }
          v= opposite(next(v,sm),sm);
        }
        while(v!=prev_hd);
        if (next_hd==boost::graph_traits<PMesh>::null_halfedge()) continue;

        //check if no constraint prevents
        //prev_hd and next_hd to be adjacent:
        do{
          v= opposite(next(v,sm),sm);
          if ( is_border(opposite(v,sm),sm) ) break;
        }
        while (v!=prev_hd);
        if (v==prev_hd) return false;
        start=v;
      }
    }


  return true;
}


} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_TEST_FACE
