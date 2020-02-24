// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Frank Da, David Cohen-Steiner, Andreas Fabri

#ifndef CGAL_AFSR_ORIENT_H
#define CGAL_AFSR_ORIENT_H

#include <CGAL/license/Advancing_front_surface_reconstruction.h>


namespace CGAL {
  namespace AFSR {


    template <class Triangulation, class TDS, class Filter>
    typename TDS::Vertex_handle
    orient(TDS& tds, const Advancing_front_surface_reconstruction<Triangulation,Filter>& surface)
    {

      typedef typename TDS::Vertex_handle Vertex_handle;
      typedef std::pair<Vertex_handle,Vertex_handle> Vh_pair;
      typedef typename TDS::Face_handle Face_handle;
      typedef typename TDS::Edge Edge;

      Triangulation& T = surface.triangulation_3();
      // create an infinite-vertex and  infinite faces with the
      // boundary edges if any.
      // return the infinite vertex if created
      Vertex_handle vinf;

      std::vector<Vertex_handle > vvh;
      if(tds.number_of_vertices() != 0)    tds.clear();
      int dim = 2;
      tds.set_dimension(dim);

      CGAL::Unique_hash_map<typename Triangulation::Vertex_handle, int> vertex_index_map(-1, T.number_of_vertices());

      int i=0;
      for (typename Triangulation::Finite_vertices_iterator v_it = T.finite_vertices_begin();
           v_it != T.finite_vertices_end();
           v_it++){
        typename CGAL::Unique_hash_map<typename Triangulation::Vertex_handle, int>::Data& d = vertex_index_map[v_it];
        if ((!v_it->is_exterior()) && d == -1){
          d = i;
          Vertex_handle vh = tds.create_vertex();
          vvh.push_back(vh);
          vh->set_point(v_it->point());
          i++;
        }
      }
      std::map<Vh_pair, Edge> edge_map;


      for(typename Triangulation::Finite_facets_iterator f_it = T.finite_facets_begin();
          f_it != T.finite_facets_end();
          f_it++)
        {
          typename Triangulation::Cell_handle n, c = (*f_it).first;
          int ni, ci = (*f_it).second;
          n = c->neighbor(ci);
          ni = n->index(c);
          int i1, i2 ,i3;

          if (c->is_selected_facet(ci))
            {
              i1 = (ci+1) & 3;
              i2 = (ci+2) & 3;
              i3 = (ci+3) & 3;

              Face_handle fh = tds.create_face(vvh[vertex_index_map[c->vertex(i1)]],
                                               vvh[vertex_index_map[c->vertex(i2)]],
                                               vvh[vertex_index_map[c->vertex(i3)]]);
              vvh[vertex_index_map[c->vertex(i1)]]->set_face(fh);
              vvh[vertex_index_map[c->vertex(i2)]]->set_face(fh);
              vvh[vertex_index_map[c->vertex(i3)]]->set_face(fh);
              for (int ih  = 0; ih < 3; ++ih) {
                tds.set_adjacency(fh, ih, edge_map);
              }
            }

          if (n->is_selected_facet(ni))
            {
              i1 = (ni+1) & 3;
              i2 = (ni+2) & 3;
              i3 = (ni+3) & 3;
              Face_handle fh = tds.create_face(vvh[vertex_index_map[n->vertex(i1)]],
                                               vvh[vertex_index_map[n->vertex(i2)]],
                                               vvh[vertex_index_map[n->vertex(i3)]]);
              vvh[vertex_index_map[n->vertex(i1)]]->set_face(fh);
              vvh[vertex_index_map[n->vertex(i2)]]->set_face(fh);
              vvh[vertex_index_map[n->vertex(i3)]]->set_face(fh);
              for (int ih  = 0; ih < 3; ++ih) {
                tds.set_adjacency(fh, ih, edge_map);
              }
            }

        }

      if ( !edge_map.empty()) {
        vinf = tds.create_vertex();
        std::map<Vh_pair, Edge> inf_edge_map;
        while (!edge_map.empty()) {
          Face_handle fh = edge_map.begin()->second.first;
          int ih = edge_map.begin()->second.second;
          Face_handle fn = tds.create_face( vinf,
                                            fh->vertex(TDS::cw(ih)),
                                            fh->vertex(TDS::ccw(ih)));
          vinf->set_face(fn);
          tds.set_adjacency(fn, 0, fh, ih);
          tds.set_adjacency(fn, 1, inf_edge_map);
          tds.set_adjacency(fn, 2, inf_edge_map);
          edge_map.erase(edge_map.begin());
        }
        CGAL_assertion(inf_edge_map.empty());
      }


      tds.reorient_faces();
      return vinf;

    }



  } // namespace AFSR
} // namespace CGAL

#endif //CGAL_AFSR_ORIENT_H
