// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory (France).


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
// Author(s) : Maxime Gimeno, Pierre Alliez
#ifndef CGAL_OUTPUT_SURFACE_FACETS_TO_FACEGRAPH_H
#define CGAL_OUTPUT_SURFACE_FACETS_TO_FACEGRAPH_H

#include <CGAL/license/Surface_mesher.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <map>

//! Gets reconstructed surface out of a SurfaceMeshComplex_2InTriangulation_3 object.
//!
//! This variant exports the surface as a FaceGraph and appends it to `graph`.
//! It requires the surface to be manifold. For this purpose,
//! you may call make_surface_mesh() with Manifold_tag or Manifold_with_boundary_tag parameter.
//!
//! @tparam C2T3 model of the SurfaceMeshComplex_2InTriangulation_3 concept.
//! @tparam FaceGraph a model of FaceGraph.
//!
//! @param c2t3 an instance of a manifold C2T3.
//! @param graph an instance of FaceGraph.
template<class C2T3, class FaceGraph>
void c2t3_to_facegraph(const C2T3& c2t3, FaceGraph& graph)
{
  typedef typename boost::property_map<FaceGraph, boost::vertex_point_t>::type VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename C2T3::Triangulation Tr;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Vertex_iterator Vertex_iterator;
  typedef typename Tr::Geom_traits::Vector_3 Vector;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;

  const Tr& tr = c2t3.triangulation();
  VertexPointMap vpmap = get(boost::vertex_point, graph);
  const typename Tr::size_type number_of_facets = c2t3.number_of_facets();
  {
    // Finite vertices coordinates.
    Finite_facets_iterator fit = tr.finite_facets_begin();
    std::set<Facet> oriented_set;
    std::stack<Facet> stack;

    CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

        while (oriented_set.size() != number_of_facets) {
      while ( fit->first->is_facet_on_surface(fit->second) == false ||
              oriented_set.find(*fit) != oriented_set.end() ||

              oriented_set.find(c2t3.opposite_facet(*fit)) !=
              oriented_set.end() ) {
        ++fit;
      }
      oriented_set.insert(*fit);
      stack.push(*fit);
      while(! stack.empty() ) {
        Facet f = stack.top();
        stack.pop();
        for(int ih = 0 ; ih < 3 ; ++ih) {
          const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
          const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));
          if( c2t3.face_status(Edge(f.first, i1, i2)) == C2T3::REGULAR ) {
            Facet fn = c2t3.neighbor(f, ih);
            if (oriented_set.find(fn) == oriented_set.end() &&
                oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
            {
              oriented_set.insert(fn);
              stack.push(fn);
            }
          } // end "if the edge is regular"
        } // end "for each neighbor of f"
      } // end "stack non empty"
    } // end "oriented_set not full"

    // Orients the whole mesh towards outside:
    // - find the facet with max z
    typename std::set<Facet>::const_iterator top_facet = oriented_set.begin();
    for(typename std::set<Facet>::const_iterator fit = oriented_set.begin();
        fit != oriented_set.end();
        ++fit)
    {
      double top_z =
          (top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 0))->point().z()
           + top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 1))->point().z()
           + top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 2))->point().z())/3.;
      double z =
          (fit->first->vertex(tr.vertex_triple_index(fit->second, 0))->point().z()
           + fit->first->vertex(tr.vertex_triple_index(fit->second, 1))->point().z()
           + fit->first->vertex(tr.vertex_triple_index(fit->second, 2))->point().z())/3.;
      if (top_z < z)
        top_facet = fit;
    }
    // - orient the facet with max z towards +Z axis
    Vertex_handle v0 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 0));
    Vertex_handle v1 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 1));
    Vertex_handle v2 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 2));
    Vector normal = cross_product(v1->point()-v0->point(), v2->point()-v1->point());
    const Vector Z(0, 0, 1);
    bool regular_orientation = (Z * normal >= 0);

    //used to set indices of vertices
    std::map<Vertex_handle, int> V;
    int inum = 0;


    //add vertices
    std::vector<typename boost::graph_traits<FaceGraph>::vertex_descriptor> vertices;
    for(Vertex_iterator vit = tr.vertices_begin();
        vit != tr.vertices_end();
        ++vit)
    {

      typename boost::graph_traits<FaceGraph>::vertex_descriptor v = add_vertex(graph);
      vertices.push_back(v);
      put(vpmap,
          v,
           Point_3(
            vit->point().x(),
            vit->point().y(),
            vit->point().z())
          );
      V.insert(std::make_pair(vit, inum++));
    }
    //add faces
    for(typename std::set<Facet>::const_iterator fit =
        oriented_set.begin();
        fit != oriented_set.end();
        ++fit)
    {
      int id0(V[fit->first->vertex(tr.vertex_triple_index(fit->second, 0))]),
          id1(regular_orientation ? V[fit->first->vertex(tr.vertex_triple_index(fit->second, 1))]
              : V[fit->first->vertex(tr.vertex_triple_index(fit->second, 2))]),
          id2(regular_orientation ? V[fit->first->vertex(tr.vertex_triple_index(fit->second, 2))]
              : V[fit->first->vertex(tr.vertex_triple_index(fit->second, 1))]);
      std::vector<typename boost::graph_traits<FaceGraph>::vertex_descriptor> face;
      face.resize(3);
      face[0] = vertices[id0];
      face[1] = vertices[id1];
      face[2] = vertices[id2];
      CGAL::Euler::add_face(face, graph);
      CGAL_assertion_code(++nb_facets);
    }
    CGAL_assertion(nb_facets == number_of_facets);
  }
}
#endif // CGAL_OUTPUT_SURFACE_FACETS_TO_FACEGRAPH_H
