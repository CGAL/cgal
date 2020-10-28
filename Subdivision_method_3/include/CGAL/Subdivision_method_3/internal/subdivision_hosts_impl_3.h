// Copyright (c) 2005-2017 GeometryFactory (France).  All Rights Reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//

#ifndef CGAL_SUBDIVISION_HOSTS_IMPL_3_H
#define CGAL_SUBDIVISION_HOSTS_IMPL_3_H

#include <CGAL/Subdivision_method_3/internal/Euler_extensions.h>

#include <CGAL/basic.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/circulator.h>
#include <CGAL/tags.h>

#include <boost/unordered_map.hpp>

#include <iterator>
#include <list>
#include <vector>

namespace CGAL {

namespace Subdivision_method_3 {

namespace internal {

template <class PolygonMesh>
void call_reserve(PolygonMesh& pm, std::size_t v, std::size_t e, std::size_t f)
{
  typedef boost::graph_traits<PolygonMesh> GT;
  reserve(pm,
          static_cast<typename GT::vertices_size_type>(v),
          static_cast<typename GT::edges_size_type>(e),
          static_cast<typename GT::faces_size_type>(f));
}

template <class Poly, class VertexPointMap, class Mask>
void PQQ_1step(Poly& p, VertexPointMap vpm, Mask mask) {
  typedef typename boost::graph_traits<Poly>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;
  typedef typename boost::graph_traits<Poly>::face_descriptor         face_descriptor;

  typedef Halfedge_around_face_circulator<Poly>  Halfedge_around_facet_circulator;

  //First back up initial vertices/faces/edges
  std::vector<vertex_descriptor> p_vertices(vertices(p).first, vertices(p).second);
  std::vector<face_descriptor> p_faces(faces(p).first, faces(p).second);
  std::vector<edge_descriptor> p_edges(edges(p).first, edges(p).second);

  std::size_t num_v = p_vertices.size();
  std::size_t num_e = p_edges.size();
  std::size_t num_f = p_faces.size();

  // Build a new vertices buffer has the following structure
  //
  // 0 1 ... e_begin ... f_begin ... (end_of_buffer)
  // 0 ... e_begin-1       : store the positions of the vertex-vertices
  // e_begin ... f_begin-1 : store the positions of the edge-vertices
  // f_begin ... (end)     : store the positions of the face-vertices
  // The index of the vertices buffer should 1-1 map to the distance
  // of the corresponding iterator to the begin of the iterator.

  // We need to reserve the memory to prevent reallocation.
  call_reserve(p,num_v+num_e+num_f, 4*2*num_e, 4*num_e/2);

  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  Point* vertex_point_buffer = new Point[num_v + num_e + num_f];
  Point* edge_point_buffer = vertex_point_buffer + num_v;
  Point* face_point_buffer = edge_point_buffer + num_e;

  int i=0;
  boost::unordered_map<vertex_descriptor,int> v_index;
  for(vertex_descriptor vh : p_vertices){
    v_index[vh]= i++;
  }

  std::vector<bool> v_onborder(num_v);
  typename std::vector<face_descriptor>::iterator fitr=p_faces.begin();
  for (std::size_t i = 0; i < num_f; i++, ++fitr)
    mask.face_node(*fitr, face_point_buffer[i]);
  {
    std::size_t i = 0;
    for(edge_descriptor ed : p_edges){
      if(is_border(ed,p)){
        halfedge_descriptor h=halfedge(ed,p);
        if (is_border(h,p)) h=opposite(h,p);
        int v = v_index[target(h,p)];
        v_onborder[v] = true;
        mask.border_node(h, edge_point_buffer[i], vertex_point_buffer[v]);

      }else{
        mask.edge_node(halfedge(ed,p), edge_point_buffer[i]);
      }
      ++i;
    }
  }

  typename std::vector<vertex_descriptor>::iterator vitr = p_vertices.begin();
  for (std::size_t i = 0; i < num_v; i++, ++vitr)
    if (!v_onborder[v_index[*vitr]]) mask.vertex_node(*vitr, vertex_point_buffer[i]);

  // Build the connectivity using insert_vertex() and insert_edge()
  // 1. insert_vertex() to all edges and set them to new positions
  // 2. insert_edge() between 2 randomly selected neighboring new inserted
  //    vertices
  // 3. insert_vertex() to the new inserted edge and set them to new positions
  // 4. insert_edge() between all other new inserted vertices of step 1 and
  //    the new inserted vertex of step 3
  // Step 1.
  typename std::vector<edge_descriptor>::iterator eitr = p_edges.begin();
  for (std::size_t i = 0; i < num_e; i++, ++eitr) {
    vertex_descriptor vh = insert_vertex(p, halfedge(*eitr,p));
    put(vpm, vh, edge_point_buffer[i]);
  }

  // TODO: the topoloy modification can be done by a template function
  //       and that gives the user a chance to create new topological masks.
  fitr = p_faces.begin();
  for (std::size_t i = 0; i < num_f; i++, ++fitr) {
    // Step 2.
    Halfedge_around_facet_circulator hcir_begin(halfedge(*fitr,p),p);
    Halfedge_around_facet_circulator hcir = hcir_begin;

    halfedge_descriptor e1 = * ++hcir; // e1 points to the newly inserted vertex
    ++hcir; // Skips one original vertex
    halfedge_descriptor e2 = * ++hcir; // points to the next newly inserted vertex
    ++hcir; // Must move the cir before inserts the new edge !!
    halfedge_descriptor newe = insert_edge(p, e1, e2);

    // Step 3.
    halfedge_descriptor newv = insert_vertex_return_edge(p, newe);
    newv = prev(opposite(newv,p),p); // change newv to the larger face and
    // still points to the newly inserted
    // vertex
    // Update the geometry data of the newly inserted face-vertices
    put(vpm, target(newv,p), face_point_buffer[i]);

    // Step 4.
    while (hcir != hcir_begin) {
      e1 = * ++hcir;
      ++hcir; // Must move the cir before inserts the new edge !!
      insert_edge(p, e1, newv);
    }
  }

  // Update the geometry data of the newly inserted vertices by the
  // vertices buffer
  vitr = p_vertices.begin();
  for (std::size_t i = 0; i < num_v; i++, ++vitr)
    put(vpm, *vitr, vertex_point_buffer[i]);

  CGAL_postcondition(CGAL::is_valid_polygon_mesh(p));
  delete []vertex_point_buffer;
}

// ======================================================================
template <class Poly,class VertexPointMap, class Mask>
void PTQ_1step(Poly& p, VertexPointMap vpm, Mask mask) {
  typedef typename boost::graph_traits<Poly>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;
  typedef typename boost::graph_traits<Poly>::face_descriptor         face_descriptor;

  typedef Halfedge_around_face_circulator<Poly>  Halfedge_around_face_circulator;

  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  //First back up initial vertices/faces/edges
  std::vector<vertex_descriptor> p_vertices(vertices(p).first, vertices(p).second);
  std::vector<face_descriptor> p_faces(faces(p).first, faces(p).second);
  std::vector<edge_descriptor> p_edges(edges(p).first, edges(p).second);

  std::size_t num_v = p_vertices.size();
  std::size_t num_e = p_edges.size();
  std::size_t num_f = p_faces.size();

  // Build a new vertices buffer has the following structure
  //
  // 0 1 ... e_begin ... f_begin ... (end_of_buffer)
  // 0 ... e_begin-1       : store the positions of the vertex-vertices
  // e_begin ... (end)     : store the positions of the edge-vertices
  // The index of the vertices buffer should 1-1 map to the distance
  // of the corresponding iterator to the begin of the iterator.

  // We need to reserve the memory to prevent reallocation.
  call_reserve(p,num_v + num_e, 2*2*num_e, 4*num_e/2);

  Point* vertex_point_buffer = new Point[num_v + num_e];
  Point* edge_point_buffer = vertex_point_buffer + num_v;

  int i=0;
  boost::unordered_map<vertex_descriptor,int> v_index;
  for(vertex_descriptor vh : p_vertices){
    v_index[vh]= i++;
  }
  std::vector<bool> v_onborder(num_v);

  {
    std::size_t i = 0;
    for(edge_descriptor ed : p_edges){
      if(! is_border(ed,p)){
        mask.edge_node(halfedge(ed,p), edge_point_buffer[i]);
      } else{
        halfedge_descriptor h = halfedge(ed,p);
        if (is_border(h, p)) h = opposite(h,p);
        int v = v_index[target(h,p)];
        v_onborder[v] = true;
        mask.border_node(h, edge_point_buffer[i], vertex_point_buffer[v]);
      }
      ++i;
    }
  }
  typename std::vector<vertex_descriptor>::iterator vitr = p_vertices.begin();
  for (std::size_t i = 0; i < num_v; i++, ++vitr)
    if (!v_onborder[i]) mask.vertex_node(*vitr, vertex_point_buffer[i]);

  // Build the connectivity using insert_vertex() and insert_edge()
  // 1. insert_vertex() to all edges and set them to new positions
  // 2. insert_edge() between 2 randomly selected neighboring new inserted
  //    vertices
  // 3. insert_vertex() to the new inserted edge and set them to new positions
  // 4. insert_edge() between all other new inserted vertices of step 1 and
  //    the new inserted vertex of step 3
  // Step 1.
  typename std::vector<edge_descriptor>::iterator eit = p_edges.begin();
  for (std::size_t i = 0; i < num_e; i++, ++eit) {
    vertex_descriptor vh = insert_vertex(p, halfedge(*eit,p));
    put(vpm,vh, edge_point_buffer[i]);
  }

  typename std::vector<face_descriptor>::iterator fitr = p_faces.begin();
  for (std::size_t i = 0; i < num_f; i++, ++fitr) {
    // Step 2.
    Halfedge_around_face_circulator hcir_begin(halfedge(*fitr,p),p);
    Halfedge_around_face_circulator hcir = hcir_begin;

    // After linsub, the facet valence = 6
    CGAL_assertion(circulator_size(hcir)==6);

    halfedge_descriptor e1 = *(++hcir);
    ++hcir;
    halfedge_descriptor e2 = *(++hcir);
    ++hcir;
    halfedge_descriptor e3 = *(++hcir);
    e2 = insert_edge(p, e1, e2);
    e3 = insert_edge(p, e2, e3);
    insert_edge(p, e3, e1);
  }

  // Update the geometry data of the newly inserted vertices by the
  // vertices buffer
  vitr = p_vertices.begin();
  for (std::size_t i = 0; i < num_v; i++, ++vitr)
    put(vpm, *vitr, vertex_point_buffer[i]);

  CGAL_postcondition(CGAL::is_valid_polygon_mesh(p));
  delete []vertex_point_buffer;
}


// ======================================================================
template <class Poly, class VertexPointMap, class Mask>
void DQQ_1step(Poly& p, VertexPointMap vpm, Mask mask) {
  typedef typename boost::graph_traits<Poly>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;
  typedef typename boost::graph_traits<Poly>::face_descriptor         face_descriptor;


  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  //First back up initial vertices/faces/edges
  std::vector<vertex_descriptor> p_vertices(vertices(p).first, vertices(p).second);
  std::vector<face_descriptor> p_faces(faces(p).first, faces(p).second);
  std::vector<edge_descriptor> p_edges(edges(p).first, edges(p).second);

  std::size_t num_v = p_vertices.size();
  std::size_t num_e = p_edges.size();
  std::size_t num_f = p_faces.size();

  std::vector<halfedge_descriptor> border_halfedges;
  for(edge_descriptor ed : p_edges){
    if(is_border(ed,p)){
      border_halfedges.push_back(halfedge(ed,p));
    }
  }
  Point* point_buffer = new Point[num_e*2];

  // Build the point_buffer
  int pi = 0;
  int vi=0;
  std::vector<bool> is_border_vertex(num_v, false);
  for(vertex_descriptor vd : p_vertices){
    for(halfedge_descriptor hd : halfedges_around_target(vd,p)){
      if (! is_border(hd,p)){
        mask.corner_node(hd, point_buffer[pi++]);
      }
      else
        is_border_vertex[vi]=true;
    }
    ++vi;
  }

  // Reserve to avoid rellocations during insertions
  call_reserve(p,num_v+num_e+num_f, 2*num_e, (2+4+2)*num_e);

  // Build the connectivity using insert_vertex() and insert_edge()
  pi = 0;
  typename std::vector<vertex_descriptor>::iterator vitr=p_vertices.begin();
  for (std::size_t i = 0; i < num_v; ++i) {
    vertex_descriptor vh = *vitr;
    ++vitr;

    Halfedge_around_target_circulator<Poly> vcir(vh,p);
    typename boost::graph_traits<Poly>::degree_size_type vn = degree(vh,p);
    for (typename boost::graph_traits<Poly>::degree_size_type j = 0; j < vn; ++j) {
      halfedge_descriptor e = *vcir;
      ++vcir;
      if (! is_border(e,p)) {
        vertex_descriptor v = insert_vertex(p, e);
        put(vpm, v, point_buffer[pi++]);
      }
    }

    vcir = Halfedge_around_target_circulator<Poly>(vh,p);
    for (typename boost::graph_traits<Poly>::vertices_size_type j = 0; j < vn; ++j) {
      if (! is_border(*vcir,p)) {
        halfedge_descriptor e1 = prev(*vcir, p);
        ++vcir;
        if (! is_border(*vcir,p)) {
          halfedge_descriptor e2 = opposite(*vcir,p);
          insert_edge(p, e1, e2);
        }
      } else {
        ++vcir;
      }
    }
  }

  typename std::vector<edge_descriptor>::iterator eitr = p_edges.begin();
  for (typename boost::graph_traits<Poly>::edges_size_type i = 0; i < num_e; ++i) {
    halfedge_descriptor eh = halfedge(*eitr,p);
    ++eitr;
    if (! is_border(edge(eh,p),p)) {
      insert_edge(p, prev(prev(eh,p),p), eh);
      eh = opposite(eh,p);
      insert_edge(p, prev(prev(eh,p),p), eh);
      Euler::join_face(eh,p);
    } else {
      if (is_border(eh,p)) {
        eh = opposite(eh,p);
        insert_edge(p, eh, prev(prev(eh,p),p));
      } else
        insert_edge(p, prev(prev(eh,p),p), eh);
    }
  }


  for(halfedge_descriptor eeh : border_halfedges){
    halfedge_descriptor eh = eeh;
    if (is_border(eh,p)){
      eh = opposite(eh,p);
    }

    halfedge_descriptor ehe = eh;
    eh = opposite(prev(eh,p),p);
    while (! is_border(eh,p)) {
      Euler::remove_face(ehe,p);
      ehe = eh;
      eh = opposite(prev(eh,p),p);
    }
    Euler::remove_face(ehe,p);
  }

  vitr = p_vertices.begin();
  for (typename boost::graph_traits<Poly>::vertices_size_type i = 0; i < num_v; ++i) {
    vertex_descriptor vh = *vitr;
    ++vitr;
    if (!is_border_vertex[i])
      Euler::remove_center_vertex(halfedge(vh,p),p);
  }

  delete []point_buffer;
}

// ======================================================================
template <class Poly, class VertexPointMap, class Mask>
void Sqrt3_1step(Poly& p, VertexPointMap vpm, Mask mask,
                 const bool refine_border = false) {
  // `refine_border` is a boolean that is meant to be true only every SECOND step
  // of the subdivision. In particular, this function makes uses of the fact
  // that there is at most a single border edge in a face, which is true if
  // the mesh is obtained from a sqrt3 subdivision before, but might otherwise
  // be wrong.

  typedef typename boost::graph_traits<Poly>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;
  typedef typename boost::graph_traits<Poly>::face_descriptor         face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  //First back up initial vertices/faces/edges
  std::vector<vertex_descriptor> p_vertices(vertices(p).first, vertices(p).second);
  std::vector<face_descriptor> p_faces(faces(p).first, faces(p).second);
  std::vector<edge_descriptor> p_edges(edges(p).first, edges(p).second);

  std::size_t num_v = p_vertices.size();
  std::size_t num_e = p_edges.size();
  std::size_t num_f = p_faces.size();

  // reserve enough size for the new points
  std::size_t new_pts_size = num_f;
  if(refine_border) {
    for(edge_descriptor ed : p_edges){
      if(is_border(ed, p))
        ++new_pts_size;
    }
  }
  Point* cpt = new Point[new_pts_size];

  // size of the subdivided mesh
  call_reserve(p,num_v + new_pts_size, (num_e + 2*num_f + new_pts_size)*2, 3*num_f);

  // keep in memory whether a face is incident to the border and, if so, which
  // halfedge corresponds to THE (there can only be one) border edge.
  std::vector<halfedge_descriptor> face_halfedge_border(num_f,
                                                        boost::graph_traits<Poly>::null_halfedge());
  std::list<std::pair<vertex_descriptor, Point> > new_positions;

  // compute the positions of new points
  std::size_t i = 0;
  std::size_t face_id = 0;
  for(face_descriptor fd : p_faces) {
    //ASSERTION_MSG(circulator_size(fitr->facet_begin())==3, "(ERROR) Non-triangle facet!");
    if(refine_border) {
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, p), p)) {
        if(is_border(opposite(hd, p),p)) {
          face_halfedge_border[face_id] = hd;
          Point vpt;
          mask.border_node(hd, cpt[i+1], cpt[i], vpt);
          new_positions.push_back(std::make_pair(target(hd,p), vpt));
          i += 2;

          // the border subdivision is only performed every second subdivision
          // step and there can thus only be one border edge per face
          break;
        }
      }

      if(face_halfedge_border[face_id] == boost::graph_traits<Poly>::null_halfedge())
        mask.face_node(fd, cpt[i++]);
    } else {
      mask.face_node(fd, cpt[i++]);
    }

    ++face_id;
  }

  // smooth the position of existing vertices
  for(vertex_descriptor vd : p_vertices){
    Point pt;
    if(!is_border(vd, p)) {
      mask.vertex_node(vd, pt);
      new_positions.push_back(std::make_pair(vd, pt));
    }
  }

  // insert the new subdividing points
  typename std::vector<face_descriptor>::iterator fit=p_faces.begin();
  for(std::size_t i=0, cpt_id=0; i < num_f; ++i, ++fit){
    face_descriptor fd = *fit;
    halfedge_descriptor hd = face_halfedge_border[i];
    if(refine_border && hd != boost::graph_traits<Poly>::null_halfedge()) {
      halfedge_descriptor hd_next = next(hd, p);
      halfedge_descriptor new_e1 = Euler::split_edge(hd, p);
      halfedge_descriptor new_e2 = Euler::split_edge(hd, p);
      put(vpm, target(new_e1, p), cpt[cpt_id++]);
      put(vpm, target(new_e2, p), cpt[cpt_id++]);
      Euler::split_face(new_e1, hd_next, p);
      Euler::split_face(new_e2, hd_next, p);
    } else {
      halfedge_descriptor center = Euler::add_center_vertex(halfedge(fd,p),p);
      put(vpm, target(center,p), cpt[cpt_id++]);
    }
  }

  // actually inserts the new positions in the vertex point property map
  typename std::list<std::pair<vertex_descriptor, Point> >::iterator it = new_positions.begin(),
      end = new_positions.end();
  for(; it!=end; ++it) {
    put(vpm, it->first, it->second);
  }

  // flip the old edges (except the border edges)
  typename std::vector<edge_descriptor>::iterator eitr = p_edges.begin();
  for (std::size_t i = 0; i < num_e; ++i) {
    halfedge_descriptor e = halfedge(*eitr,p);
    ++eitr; // move to next edge before flip since flip destroys current edge
    if (! is_border(edge(e,p),p)) {
      halfedge_descriptor h = Euler::join_face(e,p);
      Euler::split_face(prev(h,p), next(h,p),p);
    }
  }

  CGAL_postcondition(CGAL::is_valid_polygon_mesh(p));
  delete []cpt;
}

} // namespace internal
} // namespace Subdivision_method_3
} // namespace CGAL

#endif //CGAL_SUBDIVISION_HOSTS_IMPL_3_H
