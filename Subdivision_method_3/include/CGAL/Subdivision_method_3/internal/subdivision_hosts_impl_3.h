// ======================================================================
//
// Copyright (c) 2005-2017 GeometryFactory (France).  All Rights Reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//
// ======================================================================

#ifndef CGAL_SUBDIVISION_HOSTS_IMPL_3_H
#define CGAL_SUBDIVISION_HOSTS_IMPL_3_H

#include <CGAL/Subdivision_method_3/internal/Euler_extensions.h>

#include <CGAL/basic.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/circulator.h>
#include <CGAL/tags.h>

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/if.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <iterator>
#include <list>
#include <vector>

namespace CGAL {

namespace Subdivision_method_3 {

namespace internal {

template <class Poly, class VertexPointMap, class Mask>
void PQQ_1step(Poly& p, VertexPointMap vpm, Mask mask) {
  typedef typename boost::graph_traits<Poly>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;

  typedef typename boost::graph_traits<Poly>::vertex_iterator         vertex_iterator;
  typedef typename boost::graph_traits<Poly>::edge_iterator           edge_iterator;
  typedef typename boost::graph_traits<Poly>::face_iterator           face_iterator;

  typedef Halfedge_around_face_circulator<Poly>  Halfedge_around_facet_circulator;

  // Build a new vertices buffer has the following structure
  //
  // 0 1 ... e_begin ... f_begin ... (end_of_buffer)
  // 0 ... e_begin-1       : store the positions of the vertex-vertices
  // e_begin ... f_begin-1 : store the positions of the edge-vertices
  // f_begin ... (end)     : store the positions of the face-vertices
  // The index of the vertices buffer should 1-1 map to the distance
  // of the corresponding iterator to the begin of the iterator.
  typename boost::graph_traits<Poly>::vertices_size_type num_vertex = num_vertices(p);
  typename boost::graph_traits<Poly>::halfedges_size_type num_edge = num_halfedges(p)/2;
  typename boost::graph_traits<Poly>::faces_size_type num_facet = num_faces(p);

  // We need to reserve the memory to prevent reallocation.
  reserve(p,num_vertex+num_edge+num_facet, 4*2*num_edge, 4*num_edge/2);

  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  Point* vertex_point_buffer = new Point[num_vertex + num_edge + num_facet];
  Point* edge_point_buffer = vertex_point_buffer + num_vertex;
  Point* face_point_buffer = edge_point_buffer + num_edge;

  int i=0;
  boost::unordered_map<vertex_descriptor,int> v_index;
  BOOST_FOREACH(vertex_descriptor vh, vertices(p)){
    v_index[vh]= i++;
  }

  std::vector<bool> v_onborder(num_vertex);
  face_iterator fitr = faces(p).first;
  for (size_t i = 0; i < num_facet; i++, ++fitr)
    mask.face_node(*fitr, face_point_buffer[i]);


  {
    std::size_t i = 0;
    BOOST_FOREACH(edge_descriptor ed, edges(p)){
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

  vertex_iterator vitr = vertices(p).first;
  for (size_t i = 0; i < num_vertex; i++, ++vitr)
    if (!v_onborder[v_index[*vitr]]) mask.vertex_node(*vitr, vertex_point_buffer[i]);

  // Build the connectivity using insert_vertex() and insert_edge()
  // 1. insert_vertex() to all edges and set them to new positions
  // 2. insert_edge() between 2 randomly selected neighboring new inserted
  //    vertices
  // 3. insert_vertex() to the new inserted edge and set them to new positions
  // 4. insert_edge() between all other new inserted vertices of step 1 and
  //    the new inserted vertex of step 3
  // Step 1.
  edge_iterator eitr = edges(p).first;
  for (size_t i = 0; i < num_edge; i++, ++eitr) {
    vertex_descriptor vh = insert_vertex(p, halfedge(*eitr,p));
    put(vpm, vh, edge_point_buffer[i]);
  }
  fitr = faces(p).first;

  // TODO: the topoloy modification can be done by a template function
  //       and that gives the user a chance to create new topological masks.
  for (size_t i = 0; i < num_facet; i++, ++fitr) {
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
  vitr = vertices(p).first;
  for (size_t i = 0; i < num_vertex; i++, ++vitr)
    put(vpm, *vitr, vertex_point_buffer[i]);

//  CGAL_postcondition(p.is_valid());
  delete []vertex_point_buffer;
}

// ======================================================================
template <class Poly,class VertexPointMap, class Mask>
void PTQ_1step(Poly& p, VertexPointMap vpm, Mask mask) {
  typedef typename boost::graph_traits<Poly>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;

  typedef typename boost::graph_traits<Poly>::vertex_iterator         vertex_iterator;
  typedef typename boost::graph_traits<Poly>::edge_iterator           edge_iterator;
  typedef typename boost::graph_traits<Poly>::face_iterator           face_iterator;

  typedef Halfedge_around_face_circulator<Poly>  Halfedge_around_face_circulator;

  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  // Build a new vertices buffer has the following structure
  //
  // 0 1 ... e_begin ... f_begin ... (end_of_buffer)
  // 0 ... e_begin-1       : store the positions of the vertex-vertices
  // e_begin ... (end)     : store the positions of the edge-vertices
  // The index of the vertices buffer should 1-1 map to the distance
  // of the corresponding iterator to the begin of the iterator.
  typename boost::graph_traits<Poly>::vertices_size_type num_vertex = num_vertices(p);
  typename boost::graph_traits<Poly>::halfedges_size_type num_edge = num_halfedges(p)/2;
  typename boost::graph_traits<Poly>::faces_size_type num_facet = num_faces(p);

  // We need to reserve the memory to prevent reallocation.
  reserve(p,num_vertex+num_edge, 2*2*num_edge, 4*num_edge/2);

  Point* vertex_point_buffer = new Point[num_vertex + num_edge];
  Point* edge_point_buffer = vertex_point_buffer + num_vertex;

  int i=0;
  boost::unordered_map<vertex_descriptor,int> v_index;
  BOOST_FOREACH(vertex_descriptor vh, vertices(p)){
    v_index[vh]= i++;
  }
  std::vector<bool> v_onborder(num_vertex);

  {
    std::size_t i = 0;
    BOOST_FOREACH(edge_descriptor ed, edges(p)){
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
  vertex_iterator vitr = vertices(p).first;
  for (size_t i = 0; i < num_vertex; i++, ++vitr)
    if (!v_onborder[i]) mask.vertex_node(*vitr, vertex_point_buffer[i]);

  // Build the connectivity using insert_vertex() and insert_edge()
  // 1. insert_vertex() to all edges and set them to new positions
  // 2. insert_edge() between 2 randomly selected neighboring new inserted
  //    vertices
  // 3. insert_vertex() to the new inserted edge and set them to new positions
  // 4. insert_edge() between all other new inserted vertices of step 1 and
  //    the new inserted vertex of step 3
  // Step 1.
  edge_iterator eitr = edges(p).first;
  for (size_t i = 0; i < num_edge; i++, ++eitr) {
    vertex_descriptor vh = insert_vertex(p, halfedge(*eitr,p));
    put(vpm,vh, edge_point_buffer[i]);
  }
  face_iterator fitr = faces(p).first;
  for (size_t i = 0; i < num_facet; i++, ++fitr) {
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
  vitr = vertices(p).first;
  for (size_t i = 0; i < num_vertex; i++, ++vitr)
    put(vpm, *vitr, vertex_point_buffer[i]);

//  CGAL_postcondition(p.is_valid());
  delete []vertex_point_buffer;
}


// ======================================================================
template <class Poly, class VertexPointMap, class Mask>
void DQQ_1step_impl(Poly& p, VertexPointMap vpm, Mask mask, CGAL::Tag_false) {
  typedef typename boost::graph_traits<Poly>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;

  typedef typename boost::graph_traits<Poly>::vertex_iterator         vertex_iterator;
  typedef typename boost::graph_traits<Poly>::edge_iterator           edge_iterator;

  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  typename boost::graph_traits<Poly>::vertices_size_type num_v = num_vertices(p);
  typename boost::graph_traits<Poly>::halfedges_size_type num_e = num_halfedges(p)/2;
  typename boost::graph_traits<Poly>::faces_size_type num_f = num_faces(p);

  std::vector<halfedge_descriptor> border_halfedges;
  size_t num_be = 0 ;
  BOOST_FOREACH(edge_descriptor ed, edges(p)){
    if(is_border(ed,p)){
      ++num_be;
      border_halfedges.push_back(halfedge(ed,p));
    }
  }
  Point* point_buffer = new Point[num_e*2];

  // Build the point_buffer
  vertex_iterator vitr, vitr_end;
  boost::tie(vitr,vitr_end) = vertices(p);
  int pi = 0;
  BOOST_FOREACH(vertex_descriptor vd, vertices(p)){
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(vd,p)){
      if (! is_border(hd,p)){
        mask.corner_node(hd, point_buffer[pi++]);
      }
    }
  }

  // Reserve to avoid rellocations during insertions
  reserve(p,num_v+num_e+num_f, 2*num_e, (2+4+2)*num_e);

  // Build the connectivity using insert_vertex() and insert_edge()
  pi = 0;
  for (size_t i = 0; i < num_v; ++i) {
    vertex_descriptor vh = *vitr;
    ++vitr;

    Halfedge_around_target_circulator<Poly> vcir(vh,p);
    size_t vn = degree(vh,p);
    for (size_t j = 0; j < vn; ++j) {
      halfedge_descriptor e = *vcir;
      ++vcir;
      if (! is_border(e,p)) {
        vertex_descriptor v = insert_vertex(p, e);
        put(vpm, v, point_buffer[pi++]);
      }
    }

    vcir = Halfedge_around_target_circulator<Poly>(vh,p);
    for (size_t j = 0; j < vn; ++j) {
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

  edge_iterator eitr = edges(p).first;
  for (size_t i = 0; i < num_e; ++i) {
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

  // After this point, the original border edges are in front!
  BOOST_FOREACH(halfedge_descriptor eeh, border_halfedges){
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

  vitr = vertices(p).first;
  for (size_t i = 0; i < num_v-num_be; ++i) {
    vertex_descriptor vh = *vitr;
    ++vitr;
    Euler::remove_center_vertex(halfedge(vh,p),p);
  }

  delete []point_buffer;
}

template <class Poly, class VertexPointMap, class Mask>
void DQQ_1step_impl(Poly& p, VertexPointMap vpm, Mask mask, CGAL::Tag_true) {
  typedef typename boost::graph_traits<Poly>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<Poly>::face_descriptor         face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  // Note that for types like 'Surface_mesh', num_vertices() and other similar functions
  // return the TOTAL number of elements, which may include removed vertices.
  typename boost::graph_traits<Poly>::vertices_size_type num_v = num_vertices(p);
  typename boost::graph_traits<Poly>::edges_size_type num_e = num_edges(p);
  typename boost::graph_traits<Poly>::faces_size_type num_f = num_faces(p);

  // Move `p` into `moved_p`, and build the subdivided mesh from scratch in `p`.
  // This is done to make the algorithm work with CGAL::Surface_mesh,
  // even though CGAL::Surface_mesh does not insert elements at the end (due to removed elements).
  // The DooSabin subdivision of a mesh is a completely different mesh so there
  // is no additional cost to rebuild from scratch (but there is a bit from
  // using `copy_face_graph`).
  Poly moved_p;
  reserve(moved_p,num_v, num_e, num_f);

  // We must use copy_face_graph rather than an assignement operator because
  // we need the correspondence between vertex_descriptors
  boost::unordered_map<vertex_descriptor, vertex_descriptor> v2v(num_v);
  CGAL::copy_face_graph(p, moved_p, std::inserter(v2v, v2v.end()));

  VertexPointMap moved_vpm = get(vertex_point, moved_p);

  // Move the position information to the internal property map of moved_p
  typename boost::unordered_map<vertex_descriptor, vertex_descriptor>::iterator it = v2v.begin(),
      end = v2v.end();
  for(; it!=end; ++it) {
    put(moved_vpm, it->second, get(vpm, it->first));
  }

  // Temporarily change the members of the mask to `moved_p`
  mask.pmesh = &moved_p;
  mask.vpm = moved_vpm;

  clear(p);
  reserve(p,num_v+num_e+num_f, 2*num_e, (2+4+2)*num_e);

  // Correspondence between halfedges of the original mesh and some of the
  // halfedges in the subdivided mesh. Since we have halfedge_index_t,
  // we can simply use a vector!
  // Note: need to make sure the halfedge index map is initialized
  // @fixme this overwrites previous initilizations...
  helpers::init_halfedge_indices(const_cast<Poly&>(moved_p),
                                  get(boost::halfedge_index, moved_p));
  std::vector<halfedge_descriptor> old_to_new(2 * num_e);
  Property_map_binder<typename boost::property_map<Poly, boost::halfedge_index_t>::type,
                      typename Pointer_property_map<halfedge_descriptor>::type>
  hmap = bind_property_maps(get(boost::halfedge_index, moved_p),
                            make_property_map(old_to_new));

  // Build new n-faces
  BOOST_FOREACH(face_descriptor fd, faces(moved_p)) {
    halfedge_descriptor hd = halfedge(fd, moved_p);
    std::list<vertex_descriptor> vertices_of_new_face;

    // Keep the first outside; it will be used to build the correspondence
    // between old and new halfedges
    Point first_pt;
    mask.corner_node(hd, first_pt);
    vertex_descriptor first_v = add_vertex(p);
    put(vpm, first_v, first_pt);
    vertices_of_new_face.push_back(first_v);

    // Loop normally and add the rest of the vertices
    halfedge_descriptor done = hd;
    hd = next(hd, moved_p);
    while(hd != done) {
      Point pt;
      mask.corner_node(hd, pt);
      vertex_descriptor v = add_vertex(p);
      put(vpm, v, pt);
      vertices_of_new_face.push_back(v);
      hd = next(hd, moved_p);
    }

    face_descriptor new_face = Euler::add_face(vertices_of_new_face, p);

    // Find the starting halfedge in the new face that corresponds to halfedge(fd, p)
    halfedge_descriptor nf_hd = halfedge(new_face, p);
    while(target(nf_hd, p) != first_v) {
      nf_hd = next(nf_hd, p);
    }

    // Build the correspondence between old and new halfedges
    hd = halfedge(fd, moved_p);
    done = nf_hd;
    do {
      put(hmap, hd, nf_hd);
      hd = next(hd, moved_p);
      nf_hd = next(nf_hd, p);
    } while (nf_hd != done);
  }

  // Build new edge-faces
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(moved_p)) {
    if(is_border(hd, moved_p))
      continue;

    halfedge_descriptor hd_opp = opposite(hd, moved_p);
    if(is_border(hd_opp, moved_p))
      continue;

    if(hd > hd_opp)
      continue;

    halfedge_descriptor new_hd = opposite(get(hmap, hd), p);
    halfedge_descriptor new_hd_opp = opposite(get(hmap, hd_opp), p);

    boost::array<vertex_descriptor, 4> v = {{source(new_hd, p),
                                             target(new_hd, p),
                                             source(new_hd_opp, p),
                                             target(new_hd_opp, p)}};

    Euler::add_face(v, p);
  }

  // Build new vertex-faces
  BOOST_FOREACH(vertex_descriptor vd, vertices(moved_p)) {
    if(is_border(vd, moved_p))
      continue;

    halfedge_descriptor hd = halfedge(vd, moved_p);
    halfedge_descriptor new_hd = opposite(get(hmap, hd), p);
    halfedge_descriptor new_face_hd = opposite(prev(new_hd, p), p), done = new_face_hd;
    std::list<vertex_descriptor> vertices_of_new_faces;
    do {
      vertices_of_new_faces.push_back(source(new_face_hd, p));
      new_face_hd = next(new_face_hd, p);
    } while(new_face_hd != done);

    Euler::add_face(vertices_of_new_faces, p);
  }

  // Reset the members of the mask
  mask.pmesh = &p;
  mask.vpm = vpm;
}

template <class Poly, class VertexPointMap, class Mask>
void DQQ_1step(Poly& p, VertexPointMap vpm, Mask mask) {
  // Check if halfedges are index-based, which allows to use vectors instead of maps
  DQQ_1step_impl(p, vpm, mask,
                 boost::graph_has_property<Poly, boost::halfedge_index_t>());
//  CGAL_postcondition(p.is_valid());
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

  typedef typename boost::graph_traits<Poly>::edge_iterator           edge_iterator;
  typedef typename boost::graph_traits<Poly>::face_iterator           face_iterator;

  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  typename boost::graph_traits<Poly>::vertices_size_type num_v = num_vertices(p);
  typename boost::graph_traits<Poly>::halfedges_size_type num_e = num_halfedges(p)/2;
  typename boost::graph_traits<Poly>::faces_size_type num_f = num_faces(p);

  // reserve enough size for the new points
  typename boost::graph_traits<Poly>::faces_size_type new_pts_size = num_f;
  if(refine_border) {
    BOOST_FOREACH(edge_descriptor ed, edges(p)){
      if(is_border(ed, p))
        ++new_pts_size;
    }
  }
  Point* cpt = new Point[new_pts_size];

  // size of the subdivided mesh
  reserve(p,num_v + new_pts_size, (num_e + 2*num_f + new_pts_size)*2, 3*num_f);

  // keep in memory whether a face is incident to the border and, if so, which
  // halfedge corresponds to THE (there can only be one) border edge.
  std::vector<halfedge_descriptor> face_halfedge_border(num_f,
                                                        boost::graph_traits<Poly>::null_halfedge());
  std::list<std::pair<vertex_descriptor, Point> > new_positions;

  // compute the positions of new points
  std::size_t i = 0;
  std::size_t face_id = 0;
  BOOST_FOREACH (face_descriptor fd, faces(p)) {
    //ASSERTION_MSG(circulator_size(fitr->facet_begin())==3, "(ERROR) Non-triangle facet!");
    if(refine_border) {
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, p), p)) {
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
  BOOST_FOREACH(vertex_descriptor vd, vertices(p)){
    Point pt;
    if(!is_border(vd, p)) {
      mask.vertex_node(vd, pt);
      new_positions.push_back(std::make_pair(vd, pt));
    }
  }

  // insert the new subdividing points
  face_iterator b,e;
  boost::tie(b,e) = faces(p);
  for(std::size_t i=0, cpt_id=0; i < num_f; ++i, ++b){
    face_descriptor fd = *b;
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
  edge_iterator eitr = edges(p).first;
  for (size_t i = 0; i < num_e; ++i) {
    halfedge_descriptor e = halfedge(*eitr,p);
    ++eitr; // move to next edge before flip since flip destroys current edge
    if (! is_border(edge(e,p),p)) {
      halfedge_descriptor h = Euler::join_face(e,p);
      Euler::split_face(prev(h,p), next(h,p),p);
    }
  }

//  CGAL_postcondition(p.is_valid());
  delete []cpt;
}

} // namespace internal
} // namespace Subdivision_method_3
} // namespace CGAL

#endif //CGAL_SUBDIVISION_HOSTS_IMPL_3_H
