// ======================================================================
//
// Copyright (c) 2005-2011 GeometryFactory (France).  All Rights Reserved.
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

#ifndef CGAL_POLYHEDRON_SUBDIVISION_IMPL_H_02102006
#define CGAL_POLYHEDRON_SUBDIVISION_IMPL_H_02102006

#include <CGAL/basic.h>

#include <vector>

#include <CGAL/circulator.h>
#include <CGAL/Polyhedron_decorator_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>

namespace CGAL {

// ======================================================================
namespace Subdivision_method_3 {
  
  namespace Private {
  // ======================================================================
  template <class Poly, template <typename> class Mask>
  void PQQ_1step(Poly& p, Mask<Poly> mask) {
    typedef Polyhedron_decorator_3<Poly>           PD;

    typedef typename boost::graph_traits<Poly>::vertex_descriptor           vertex_descriptor;
    typedef typename boost::graph_traits<Poly>::halfedge_descriptor         halfedge_descriptor;
    typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;

    typedef typename boost::graph_traits<Poly>::vertex_iterator         vertex_iterator;
    typedef typename boost::graph_traits<Poly>::edge_iterator           edge_iterator;
    typedef typename boost::graph_traits<Poly>::face_iterator          face_iterator;

    typedef Halfedge_around_face_circulator<Poly>  Halfedge_around_facet_circulator;

    typedef typename boost::property_map<Poly, vertex_point_t>::type Vertex_pmap;
    typedef typename boost::property_traits<Vertex_pmap>::value_type Point;
        
    Vertex_pmap vpm = get(CGAL::vertex_point, p);

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

    // If Polyhedron is using vector, we need to reserve the memory to prevent 
    // the CGAL_assertion.
    // This function for polyhedron using list is VOID.
    p.reserve(num_vertex+num_edge+num_facet, 4*2*num_edge, 4*num_edge/2);

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
      mask.facet_node(*fitr, face_point_buffer[i]);


    {
      std::size_t i = 0;
      BOOST_FOREACH(edge_descriptor ed, edges(p)){
        if(is_border(ed,p)){
          int v = v_index[target(ed,p)];
          v_onborder[v] = true;
          mask.border_node(halfedge(ed,p), edge_point_buffer[i], vertex_point_buffer[v]);
          
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
      vertex_descriptor vh = PD::insert_vertex(p, halfedge(*eitr,p));
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
      halfedge_descriptor newe = PD::insert_edge(p, e1, e2);

      // Step 3.
      halfedge_descriptor newv = PD::insert_vertex_return_edge(p, newe);
      newv = prev(opposite(newv,p),p); // change newv to the larger face and 
      // still points to the newly inserted 
      // vertex
      // Update the geometry data of the newly inserted face-vertices
      put(vpm, target(newv,p), face_point_buffer[i]);

      // Step 4.
      while (hcir != hcir_begin) {
        e1 = * ++hcir;
        ++hcir; // Must move the cir before inserts the new edge !!
        PD::insert_edge(p, e1, newv); 
      }
    }

    // Update the geometry data of the newly inserted vertices by the 
    // vertices buffer
    vitr = vertices(p).first;
    for (size_t i = 0; i < num_vertex; i++, ++vitr) 
      put(vpm, *vitr, vertex_point_buffer[i]);

    delete []vertex_point_buffer;
  }

  // ======================================================================
  template <class Poly, template <typename> class Mask>
  void PTQ_1step(Poly& p, Mask<Poly> mask) {

    typedef Polyhedron_decorator_3<Poly>           PD;

    typedef typename boost::graph_traits<Poly>::vertex_descriptor           vertex_descriptor;
    typedef typename boost::graph_traits<Poly>::halfedge_descriptor         halfedge_descriptor;
    typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;

    typedef typename boost::graph_traits<Poly>::vertex_iterator         vertex_iterator;
    typedef typename boost::graph_traits<Poly>::edge_iterator           edge_iterator;
    typedef typename boost::graph_traits<Poly>::face_iterator          face_iterator;

    typedef Halfedge_around_face_circulator<Poly>  Halfedge_around_face_circulator;

    typedef typename boost::property_map<Poly, vertex_point_t>::type Vertex_pmap;
    typedef typename boost::property_traits<Vertex_pmap>::value_type Point;

    Vertex_pmap vpm = get(CGAL::vertex_point, p);

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

    // If Polyhedron is using vector, we need to reserve the memory to prevent 
    // the CGAL_assertion.
    // This function for polyhedron using list is VOID.
    p.reserve(num_vertex+num_edge, 2*2*num_edge, 4*num_edge/2);

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
          int v = v_index[target(ed,p)];
          v_onborder[v] = true;
          mask.border_node(halfedge(ed,p), edge_point_buffer[i], vertex_point_buffer[v]);
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
      vertex_descriptor vh = PD::insert_vertex(p, halfedge(*eitr,p));
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
      e2 = PD::insert_edge(p, e1, e2);
      e3 = PD::insert_edge(p, e2, e3);
      PD::insert_edge(p, e3, e1);
    }

    // Update the geometry data of the newly inserted vertices by the 
    // vertices buffer
    vitr = vertices(p).first;
    for (size_t i = 0; i < num_vertex; i++, ++vitr)
      put(vpm, *vitr, vertex_point_buffer[i]);

    delete []vertex_point_buffer;
    }


  // ======================================================================
//#define CGAL_EULER_DQQ_SPLITTING
//#define CGAL_EULER_DQQ_TILTING   // Tilting is faster
  template <class Poly, template <typename> class Mask>
  void DQQ_1step(Poly& p, Mask<Poly> mask) {

    typedef Polyhedron_decorator_3<Poly>           PD;

    typedef typename boost::graph_traits<Poly>::vertex_descriptor           vertex_descriptor;
    typedef typename boost::graph_traits<Poly>::halfedge_descriptor         halfedge_descriptor;
    typedef typename boost::graph_traits<Poly>::edge_descriptor         edge_descriptor;

    typedef typename boost::graph_traits<Poly>::vertex_iterator         vertex_iterator;
    typedef typename boost::graph_traits<Poly>::edge_iterator           edge_iterator;
    typedef typename boost::graph_traits<Poly>::face_iterator          face_iterator;

    typedef Halfedge_around_face_circulator<Poly>  Halfedge_around_face_circulator;

    typedef typename boost::property_map<Poly, vertex_point_t>::type Vertex_pmap;
    typedef typename boost::property_traits<Vertex_pmap>::value_type Point;

    Vertex_pmap vpm = get(CGAL::vertex_point, p);

    typename boost::graph_traits<Poly>::vertices_size_type num_v = num_vertices(p);
    typename boost::graph_traits<Poly>::halfedges_size_type num_e = num_halfedges(p)/2;
    typename boost::graph_traits<Poly>::faces_size_type num_f = num_faces(p);

    size_t num_be ;// AF= p.size_of_border_edges();

    Point* point_buffer = new Point[num_e*2];

    //
#ifdef CGAL_EULER_DQQ_SPLITTING
    //
    // Splitting

    //! Splitting is not implemented to support border

    // build the point_buffer
    Facet_iterator fitr, fitr_end = p.facets_end();
    int pi = 0;
    for (fitr = p.facets_begin(); fitr != fitr_end; ++fitr) {
      Halfedge_around_face_circulator cir = fitr->facet_begin();
      do {
        mask.corner_node(cir, point_buffer[pi++]);
      } while (--cir != fitr->facet_begin());
    }

    // If Polyhedron is using vector, we need to reserve the memory to prevent 
    // the CGAL_assertion. This function for polyhedron using list is VOID.
    p.reserve(num_v+num_e+num_f, 2*num_e, (2+4+2)*num_e);

    // Build the connectivity using insert_vertex() and insert_edge()
    // 1. create barycentric centers of each facet
    fitr = p.facets_begin();
    pi = 0;
    for (size_t i = 0; i < num_f; i++) {
      Facet_handle fh = fitr;
      ++fitr;
      Vertex_handle vh = (p.create_center_vertex(fh->facet_begin()))->vertex();

      // 1.1 add vertex on each new edges
      Halfedge_around_vertex_circulator vcir = vh->vertex_begin();  
      int vn = circulator_size(vcir);
      for (int j = 0; j < vn; ++j) {
        Halfedge_handle e = vcir;
        ++vcir;
        Vertex_handle v = PD::insert_vertex(p, e);
        v->point() = point_buffer[pi++];
      }
      // 1.2 connect new vertices surround each barycentric center
      for (int j = 0; j < vn; ++j) {
        Halfedge_handle e1 = vcir->prev();
        ++vcir;
        Halfedge_handle e2 = vcir->opposite();
        PD::insert_edge(p, e1, e2);
      }
      // 1.3 remove the barycentric centers
      p.erase_center_vertex(vcir);
    }

    // 2. remove old edges
    Edge_iterator eitr = p.edges_begin();
    for (size_t i = 0; i < num_e; ++i) {
      Halfedge_handle eh = eitr;
      ++eitr;
      p.join_facet(eh);
    }

    // 3. connect new vertices surround old vertices and then remove 
    //    old vertices.
    vertex_iterator vitr = p.vertices_begin();
    for (size_t i = 0; i < num_v; ++i) {
      Halfedge_around_vertex_circulator vcir = vitr->vertex_begin();  
      int vn = circulator_size(vcir);
      for (int j = 0; j < vn; ++j) {
        Halfedge_handle e1 = vcir->prev();
        ++vcir;
        Halfedge_handle e2 = vcir->opposite();
        PD::insert_edge(p, e1, e2);
      }
      ++vitr;
      p.erase_center_vertex(vcir);    
    }

    //
#else
    //
    // Tilting

    // build the point_buffer
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

    // If Polyhedron is using vector, we need to reserve the memory to prevent 
    // the CGAL_assertion. This function for polyhedron using list is VOID.
    p.reserve(num_v+num_e+num_f, 2*num_e, (2+4+2)*num_e);

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
          vertex_descriptor v = PD::insert_vertex(p, e);
          put(vpm, v, point_buffer[pi++]);
        }
      }

      vcir = Halfedge_around_target_circulator<Poly>(vh,p);
      for (size_t j = 0; j < vn; ++j) {
        if (! is_border(*vcir,p)) {
          halfedge_descriptor e1 = * CGAL::cpp11::prev(vcir);
          ++vcir;
          if (! is_border(*vcir,p)) {
            halfedge_descriptor e2 = opposite(*vcir,p);
            PD::insert_edge(p, e1, e2);
          }
        } else ++vcir;
      }
      //p.erase_center_vertex(vh->vertex_begin());
    }

    edge_iterator eitr = edges(p).first;
    for (size_t i = 0; i < num_e; ++i) {
      halfedge_descriptor eh = halfedge(*eitr,p);
      ++eitr;
      if (! is_border(edge(eh,p),p)) {
        PD::insert_edge(p, prev(prev(eh,p),p), eh);
        eh = opposite(eh,p);
        PD::insert_edge(p, prev(prev(eh,p),p), eh);
        Euler::join_face(eh,p);
      } else {
        if (is_border(eh,p)) {
          eh = opposite(eh,p);
          PD::insert_edge(p, eh, prev(prev(eh,p),p));
        } else 
          PD::insert_edge(p, prev(prev(eh,p),p), eh);
      }
    }

    // after this point, the original border edges are in front!
    eitr = edges(p).first;
    for (size_t i = 0; i < num_be; ++i) {
      halfedge_descriptor eh = halfedge(*eitr,p);
      ++eitr;

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

#endif //CGAL_EULER_DQQ_SPLITTING

    delete []point_buffer;
  }

  // ======================================================================
  template <class Poly, template <typename> class Mask>
  void Sqrt3_1step(Poly& p, Mask<Poly> mask) {

    typedef typename boost::graph_traits<Poly>::vertex_descriptor         vertex_descriptor;
    typedef typename boost::graph_traits<Poly>::halfedge_descriptor         halfedge_descriptor;
    typedef typename boost::graph_traits<Poly>::face_descriptor         face_descriptor;

    typedef typename boost::graph_traits<Poly>::vertex_iterator         vertex_iterator;
    typedef typename boost::graph_traits<Poly>::edge_iterator           edge_iterator;
    typedef typename boost::graph_traits<Poly>::face_iterator          face_iterator;

    typedef typename boost::property_map<Poly, vertex_point_t>::type Vertex_pmap;
    typedef typename boost::property_traits<Vertex_pmap>::value_type Point;

    Vertex_pmap vpm = get(CGAL::vertex_point, p);


    typename boost::graph_traits<Poly>::vertices_size_type num_v = num_vertices(p);
    typename boost::graph_traits<Poly>::halfedges_size_type num_e = num_halfedges(p)/2;
    typename boost::graph_traits<Poly>::faces_size_type num_f = num_faces(p);

    p.reserve(num_v+num_f, (num_e+3*num_f)*2, 3*num_f);

    // prepare the smoothed center points
    Point* cpt = new Point[num_f]; 
    
    std::size_t i = 0;
    BOOST_FOREACH (face_descriptor fd, faces(p)) {
      //ASSERTION_MSG(circulator_size(fitr->facet_begin())==3, "(ERROR) Non-triangle facet!");
      mask.facet_node(fd, cpt[i++]);
    }

    // smooth the vertex points
    BOOST_FOREACH(vertex_descriptor vd, vertices(p)){
      Point p;
      mask.vertex_node(vd,p);
      put(vpm,vd,p);
    }


    // insert the facet points
    face_iterator b,e;
    boost::tie(b,e) = faces(p);
    for(std::size_t i=0 ; i < num_f; ++i, ++b){
      face_descriptor fd = *b;
      halfedge_descriptor center = Euler::add_center_vertex(halfedge(fd,p),p);
      put(vpm, target(center,p), cpt[i]);
    }

    delete []cpt;

    // flip the old edges except the border edges
    edge_iterator eitr = edges(p).first;
    for (size_t i = 0; i < num_e; ++i) {
      halfedge_descriptor e = halfedge(*eitr,p);
      ++eitr; // move to next edge before flip since flip destroys current edge
      if (! is_border(edge(e,p),p)) {
        halfedge_descriptor h = Euler::join_face(e,p);
        Euler::split_face(prev(h,p), next(h,p),p);
      }
    }

    // TODO: border ...

    CGAL_postcondition(p.is_valid());
  }
  }
}

} //namespace CGAL

#endif //CGAL_POLYHEDRON_SUBDIVISION_H_01292002
