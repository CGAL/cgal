// Copyright (c) 2010, 2012 GeometryFactory Sarl (France).
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
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_MESH_3_DETECT_POLYLINES_IN_POLYHEDRA_H
#define CGAL_MESH_3_DETECT_POLYLINES_IN_POLYHEDRA_H

#include <CGAL/Mesh_3/Detect_polylines_in_polyhedra_fwd.h>
#include <CGAL/Has_timestamp.h>
#include <CGAL/Default.h>

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/mpl/if.hpp>

namespace CGAL { namespace Mesh_3 {

template <typename Handle>
struct CGAL_with_time_stamp
{
public:
  static bool less(Handle h1, Handle h2)
  {
    return h1->time_stamp() < h2->time_stamp();
  }
};

template <typename Handle>
struct CGAL_no_time_stamp
{
public:
  static bool less(Handle h1, Handle h2)
  {
    return &*h1 < &*h2;
  }
};

struct Detect_polyline_less
{
  template<typename Handle>
  bool operator()(const Handle& h1, const Handle& h2) const
  {
    typedef typename std::iterator_traits<Handle>::value_type Type;
    typedef typename boost::mpl::if_c<
        CGAL::internal::Has_timestamp<Type>::value,
        CGAL_with_time_stamp<Handle>,
        CGAL_no_time_stamp<Handle> >::type    Comparator;
    return Comparator::less(h1, h2);
  }
};

template <typename Polyhedron>
struct Detect_polylines {
  typedef typename Polyhedron::Traits Geom_traits;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::size_type size_type;

  typedef std::set<Vertex_handle, Detect_polyline_less> Vertices_set;
  typedef std::map<Vertex_handle, 
                   size_type,                    
                   Detect_polyline_less> Vertices_counter;

  typedef std::set<Halfedge_handle, Detect_polyline_less> Feature_edges_set;

  Feature_edges_set edges_to_consider;
  Vertices_set corner_vertices;

  // typedef std::vector<Point_3> Polyline_and_context;

  typedef typename Polyhedron::Vertex Polyhedron_vertex;
  typedef typename Polyhedron_vertex::Set_of_indices Set_of_indices;

  template <typename T>
  static 
  void display_index(std::ostream& stream, const T& x)
  {
    stream << x;
  }

  template <typename T, typename U>
  static 
  void display_index(std::ostream& stream, const std::pair<T,U>& p)
  {
    stream << p.first << "+" << p.second;
  }

  static 
  void display_set(std::ostream& stream, Set_of_indices set) {
    stream << "( ";
    BOOST_FOREACH(typename Set_of_indices::value_type i, set) {
      display_index(stream, i);
      stream << " ";
    }
    stream << ")";
  }

  static Set_of_indices
  edge_indices(const Halfedge_handle he) {
    Set_of_indices set_of_indices;
    const Set_of_indices& source_set = 
      he->opposite()->vertex()->incident_patches_ids_set();
    const Set_of_indices& target_set = 
      he->vertex()->incident_patches_ids_set();
    std::set_intersection(source_set.begin(), source_set.end(),
                          target_set.begin(), target_set.end(),
                          std::inserter(set_of_indices,
                                        set_of_indices.begin()));
    if(set_of_indices.empty()) {
      std::cerr << "Indices set of following edge is empty:\n";
      std::cerr << "  " << he->opposite()->vertex()->point()
                << " ";
      display_set(std::cerr, source_set);
      std::cerr << "\n";
      std::cerr << "  " << he->vertex()->point()
                << " ";
      display_set(std::cerr, target_set);
      std::cerr << "\n";
    }
    return set_of_indices;
  }

  static Halfedge_handle canonical(Halfedge_handle he)
  {
    const Halfedge_handle& op = he->opposite();
    if(Detect_polyline_less()(he, op))
      return he;
    else 
      return op;
  }

  static bool is_feature(const Halfedge_handle he) {
    return 
      he->is_feature_edge() || he->opposite()->is_feature_edge();
  }

  /** Follow a polyline or a polygon, from the halfedge he. */
  template <typename Polyline_and_context, typename Polylines_output_iterator>
  Polylines_output_iterator
  follow_half_edge(const Halfedge_handle he,
                   Polylines_output_iterator polylines_out, 
                   Polyline_and_context = Polyline_and_context())
  {
    typename Feature_edges_set::iterator it = 
      edges_to_consider.find(canonical(he));
    if(it == edges_to_consider.end()) {
      return polylines_out;
    }

    Polyline_and_context polyline;
    polyline.polyline_content.push_back(he->opposite()->vertex()->point());

    Halfedge_handle current_he = he;

    Set_of_indices set_of_indices_of_current_edge
      = edge_indices(current_he);

    do {
      CGAL_assertion(!set_of_indices_of_current_edge.empty());
      CGAL_assertion(is_feature(current_he));
      CGAL_assertion_code(const size_type n = )
        edges_to_consider.erase(canonical(current_he));
      CGAL_assertion(n > 0);
      Vertex_handle v = current_he->vertex();
      polyline.polyline_content.push_back(v->point());
      // std::cerr << v->point() << std::endl;
      if(corner_vertices.count(v) > 0) break;
      typename Polyhedron::Halfedge_around_vertex_circulator 
        loop_he = v->vertex_begin();
      ++loop_he;
      // CGAL_assertion((&*loop_he) != (&*current_he) );
      while((&*loop_he) == (&*current_he) ||
            (!is_feature(loop_he)) ) {
        ++loop_he;
        // CGAL_assertion((&*loop_he) != (&*current_he) );
      }

      Set_of_indices set_of_indices_of_next_edge = 
        edge_indices(loop_he);

      if(! (set_of_indices_of_next_edge.size() == 
            set_of_indices_of_current_edge.size() 
            &&
            std::equal(set_of_indices_of_next_edge.begin(),
                       set_of_indices_of_next_edge.end(),
                       set_of_indices_of_current_edge.begin())) ) 
      {
        // the vertex is a special vertex, a new corner
#ifdef CGAL_MESH_3_PROTECTION_DEBUG
        std::cerr << "New corner vertex " << v->point() << std::endl;
        std::cerr << "  indices were: ";
        BOOST_FOREACH(typename Set_of_indices::value_type i,
                      set_of_indices_of_current_edge) {
          std::cerr << i << " ";
        }
        std::cerr << "\n           now: ";
        BOOST_FOREACH(typename Set_of_indices::value_type i,
                      set_of_indices_of_next_edge) {
          std::cerr << i << " ";
        }
        std::cerr << "\n";
#endif
        ++v->nb_of_feature_edges;
        corner_vertices.insert(v);
        polyline.context.adjacent_patches_ids=set_of_indices_of_current_edge;
        *polylines_out++ = polyline;
        polyline.polyline_content.clear();
        polyline.polyline_content.push_back(loop_he->vertex()->point());
        set_of_indices_of_current_edge = set_of_indices_of_next_edge;
      }

      current_he = loop_he->opposite();
    } while(current_he != he );

    polyline.context.adjacent_patches_ids=set_of_indices_of_current_edge;
    *polylines_out++ = polyline;
    return polylines_out;
  }

  /** Loop around a corner vertex, and try to follow a polyline of feature
      edges, from each incident edge. */
  template <typename Polyline_and_context, typename Polylines_output_iterator>
  Polylines_output_iterator 
  loop_around_corner(const Vertex_handle v,
                     Polylines_output_iterator polylines_out, 
                     Polyline_and_context empty_polyline = 
                     Polyline_and_context() )
  {
    typename Polyhedron::Halfedge_around_vertex_circulator 
      he = v->vertex_begin(), end(he);
    do {
      CGAL_assertion(he->vertex() == v);
      polylines_out = follow_half_edge(he->opposite(), 
                                       polylines_out,
                                       empty_polyline);
      ++he;
    } while(he != end);
    return polylines_out;
  }

  /** For a non-corner vertex v (that is incident to two feature edges),
      measure the angle between the two edges, and mark the vertex as corner
      edge, if the angle is < 120°. **/
  static bool measure_angle(const Vertex_handle v)
  {
    Halfedge_handle e1;
    Halfedge_handle e2;
    typename Polyhedron::Halfedge_around_vertex_circulator he =
      v->vertex_begin(), end(he);
    // std::cerr << "measure_handle(" << (void*)(&*v)
    //           << " = " << v->point() << ")";
    bool first = true;
    bool done = false;
    do {
      CGAL_assertion(he->vertex() == v);
      // std::cerr << he->opposite()->vertex()->point() << std::endl;
      if(is_feature(he)) {
        if(first) {
          e1 = he;
          first = false;
        }
        else {
          if(done) {
            std::cerr << v->point() << " should be a corner!\n"
                      << "  Too many adjacent feature edges!\n";
            return false;
          }
          e2 = he;
          done = true;
        }
        // std::cerr << "x";
      }
      // else
        // std::cerr << ".";
      ++he;
    } while(he != end);
    if(!done) {
      std::cerr << v->point() << " should be a corner!\n"
                << "  Not enough adjacent feature edge!\n";
      return false;
    }
      // std::cerr << "\n";
    const Point_3 pv = v->point();
    const Point_3 pa = e1->opposite()->vertex()->point();
    const Point_3 pb = e2->opposite()->vertex()->point();
    const typename Geom_traits::Vector_3 av = pv - pa;
    const typename Geom_traits::Vector_3 bv = pv - pb;
    const typename Geom_traits::FT sc_prod = av * bv;
    if( sc_prod >= 0 ||
        (sc_prod < 0 && 
         CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
    {
      // std::cerr << "Corner (" << pa << ", " << pv
      //           << ", " << pb << ")\n";
      return true;
    }
    else {
      return false;
    }
  }

  template <typename Polyline_and_context, 
            typename Polylines_output_iterator>
  Polylines_output_iterator
  operator()(Polyhedron* pMesh, 
             Polylines_output_iterator out_it, 
             Polyline_and_context empty_polyline)
  {
    // That call orders the set of edges of the polyhedron, so that the
    // feature edges are at the end of the sequence of edges.
    // pMesh->normalize_border();
    Vertices_counter feature_vertices;

    // Iterate over all edges, and find out which vertices are corner
    // vertices (more than two incident feature edges).
    for(typename Polyhedron::Edge_iterator 
          eit = pMesh->edges_begin (),
          end = pMesh->edges_end();
        eit != end; ++eit)
    {
      if(!eit->is_feature_edge()) continue;
      edges_to_consider.insert(canonical(eit));
      typename Polyhedron::Vertex_handle v = eit->vertex();
      for(unsigned i = 0; i < 2; ++i) {
        if(++feature_vertices[v] == 3)
          corner_vertices.insert(v);
        v = eit->opposite()->vertex();
      }
    }

    for(typename Polyhedron::Vertex_iterator 
          vit = pMesh->vertices_begin (),
          end = pMesh->vertices_end();
        vit != end; ++vit)
    {
      if(feature_vertices.count(vit) !=0 && 
         feature_vertices[vit] == 1) {
        corner_vertices.insert(vit);
      }
    }

#ifdef CGAL_MESH_3_PROTECTION_DEBUG
    std::cerr << "Corner vertices: " << corner_vertices.size() << std::endl;
    std::cerr << "Feature vertices: " << feature_vertices.size() << std::endl;
#endif
    
    // // Iterate over non-corner feature vertices, and measure the angle.
    for(typename Vertices_counter::iterator it = feature_vertices.begin(),
        end = feature_vertices.end(); it != end; ++it)
    {
      const Vertex_handle v = it->first;
      if(corner_vertices.count(v) == 0) {
        CGAL_assertion(it->second == 2);
        if(measure_angle(v)) {
          corner_vertices.insert(v);
        }
      }
    }
#ifdef CGAL_MESH_3_PROTECTION_DEBUG
    std::cerr << "New corner vertices: "
              << corner_vertices.size() << std::endl;
#endif

    // Follow the polylines...
    for(typename Vertices_set::iterator it = corner_vertices.begin(),
        end = corner_vertices.end(); it != end; ++it)
    {
      out_it = loop_around_corner(*it, out_it, empty_polyline);
    }

    // ... and the cycles.
    while(! edges_to_consider.empty() ) {
      out_it = follow_half_edge(*edges_to_consider.begin(),
                                out_it,
                                empty_polyline);
    }

    return out_it;
  }
};

template <typename Polyhedron, 
          typename Polyline_and_context, 
          typename Polylines_output_iterator>
Polylines_output_iterator
detect_polylines(Polyhedron* pMesh, 
                 Polylines_output_iterator out_it) {

  Detect_polylines<Polyhedron> go;
  Polyline_and_context empty_polyline;
  return go(pMesh, out_it, empty_polyline);
}

} // end namespace CGAL::Mesh_3
} // end namespace CGAL


#endif // CGAL_MESH_3_DETECT_POLYLINES_IN_POLYHEDRA_H
