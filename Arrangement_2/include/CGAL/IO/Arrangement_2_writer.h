// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Michal Meyerovitch <gorgymic@post.tau.ac.il>
//                 Ron Wein           <wein@post.tau.ac.il>
//                 Efi Fogel          <efif@post.tau.ac.il>
//                 (based on old version by Ester Ezra)

#ifndef CGAL_IO_ARRANGEMENT_2_WRITER_H
#define CGAL_IO_ARRANGEMENT_2_WRITER_H

/*! \file
 * The header file for the Arrangement_2_writer<Arrangement> class.
 */

#include <CGAL/Arr_accessor.h>
#include <map>

CGAL_BEGIN_NAMESPACE

/*! \class
 * An auxiliary class for writing an arrangement to an output stream.
 */
template <class Arrangement_>
class Arrangement_2_writer
{
public:

  typedef Arrangement_                                    Arrangement_2;
  typedef Arrangement_2_writer<Arrangement_2>             Self;

protected:

  typedef typename Arrangement_2::Size                    Size;
  
  typedef CGAL::Arr_accessor<Arrangement_2>            Arr_accessor;
  typedef typename Arr_accessor::All_vertex_const_iterator
                                                       Vertex_const_iterator;
  typedef typename Arr_accessor::All_edge_const_iterator
                                                       Edge_const_iterator;
  typedef typename Arrangement_2::Face_const_iterator 
                                                       Face_const_iterator;

  typedef typename Arrangement_2::Vertex_const_handle 
                                                      Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle
                                                      Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle   
                                                      Face_const_handle;
 
  typedef typename Arrangement_2::Hole_const_iterator
                                             Hole_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                             Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Isolated_vertex_const_iterator
                                             Isolated_vertex_const_iterator;
  
  typedef typename Arrangement_2::Vertex              Vertex;
  typedef typename Arrangement_2::Halfedge            Halfedge;
  typedef std::map<const Vertex*, int>                Vertex_index_map;
  typedef std::map<const Halfedge*, int>              Halfedge_index_map;

  // Data memebrs:
  const Arrangement_2&   m_arr;
  const Arr_accessor     m_arr_access;
  int                    m_curr_v;
  Vertex_index_map       m_v_index;
  int                    m_curr_he;
  Halfedge_index_map     m_he_index;

private:

  // Copy constructor and assignment operator - not supported.
  Arrangement_2_writer (const Self& );
  Self& operator= (const Self& );

public:

  /*! Constructor. */
  Arrangement_2_writer (const Arrangement_2& arr) :
    m_arr (arr),
    m_arr_access (const_cast<Arrangement_2&>(arr)),
    m_curr_v (0),
    m_curr_he (0)
  {}

  /*! Destructor. */
  virtual ~Arrangement_2_writer()
  {}

  /*! Write the arrangement. */
  template <class Formatter>
  void operator() (Formatter& formatter)
  {
    formatter.write_arrangement_begin();
    formatter.write_size ("number_of_vertices",
                          m_arr.number_of_vertices() +
                          m_arr.number_of_vertices_at_infinity());
    formatter.write_size ("number_of_edges",
                          m_arr.number_of_edges() +
                          m_arr.number_of_vertices_at_infinity() + 4);
    formatter.write_size ("number_of_faces", m_arr.number_of_faces() + 1);

    // Reset indices.
    m_curr_v = 0;
    m_curr_he = 0;

    // Write the vertices.
    formatter.write_vertices_begin();
    Vertex_const_iterator  vit;
    for (vit = m_arr_access.all_vertices_begin(); 
         vit != m_arr_access.all_vertices_end(); ++vit)
    {
      _write_vertex (formatter, vit);
    }
    formatter.write_vertices_end();

    // Write the edges.
    formatter.write_edges_begin();
    Edge_const_iterator    eit;
    for (eit = m_arr_access.all_edges_begin();
         eit != m_arr_access.all_edges_end(); ++eit)
    {
      _write_edge (formatter, eit);
    }
    formatter.write_edges_end();

    // Write the faces (the fictitious face first).
    formatter.write_faces_begin();

    _write_face(formatter, m_arr_access.fictitious_face());

    Face_const_iterator    fit;
    for (fit = m_arr.faces_begin(); fit != m_arr.faces_end(); ++fit)
      _write_face(formatter, fit);
    formatter.write_faces_end();

    formatter.write_arrangement_end();
  }

protected:

  /*! Write a vertex. */
  template <class Formatter>
  void _write_vertex (Formatter& formatter, Vertex_const_iterator vit)
  {
    // Map the current vertex to its index.
    Vertex_const_handle  v = vit;

    m_v_index[&(*vit)] = m_curr_v;
    ++m_curr_v;

    // Write the vertex.
    formatter.write_vertex_begin();
    formatter.write_vertex_index (static_cast<int> (v->boundary_in_x()));
    formatter.write_vertex_index (static_cast<int> (v->boundary_in_y()));
    if (! v->is_at_infinity())
      formatter.write_point (v->point()); // Write the associated point.
    formatter.write_vertex_data (v);    // Write additional user-defined data.
    formatter.write_vertex_end();

    return;
  }

  /*! Write an edge (a pair of halfedges). */
  template <class Formatter>
  void _write_edge (Formatter& formatter, Edge_const_iterator hit)
  {
    // Map the halfedge and its twin to their indices.
    Halfedge_const_handle  he = hit;
    Halfedge_const_handle  he_t = he->twin();

    m_he_index[&(*he)] = m_curr_he;
    ++m_curr_he;
    m_he_index[&(*he_t)] = m_curr_he;
    ++m_curr_he;

    // Write the edge.
    formatter.write_edge_begin ();
    formatter.write_vertex_index (_get_index(he->source()));
    formatter.write_vertex_index (_get_index(he->target()));
    
    if (he->direction() == SMALLER)
      formatter.write_vertex_index (0);
    else
      formatter.write_vertex_index (1);
      
    if (! he->is_fictitious())
    {
      // Write the associated curve.
      formatter.write_x_monotone_curve (he->curve()); 

      // Write additional user-defined data.
      formatter.write_halfedge_data (he);
      formatter.write_halfedge_data (he_t);
    }
    formatter.write_edge_end ();

    return;
  }

  /*! Write a face. */
  template <class Formatter>
  void _write_face (Formatter& formatter, Face_const_iterator fit) const
  {
    Face_const_handle     f = fit;

    formatter.write_face_begin();

    // Write the outer CCB.
    formatter.write_outer_ccb_begin();
    if (f->is_fictitious())
    {
      formatter.write_size ("halfedges_on_outer_CCB", 0);
    }
    else
    {
      Ccb_halfedge_const_circulator  out_ccb = f->outer_ccb();
      const std::size_t              n = _circulator_size (out_ccb);

      formatter.write_size ("halfedges_on_outer_CCB", n);
      _write_ccb (formatter, out_ccb);
    }
    formatter.write_outer_ccb_end();

    // Write the holes inside the face.
    formatter.write_holes_begin();
    const std::size_t    n_holes = std::distance (f->holes_begin(),
                                                  f->holes_end());
    formatter.write_size ("number_of_holes", n_holes);

    Hole_const_iterator  hole_it;
    for (hole_it = f->holes_begin(); hole_it != f->holes_end(); ++hole_it)
    {
      Ccb_halfedge_const_circulator  in_ccb = (*hole_it);      
      const std::size_t              n = _circulator_size (in_ccb);

      formatter.write_size ("halfedges_on_inner_CCB", n);
      _write_ccb (formatter, in_ccb);      
    }
    formatter.write_holes_end();

    // Write the isolated vertices inside the face.
    formatter.write_isolated_vertices_begin();
    std::size_t  n_isolated = std::distance (f->isolated_vertices_begin(),
                                             f->isolated_vertices_end());
    formatter.write_size ("number_of_isolated_vertices", n_isolated);

    Isolated_vertex_const_iterator  iso_vit;
    for (iso_vit = f->isolated_vertices_begin();
         iso_vit != f->isolated_vertices_end(); ++iso_vit)
    {
      formatter.write_vertex_index (_get_index(iso_vit));
    }
    formatter.write_isolated_vertices_end();
    
    // Write additional user-defined data associated with the face.
    if (! f->is_fictitious())
      formatter.write_face_data (f);

    formatter.write_face_end();

    return;
  }

  /*! Write the edges along a given CCB. */
  template <class Formatter>   
  void _write_ccb (Formatter& formatter,
                   Ccb_halfedge_const_circulator circ) const      
  {
    Ccb_halfedge_const_circulator  curr = circ;

    formatter.write_ccb_halfedges_begin();
    do {
      formatter.write_halfedge_index (_get_index(curr));
      ++curr;
    } while (curr != circ);
    formatter.write_ccb_halfedges_end();

    return;
  }
  
  /*! Get the mapped index of a given vertex. */
  int _get_index (Vertex_const_handle v) const
  {
    if (v == m_arr_access.bottom_left_fictitious_vertex())
      return (-1);
    else if (v == m_arr_access.top_left_fictitious_vertex())
      return (-2);
    else if (v == m_arr_access.bottom_right_fictitious_vertex())
      return (-3);
    else if (v == m_arr_access.top_right_fictitious_vertex())
      return (-4);
    
    typename Vertex_index_map::const_iterator  pos = m_v_index.find (&(*v));

    CGAL_assertion (pos != m_v_index.end());
    return (pos->second);
  }

  /*! Get the mapped index of a given halfegde. */
  int _get_index (Halfedge_const_handle he) const
  {
    typename Halfedge_index_map::const_iterator
                                               pos = m_he_index.find (&(*he));

    CGAL_assertion (pos != m_he_index.end());
    return (pos->second);
  }

  /*! Get the number of edges along a given CCB. */
  std::size_t _circulator_size (Ccb_halfedge_const_circulator circ) const
  {
    CGAL_assertion (circ != CGAL_CIRC_NULL);

    std::size_t                    n = 0;
    Ccb_halfedge_const_circulator  curr = circ;

    do {
      ++n;
      ++curr;
    } while(curr != circ);

    return (n);
  }  
};

CGAL_END_NAMESPACE

#endif // CGAL_IO_ARRANGEMENT_2_WRITER_H 
