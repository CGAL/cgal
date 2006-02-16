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
//                 (based on old version by Ester Ezra)
#ifndef CGAL_IO_ARRANGEMENT_2_WRITER_H
#define CGAL_IO_ARRANGEMENT_2_WRITER_H

/*! \file
 * The header file for the Arrangement_2_writer<Arrangement> class.
 */

#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <algorithm>

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
  
  typedef typename Arrangement_2::Vertex_const_iterator
                                                      Vertex_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator
                                                      Halfedge_const_iterator;
  typedef typename Arrangement_2::Edge_const_iterator
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
  
  typedef Inverse_index<Halfedge_const_iterator>      Halfedges_index;
  typedef Inverse_index<Vertex_const_iterator>        Vertices_index;

  // Data memebrs:
  const Arrangement_2&   m_arr;
  Halfedges_index*       m_he_index;
  Vertices_index*        m_v_index;

private:

  // Copy constructor and assignment operator - not supported.
  Arrangement_2_writer (const Self& );
  Self& operator= (const Self& );

public:

  /*! Constructor. */
  Arrangement_2_writer (const Arrangement_2& arr) :
    m_arr(arr)
  {
    m_he_index = new Halfedges_index (arr.halfedges_begin(),
                                      arr.halfedges_end());
    m_v_index = new Vertices_index (arr.vertices_begin(),
                                    arr.vertices_end());
  }

  /*! Destructor. */
  virtual ~Arrangement_2_writer()
  {
    if (m_he_index)
      delete m_he_index;
    if (m_v_index)
      delete m_v_index;
  }

  /*! Write the arrangement. */
  template <class Formatter>
  void operator() (Formatter& formatter) const
  {
    formatter.write_arrangement_begin();
    formatter.write_size ("number_of_vertices", m_arr.number_of_vertices());
    formatter.write_size ("number_of_edges", m_arr.number_of_edges());
    formatter.write_size ("number_of_faces", m_arr.number_of_faces());

    // Write the vertices.
    formatter.write_vertices_begin();
    Vertex_const_iterator  vit;
    for (vit = m_arr.vertices_begin(); vit != m_arr.vertices_end(); ++vit)
      _write_vertex (formatter, vit);
    formatter.write_vertices_end();

    // Write the edges.
    formatter.write_edges_begin();
    Edge_const_iterator    eit;
    for (eit = m_arr.edges_begin(); eit != m_arr.edges_end(); ++eit)
      _write_edge (formatter, eit);
    formatter.write_edges_end();

    // Write the faces.
    formatter.write_faces_begin();
    Face_const_iterator    fit;
    for (fit = m_arr.faces_begin(); fit != m_arr.faces_end(); ++fit)
      _write_face(formatter, fit);
    formatter.write_faces_end();

    formatter.write_arrangement_end();
  }

protected:

  /*! Write a vertex. */
  template <class Formatter>
  void _write_vertex (Formatter& formatter, Vertex_const_iterator vit) const
  {
    Vertex_const_handle  v = vit;
    formatter.write_vertex_begin();
    formatter.write_point (v->point()); // Write the associated point.
    formatter.write_vertex_data (v);    // Write additional user-defined data.
    formatter.write_vertex_end();

    return;
  }

  /*! Write an edge (a pair of halfedges). */
  template <class Formatter>
  void _write_edge (Formatter& formatter, Edge_const_iterator hit) const
  {
    Halfedge_const_handle he = hit;

    formatter.write_edge_begin ();
    formatter.write_vertex_index (_get_index(he->source()));
    formatter.write_vertex_index (_get_index(he->target()));
    
    if (he->direction() == SMALLER)
      formatter.write_vertex_index (0);
    else
      formatter.write_vertex_index (1);
      
    formatter.write_x_monotone_curve (he->curve()); 
                                         // Write the associated curve.
    formatter.write_halfedge_data (he);  // Write additional user-defined data.
    formatter.write_halfedge_data (he->twin());
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
    if (f->is_unbounded())
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
  std::size_t _get_index (Vertex_const_handle v) const
  {
    return (*m_v_index)[v];
  }

  /*! Get the mapped index of a given halfegde. */
  std::size_t _get_index (Halfedge_const_handle he) const
  {
    return (*m_he_index)[he];
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
