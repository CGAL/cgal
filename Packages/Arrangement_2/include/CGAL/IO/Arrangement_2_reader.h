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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michal Meyerovitch <gorgymic@post.tau.ac.il>
//                 Ron Wein           <wein@post.tau.ac.il>
//                 (based on old version by Ester Ezra)
#ifndef CGAL_IO_ARRANGEMENT_2_READER_H
#define CGAL_IO_ARRANGEMENT_2_READER_H

/*! \file
 * The header file for the Arrangement_2_reader<Arrangement> class.
 */

#include <CGAL/IO/Arrangement_2_ascii_formatter.h>
#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <algorithm>
#include <iostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * An auxiliary class for reading an arrangement from an input stream.
 */
template <class Arrangement_>
class Arrangement_2_reader
{
public:

  typedef Arrangement_                                    Arrangement_2;
  typedef Arrangement_2_reader<Arrangement_2>             Self;

protected:
 
  typedef typename Arrangement_2::Size                    Size;
  typedef typename Arrangement_2::Vertex_iterator         Vertex_iterator;
  typedef typename Arrangement_2::Halfedge_iterator       Halfedge_iterator;
  typedef typename Arrangement_2::Face_iterator           Face_iterator;

  typedef typename Arrangement_2::Vertex_handle           Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement_2::Face_handle             Face_handle;

  typedef typename Arrangement_2::Dcel                    Dcel;  

  typedef typename Arrangement_2::Traits_2                Traits_2;
  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits_2::Point_2                      Point_2;

protected:

  typedef typename Dcel::Vertex                           DVertex;
  typedef typename Dcel::Halfedge                         DHalfedge;
  typedef typename Dcel::Face                             DFace;

  typedef typename Arrangement_2::Stored_point_2          Stored_point_2;
  typedef typename Arrangement_2::Stored_curve_2          Stored_curve_2;

  // Data members:
  Arrangement_2&           m_arr;
  std::vector<DVertex*>    m_vertices;
  std::vector<DHalfedge*>  m_halfedges;

private:

  // Copy constructor and assignment operator - not supported.
  Arrangement_2_reader (const Self& );
  Self& operator= (const Self& );

public:

  /*! Constructor. */
  Arrangement_2_reader (Arrangement_2& arr) :
    m_arr(arr)
  {}

  /*! Destructor. */
  virtual ~Arrangement_2_reader ()
  {}

  /*! Read the arrangement. */
  template <class Formatter>
  void operator()(Formatter& formatter)
  {
    // Clear the exisiting arrangement so it contains only an unbounded face.
    m_arr.clear();

    // Read the arrangement dimensions.
    formatter.read_arr_begin();

    const Size  number_of_vertices = formatter.read_size("number_of_vertices");
    const Size  number_of_halfedges = 2*formatter.read_size("number_of_edges");
    const Size  number_of_faces = formatter.read_size("number_of_faces");
    Size        k;
    
    // Read the DCEL vertices and store them in the vertices vector.
    formatter.read_vertices_begin();

    m_vertices.resize (number_of_vertices);
    for (k = 0; k < number_of_vertices; k++)
      m_vertices[k] = _read_vertex (formatter);

    formatter.read_vertices_end();

    // Read the DCEL halfedges and store them in the halfedges vector.
    DHalfedge   *he = NULL;

    formatter.read_edges_begin();

    m_halfedges.resize (number_of_halfedges);
    for (k = 0; k < number_of_halfedges; k += 2)
    { 
      he = _read_edge (formatter);
      m_halfedges[k] = he;
      m_halfedges[k + 1] = he->opposite();
    }
    formatter.read_edges_end();

    // Read the DCEL faces.
    formatter.read_faces_begin();
    for (k = 0; k < number_of_faces; k++)
      _read_face (formatter);
    formatter.read_faces_end();

    formatter.read_arr_end();
    return;
  }

protected:

  /*! Read a DCEL vertex. */
  template <class Formatter>
  DVertex* _read_vertex (Formatter& formatter)
  {
    formatter.read_vertex_begin();

    // Allocate a new DCEL vertex and associate it with the point we read.
    DVertex*        new_v = m_arr.dcel.new_vertex();
    Point_2         p;

    formatter.read_point (p);

    Stored_point_2 *p_p = m_arr._new_point (p);
    new_v->set_point(p_p);

    // Read any auxiliary data associated with the vertex.
    formatter.read_vertex_data (Vertex_handle (new_v));

    formatter.read_vertex_end();
    return (new_v);
  }
 
  /*! Read a DCEL edge (a pair of twin halfedges). */
  template <class Formatter>
  DHalfedge* _read_edge (Formatter& formatter)
  {
    formatter.read_edge_begin();

    // Allocate a pair of new DCEL halfegdes and associate them with the
    // x-monotone curve we read.
    DHalfedge          *new_he = m_arr.dcel.new_edge();
    std::size_t         source_idx = formatter.read_vertex_index();
    std::size_t         target_idx = formatter.read_vertex_index();
    std::size_t         direction = formatter.read_vertex_index();
    X_monotone_curve_2  cv;

    formatter.read_curve (cv);

    Stored_curve_2     *p_cv = m_arr._new_curve (cv);
    new_he->set_curve(p_cv);

    // Set the cross pointers between the twin halfedges and the end vertices.
    m_vertices[target_idx]->set_halfedge (new_he);
    new_he->set_vertex (m_vertices[target_idx]);
    
    m_vertices[source_idx]->set_halfedge (new_he->opposite());
    new_he->opposite()->set_vertex (m_vertices[source_idx]);
   
    // Set the directionf of the halfedges.
    if (direction == 0)
    {
      new_he->set_direction (SMALLER);
    }
    else
    {
      CGAL_assertion (direction == 1);
      new_he->set_direction (LARGER);
    }

    // Read any auxiliary data associated with the halfedges.
    formatter.read_halfedge_data (Halfedge_handle (new_he));
    formatter.read_halfedge_data (Halfedge_handle ((new_he->opposite())));

    formatter.read_edge_end();
    return (new_he);
  }

  /*! Read a DCEL face. */
  template <class Formatter>
  void _read_face(Formatter& formatter)
  {
    formatter.read_face_begin();

    // Try reading the outer CCB of the face.
    formatter.read_outer_ccb_begin();
    Size       outer_size = formatter.read_size ("halfedges_on_outer_CCB");
    DFace     *new_f = NULL;
    DHalfedge *he;

    if (outer_size == 0)
    {
      // We currently read the unbounded face.
      new_f = m_arr.un_face;
    }
    else
    {
      // Allocate a new DCEL face and read its outer CCB.
      new_f = m_arr.dcel.new_face();
      he = _read_ccb (formatter, new_f, outer_size);
      new_f->set_halfedge (he);
    }
    formatter.read_outer_ccb_end();

    // Read the holes inside the face.
    formatter.read_holes_begin();

    const Size  n_holes = formatter.read_size ("number_of_holes");
    Size        inner_size;
    Size        k;

    for (k = 0; k < n_holes; k++)
    {
      // Read the current hole.
      inner_size = formatter.read_size ("halfedges_on_inner_CCB");
      he = _read_ccb (formatter, new_f, inner_size);
      new_f->add_hole (he);
    }
    formatter.read_holes_end();

    // Read the isolated vertices inside the face.
    formatter.read_isolated_vertices_begin();

    Size         n_isolated_vertices = 
                          formatter.read_size ("number_of_isolated_vertices");
    std::size_t  v_idx;
    DVertex*     iso_v;

    for (k = 0; k < n_isolated_vertices; k++)
    {
      // Read the current isolated vertices.
      v_idx = formatter.read_vertex_index ();
      iso_v = m_vertices[v_idx];
      iso_v->set_face (new_f);
      new_f->add_isolated_vertex (iso_v);
    }
    formatter.read_isolated_vertices_end();

    m_arr.n_iso_verts += n_isolated_vertices;

    // Read any auxiliary data associated with the face.
    formatter.read_face_data (Face_handle (new_f));
    formatter.read_face_end();

    return;
  }

  /*!
   * Read a circular boundary of a conncted component.
   * \param formatter The formatter.
   * \param f The incident DCEL face.
   * \param boundary_size The number of halfedges along the boundary.
   * \return A pointer to the first halfedge read.
   */
  template <class Formatter>
  DHalfedge* _read_ccb (Formatter& formatter, 
                        DFace *f,
                        Size boundary_size)
  {
    formatter.read_ccb_halfedges_begin();
 
    // Find the first halfedge, and set its incident face.
    std::size_t   first_idx = formatter.read_halfedge_index();
    DHalfedge    *first_he = m_halfedges [first_idx];

    first_he->set_face (f);

    // Read the rest of the halfedge along the boundary.
    std::size_t   curr_idx;
    DHalfedge    *prev_he = first_he;    
    DHalfedge    *curr_he;
    Size          k;

    for (k = 1; k < boundary_size; k++)
    {
      curr_idx = formatter.read_halfedge_index();
      curr_he = m_halfedges[curr_idx];

      // Connect the previous halfedge and the current one.
      prev_he->set_next (curr_he);

      // Set the incident face.
      curr_he->set_face (f);

      prev_he = curr_he;
    }

    // Close the circular list be connecting the first and the last halfedges.
    prev_he->set_next (first_he);

    formatter.read_ccb_halfedges_end();

    // Return the first halfedge.
    return (first_he);
  }
   
};

CGAL_END_NAMESPACE

#endif // CGAL_IO_ARRANGEMENT_2_READER_H 
