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
#ifndef CGAL_IO_ARRANGEMENT_2_READER_H
#define CGAL_IO_ARRANGEMENT_2_READER_H

/*! \file
 * The header file for the Arrangement_2_reader<Arrangement> class.
 */

#include <CGAL/Arr_accessor.h>
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

  typedef CGAL::Arr_accessor<Arrangement_2>               Arr_accessor;
  typedef typename Arr_accessor::Dcel_vertex              DVertex;
  typedef typename Arr_accessor::Dcel_halfedge            DHalfedge;
  typedef typename Arr_accessor::Dcel_face                DFace;
  typedef typename Arr_accessor::Dcel_hole                DHole;
  typedef typename Arr_accessor::Dcel_isolated_vertex     DIso_vert;
  
  // Data members:
  Arrangement_2&           m_arr;
  Arr_accessor             m_arr_access;
  Point_2                  m_point;
  std::vector<DVertex*>    m_vertices;
  X_monotone_curve_2       m_curve;
  std::vector<DHalfedge*>  m_halfedges;
  DVertex                 *v_bl;
  DVertex                 *v_tl;
  DVertex                 *v_br;
  DVertex                 *v_tr;

private:

  // Copy constructor and assignment operator - not supported.
  Arrangement_2_reader (const Self& );
  Self& operator= (const Self& );

public:

  /*! Constructor. */
  Arrangement_2_reader (Arrangement_2& arr) :
    m_arr (arr),
    m_arr_access (arr),
    v_bl (NULL), v_tl (NULL), v_br (NULL), v_tr (NULL)
  {}

  /*! Destructor. */
  virtual ~Arrangement_2_reader ()
  {}

  /*! Read the arrangement. */
  template <class Formatter>
  void operator()(Formatter& formatter)
  {
    // Clear the exisiting arrangement so it contains no DCEL features.
    m_arr_access.clear_all();

    // Read the arrangement dimensions.
    formatter.read_arrangement_begin();

    const Size  number_of_vertices = formatter.read_size("number_of_vertices");
    const Size  number_of_halfedges = 2*formatter.read_size("number_of_edges");
    const Size  number_of_faces = formatter.read_size("number_of_faces");
    Size        k;
    
    // Create the four fictitious DCEL vertices.
    v_bl =  m_arr_access.new_vertex_at_infinity (MINUS_INFINITY,
                                                 MINUS_INFINITY);
    v_tl =  m_arr_access.new_vertex_at_infinity (MINUS_INFINITY,
                                                 PLUS_INFINITY);
    v_br =  m_arr_access.new_vertex_at_infinity (PLUS_INFINITY,
                                                 MINUS_INFINITY);
    v_tr =  m_arr_access.new_vertex_at_infinity (PLUS_INFINITY,
                                                 PLUS_INFINITY);

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

    formatter.read_arrangement_end();
    return;
  }

protected:

  /*! Read a DCEL vertex. */
  template <class Formatter>
  DVertex* _read_vertex (Formatter& formatter)
  {
    formatter.read_vertex_begin();

    // Read the infinity types.
    Boundary_type   inf_x = Boundary_type (formatter.read_vertex_index());
    Boundary_type   inf_y = Boundary_type (formatter.read_vertex_index());
    DVertex        *new_v;

    if (inf_x == NO_BOUNDARY && inf_y == NO_BOUNDARY)
    {
      // Read the point associated with the vertex.
      formatter.read_point (m_point);

      // Allocate a new DCEL vertex and associate it with this point.
      new_v = m_arr_access.new_vertex (m_point);

      // Read any auxiliary data associated with the vertex.
      formatter.read_vertex_data (Vertex_handle (new_v));
    }
    else
    {
      // Allocate a vertex at infinity.
      new_v = m_arr_access.new_vertex_at_infinity (inf_x, inf_y);
    }

    formatter.read_vertex_end();
    return (new_v);
  }
 
  /*! Read a DCEL edge (a pair of twin halfedges). */
  template <class Formatter>
  DHalfedge* _read_edge (Formatter& formatter)
  {
    formatter.read_edge_begin();

    // Read the indices of the end-vertices and the edge direction.
    int                 source_idx = formatter.read_vertex_index();
    int                 target_idx = formatter.read_vertex_index();
    int                 direction = formatter.read_vertex_index();
    DHalfedge          *new_he;
    DVertex            *src_v;
    DVertex            *trg_v;

    if (source_idx == -1)
      src_v = v_bl;
    else if (source_idx == -2)
      src_v = v_tl;
    else if (source_idx == -3)
      src_v = v_br;
    else if (source_idx == -4)
      src_v = v_tr;
    else
      src_v = m_vertices[source_idx];

    if (target_idx == -1)
      trg_v = v_bl;
    else if (target_idx == -2)
      trg_v = v_tl;
    else if (target_idx == -3)
      trg_v = v_br;
    else if (target_idx == -4)
      trg_v = v_tr;
    else
      trg_v = m_vertices[target_idx];

    if (source_idx >= 0 || target_idx >= 0)
    {
      // Read the x-monotone curve associated with the edge. 
      formatter.read_x_monotone_curve (m_curve);

      // Allocate a pair of new DCEL halfegdes and associate them with the
      // x-monotone curve we read.
      new_he = m_arr_access.new_edge (m_curve);
    }
    else
    {
      // Allocate a new fictitious edge.
      new_he = m_arr_access.new_fictitious_edge();
    }

    // Set the cross pointers between the twin halfedges and the end vertices.
    trg_v->set_halfedge (new_he);
    new_he->set_vertex (trg_v);
    
    src_v->set_halfedge (new_he->opposite());
    new_he->opposite()->set_vertex (src_v);
   
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
    if (source_idx >= 0 || target_idx >= 0)
    {
      formatter.read_halfedge_data (Halfedge_handle (new_he));
      formatter.read_halfedge_data (Halfedge_handle ((new_he->opposite())));
    }

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
    const Size  outer_size = formatter.read_size ("halfedges_on_outer_CCB");
    DFace      *new_f = NULL;
    DHalfedge  *he;

    if (outer_size == 0)
    {
      // Allocate the fictitious DCEL face.
      new_f = m_arr_access.new_face();
      new_f->set_halfedge (NULL);
    }
    else
    {
      // Allocate a new DCEL face and read its outer CCB.
      new_f = m_arr_access.new_face();
      he = _read_ccb (formatter, new_f, outer_size, NULL);
      new_f->set_halfedge (he);
    }
    formatter.read_outer_ccb_end();

    // Read the holes inside the face.
    formatter.read_holes_begin();

    DHole      *new_hole;
    const Size  n_holes = formatter.read_size ("number_of_holes");
    Size        inner_size;
    Size        k;

    for (k = 0; k < n_holes; k++)
    {
      // Allocate a new hole record and set its incident face.
      new_hole = m_arr_access.new_hole();
      new_hole->set_face (new_f);

      // Read the current hole.
      inner_size = formatter.read_size ("halfedges_on_inner_CCB");
      he = _read_ccb (formatter, new_f, inner_size, new_hole);
      new_hole->set_iterator (new_f->add_hole (he));
    }
    formatter.read_holes_end();

    // Read the isolated vertices inside the face.
    formatter.read_isolated_vertices_begin();

    DIso_vert   *new_iso_vert;
    Size         n_isolated_vertices = 
                          formatter.read_size ("number_of_isolated_vertices");
    std::size_t  v_idx;
    DVertex*     iso_v;

    for (k = 0; k < n_isolated_vertices; k++)
    {
      // Allocate a new isolated vertex record and set its incident face.
      new_iso_vert = m_arr_access.new_isolated_vertex();
      new_iso_vert->set_face (new_f);

      // Read the current isolated vertex.
      v_idx = formatter.read_vertex_index ();
      iso_v = m_vertices[v_idx];
      iso_v->set_isolated_vertex (new_iso_vert);
      new_iso_vert->set_iterator (new_f->add_isolated_vertex (iso_v));
    }
    formatter.read_isolated_vertices_end();

    // Read any auxiliary data associated with the face.
    if (outer_size != 0)
      formatter.read_face_data (Face_handle (new_f));
    
    formatter.read_face_end();

    return;
  }

  /*!
   * Read a circular boundary of a conncted component.
   * \param formatter The formatter.
   * \param f The incident DCEL face.
   * \param boundary_size The number of halfedges along the boundary.
   * \param p_hole If NULL, the CCB corresponds to the outer boundary of f;
   *               otherwise, it corresponds to an inner component (hole).
   * \return A pointer to the first halfedge read.
   */
  template <class Formatter>
  DHalfedge* _read_ccb (Formatter& formatter, 
                        DFace *f,
                        Size boundary_size,
                        DHole *p_hole)
  {
    formatter.read_ccb_halfedges_begin();
 
    // Find the first halfedge, and set its incident face.
    std::size_t   first_idx = formatter.read_halfedge_index();
    DHalfedge    *first_he = m_halfedges [first_idx];

    if (p_hole == NULL)
      first_he->set_face (f);
    else
      first_he->set_hole (p_hole);

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
      if (p_hole == NULL)
        curr_he->set_face (f);
      else
        curr_he->set_hole (p_hole);

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
