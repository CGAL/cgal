// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>
//                 (based on old version by Michal Meyerovitch and Ester Ezra)
//
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

namespace CGAL {

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
    typedef typename Arrangement_2::Dcel                    Dcel;  

    typedef typename Arrangement_2::X_monotone_curve_2      X_monotone_curve_2;
    typedef typename Arrangement_2::Point_2                 Point_2;

    typedef typename Arrangement_2::Vertex_handle           Vertex_handle;
    typedef typename Arrangement_2::Halfedge_handle         Halfedge_handle;
    typedef typename Arrangement_2::Face_handle             Face_handle;

    typedef CGAL::Arr_accessor<Arrangement_2>               Arr_accessor;
    typedef typename Arr_accessor::Dcel_vertex              DVertex;
    typedef typename Arr_accessor::Dcel_halfedge            DHalfedge;
    typedef typename Arr_accessor::Dcel_face                DFace;
    typedef typename Arr_accessor::Dcel_outer_ccb           DOuter_ccb;
    typedef typename Arr_accessor::Dcel_inner_ccb           DInner_ccb;
    typedef typename Arr_accessor::Dcel_isolated_vertex     DIso_vert;
  
    // Data members:
    Arrangement_2&           m_arr;
    Arr_accessor             m_arr_access;
    Point_2                  m_point;
    std::vector<DVertex*>    m_vertices;
    X_monotone_curve_2       m_curve;
    std::vector<DHalfedge*>  m_halfedges;

  private:

    // Copy constructor and assignment operator - not supported.
    Arrangement_2_reader(const Self&);
    Self& operator=(const Self&);

  public:

    /*! Constructor. */
    Arrangement_2_reader(Arrangement_2& arr) :
      m_arr(arr),
      m_arr_access(arr)
    {}

    /*! Destructor. */
    virtual ~Arrangement_2_reader()
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

      // std::cout << number_of_vertices << std::endl;
      // std::cout << number_of_halfedges << std::endl;
      // std::cout << number_of_faces << std::endl;
    
      // Read the DCEL vertices and store them in the vertices vector.
      formatter.read_vertices_begin();

      m_vertices.resize(number_of_vertices);
      for (k = 0; k < number_of_vertices; k++)
        m_vertices[k] = _read_vertex(formatter);

      formatter.read_vertices_end();

      // Read the DCEL halfedges and store them in the halfedges vector.
      DHalfedge* he = NULL;
      formatter.read_edges_begin();

      m_halfedges.resize(number_of_halfedges);
      for (k = 0; k < number_of_halfedges; k += 2)
      { 
        he = _read_edge(formatter);
        m_halfedges[k] = he;
        m_halfedges[k + 1] = he->opposite();
      }
      formatter.read_edges_end();

      // Read the DCEL faces.
      formatter.read_faces_begin();
      for (k = 0; k < number_of_faces; k++)
        _read_face(formatter);
      formatter.read_faces_end();

      formatter.read_arrangement_end();

      // Use the accessor an update the topology-traits properties with the
      // new DCEL we have just read.
      m_arr_access.dcel_updated();
    }

  protected:

    /*! Read a DCEL vertex. */
    template <class Formatter>
    DVertex* _read_vertex(Formatter& formatter)
    {
      formatter.read_vertex_begin();

      // Read the boundary conditions.
      Arr_parameter_space   ps_x =
        Arr_parameter_space(formatter.read_vertex_index());
      Arr_parameter_space   ps_y =
        Arr_parameter_space(formatter.read_vertex_index());
      int       has_point = formatter.read_vertex_index();
      DVertex*  new_v;

      if (has_point)
      {
        // Read the point associated with the vertex.
        formatter.read_point(m_point);

        // Allocate a new DCEL vertex and associate it with this point.
        new_v = m_arr_access.new_vertex(&m_point, ps_x, ps_y);

        // Read any auxiliary data associated with the vertex.
        formatter.read_vertex_data(Vertex_handle(new_v));
      }
      else
      {
        // Allocate a vertex at infinity.
        new_v = m_arr_access.new_vertex(NULL, ps_x, ps_y);
      }

      formatter.read_vertex_end();
      return (new_v);
    }
 
    /*! Read a DCEL edge (a pair of twin halfedges). */
    template <class Formatter>
    DHalfedge* _read_edge(Formatter& formatter)
    {
      formatter.read_edge_begin();

      // Read the indices of the end-vertices and the edge direction.
      int         source_idx = formatter.read_vertex_index();
      int         target_idx = formatter.read_vertex_index();
      int         direction = formatter.read_vertex_index();
      int         has_curve = formatter.read_vertex_index();
      DHalfedge*  new_he;
      DVertex*    src_v = m_vertices[source_idx];
      DVertex*    trg_v = m_vertices[target_idx];

      if (has_curve)
      {
        // Read the x-monotone curve associated with the edge. 
        formatter.read_x_monotone_curve(m_curve);

        // Allocate a pair of new DCEL halfegdes and associate them with the
        // x-monotone curve we read.
        new_he = m_arr_access.new_edge(&m_curve);
      }
      else
      {
        // Allocate a new fictitious edge.
        new_he = m_arr_access.new_edge(NULL);
      }

      // Set the cross pointers between the twin halfedges and the end vertices.
      trg_v->set_halfedge(new_he);
      new_he->set_vertex(trg_v);
    
      src_v->set_halfedge(new_he->opposite());
      new_he->opposite()->set_vertex(src_v);
   
      // Set the direction of the halfedges.
      if (direction == 0)
      {
        new_he->set_direction(ARR_LEFT_TO_RIGHT);
      }
      else
      {
        CGAL_assertion(direction == 1);
        new_he->set_direction(ARR_RIGHT_TO_LEFT);
      }

      // Read any auxiliary data associated with the halfedges.
      if (has_curve)
      {
        formatter.read_halfedge_data(Halfedge_handle(new_he));
        formatter.read_halfedge_data(Halfedge_handle((new_he->opposite())));
      }

      formatter.read_edge_end();
      return (new_he);
    }

    /*! Read a DCEL face. */
    template <class Formatter>
    void _read_face(Formatter& formatter)
    {
      formatter.read_face_begin();

      // Allocate a new face and determine whether it is unbounded and wether it
      // is valid (non-fictitious).
      DFace*      new_f = m_arr_access.new_face();
      const bool  is_unbounded = (formatter.read_vertex_index() != 0);
      const bool  is_valid = (formatter.read_vertex_index() != 0);

      new_f->set_unbounded(is_unbounded);
      new_f->set_fictitious(! is_valid);

      // Read the outer CCBs of the face.
      formatter.read_outer_ccbs_begin();

      DOuter_ccb* new_occb;
      const Size  n_occbs = formatter.read_size("number_of_outer_ccbs");
      DHalfedge*  he;
      Size        n, k;

      for (k = 0; k < n_occbs; k++)
      {
        // Allocate a new outer CCB record and set its incident face.
        new_occb = m_arr_access.new_outer_ccb();
        new_occb->set_face(new_f);

        // Read the current outer CCB.
        n = formatter.read_size("halfedges_on_outer_ccb");
        he = _read_ccb(formatter, n, new_occb, NULL);
        new_f->add_outer_ccb(new_occb, he);
      }
      formatter.read_outer_ccbs_end();

      // Read the inner CCBs of the face.
      formatter.read_inner_ccbs_begin();

      DInner_ccb*  new_iccb;
      const Size   n_iccbs = formatter.read_size("number_of_inner_ccbs");
      for (k = 0; k < n_iccbs; k++) {
        // Allocate a new inner CCB record and set its incident face.
        new_iccb = m_arr_access.new_inner_ccb();
        new_iccb->set_face(new_f);

        // Read the current inner CCB.
        n = formatter.read_size("halfedges_on_inner_ccb");
        he = _read_ccb(formatter, n, NULL, new_iccb);
        new_f->add_inner_ccb(new_iccb, he);
      }
      formatter.read_inner_ccbs_end();

      // Read the isolated vertices inside the face.
      Size n_isolated_vertices = 
        formatter.read_size("number_of_isolated_vertices");
      if (n_isolated_vertices) {
        formatter.read_isolated_vertices_begin();
        Size k;
        for (k = 0; k < n_isolated_vertices; k++) {
          // Allocate a new isolated vertex record and set its incident face.
          DIso_vert* new_iso_vert = m_arr_access.new_isolated_vertex();
          new_iso_vert->set_face(new_f);

          // Read the current isolated vertex.
          std::size_t v_idx = formatter.read_vertex_index();
          DVertex* iso_v = m_vertices[v_idx];
          iso_v->set_isolated_vertex(new_iso_vert);
          new_f->add_isolated_vertex(new_iso_vert, iso_v);
        }
        formatter.read_isolated_vertices_end();
      }
    
      // Read any auxiliary data associated with the face.
      if (is_valid)
        formatter.read_face_data(Face_handle(new_f));

      formatter.read_face_end();
    }

    /*!
     * Read a circular boundary of a conncted component.
     * \param formatter The formatter.
     * \param boundary_size The number of halfedges along the boundary.
     * \param p_outer The outer CCB.
     * \param p_inner The inner CCB.
     * \pre p_outer is valid and p_inner is NULL, or vice versa.
     * \return A pointer to the first halfedge read.
     */
    template <class Formatter>
    DHalfedge* _read_ccb(Formatter& formatter, 
                         Size boundary_size,
                         DOuter_ccb* p_outer,
                         DInner_ccb* p_inner)
    {
      CGAL_assertion((p_outer != NULL && p_inner == NULL) ||
                     (p_outer == NULL && p_inner != NULL));

      formatter.read_ccb_halfedges_begin();
 
      // Find the first halfedge, and set its CCB.
      std::size_t  first_idx = formatter.read_halfedge_index();
      DHalfedge*   first_he = m_halfedges [first_idx];

      if (p_outer != NULL)
        first_he->set_outer_ccb(p_outer);
      else
        first_he->set_inner_ccb(p_inner);

      // Read the rest of the halfedge along the boundary.
      std::size_t  curr_idx;
      DHalfedge*   prev_he = first_he;    
      DHalfedge*   curr_he;
      Size         k;

      for (k = 1; k < boundary_size; k++)
      {
        curr_idx = formatter.read_halfedge_index();
        curr_he = m_halfedges[curr_idx];

        // Connect the previous halfedge and the current one.
        prev_he->set_next(curr_he);

        // Set the CCB.
        if (p_outer != NULL)
          curr_he->set_outer_ccb(p_outer);
        else
          curr_he->set_inner_ccb(p_inner);

        prev_he = curr_he;
      }

      // Close the circular list be connecting the first and the last halfedges.
      prev_he->set_next(first_he);

      formatter.read_ccb_halfedges_end();

      // Return the first halfedge.
      return (first_he);
    }
  };

} //namespace CGAL

#endif // CGAL_IO_ARRANGEMENT_2_READER_H 
