// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>
//                 (based on old version by Michal Meyerovitch and Ester Ezra)
//
#ifndef CGAL_IO_ARRANGEMENT_2_WRITER_H
#define CGAL_IO_ARRANGEMENT_2_WRITER_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * The header file for the Arrangement_2_writer<Arrangement> class.
 */

#include <CGAL/Arr_accessor.h>
#include <map>

namespace CGAL {

  /*! \class
   * An auxiliary class for writing an arrangement to an output stream.
   */
  template <class Arrangement_>
  class Arrangement_2_writer
  {
  public:

    typedef Arrangement_                                  Arrangement_2;
    typedef Arrangement_2_writer<Arrangement_2>           Self;

  protected:

    typedef typename Arrangement_2::Size                  Size;
    typedef typename Arrangement_2::Dcel                  Dcel;

    typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
    typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

    typedef CGAL::Arr_accessor<Arrangement_2>             Arr_accessor;
    typedef typename Arr_accessor::Dcel_vertex_iterator   Vertex_const_iterator;
    typedef typename Arr_accessor::Dcel_edge_iterator     Edge_const_iterator;
    typedef typename Arr_accessor::Dcel_face_iterator     Face_const_iterator;

    typedef typename Arr_accessor::Dcel_outer_ccb_iterator
      Outer_ccb_iterator;
    typedef typename Arr_accessor::Dcel_inner_ccb_iterator
      Inner_ccb_iterator;
    typedef typename Arr_accessor::Dcel_iso_vertex_iterator
      Isolated_vertex_iterator;

    typedef typename Arr_accessor::Dcel_vertex            DVertex;
    typedef typename Arr_accessor::Dcel_halfedge          DHalfedge;
    typedef typename Arr_accessor::Dcel_face              DFace;
    typedef std::map<const DVertex*, int>                 Vertex_index_map;
    typedef std::map<const DHalfedge*, int>               Halfedge_index_map;

    // Data memebrs:
    const Arrangement_2&   m_arr;
    const Dcel*            m_dcel;
    int                    m_curr_v;
    Vertex_index_map       m_v_index;
    int                    m_curr_he;
    Halfedge_index_map     m_he_index;

  private:

    // Copy constructor and assignment operator - not supported.
    Arrangement_2_writer(const Self&);
    Self& operator= (const Self&);

  public:

    /*! Constructor. */
    Arrangement_2_writer(const Arrangement_2& arr) :
      m_arr(arr),
      m_dcel(NULL),
      m_curr_v(0),
      m_curr_he(0)
    {
      const Arr_accessor     arr_access(const_cast<Arrangement_2&>(arr));
      m_dcel = &(arr_access.dcel());
    }

    /*! Destructor. */
    virtual ~Arrangement_2_writer()
    {}

    /*! Write the arrangement. */
    template <class Formatter>
    void operator()(Formatter& formatter)
    {
      formatter.write_arrangement_begin();
      formatter.write_size("number_of_vertices",
                           m_dcel->size_of_vertices());
      formatter.write_size("number_of_edges",
                           m_dcel->size_of_halfedges() / 2);
      formatter.write_size("number_of_faces",
                           m_dcel->size_of_faces());

      // Reset indices.
      m_curr_v = 0;
      m_curr_he = 0;

      // Write the vertices.
      formatter.write_vertices_begin();
      Vertex_const_iterator  vit;
      for (vit = m_dcel->vertices_begin(); vit != m_dcel->vertices_end(); ++vit)
      {
        _write_vertex(formatter, vit);
      }
      formatter.write_vertices_end();

      // Write the edges.
      formatter.write_edges_begin();
      Edge_const_iterator    eit;
      for (eit = m_dcel->edges_begin(); eit != m_dcel->edges_end(); ++eit)
      {
        _write_edge(formatter, eit);
      }
      formatter.write_edges_end();

      // Write the faces (the fictitious face first).
      formatter.write_faces_begin();

      Face_const_iterator    fit;
      for (fit = m_dcel->faces_begin(); fit != m_dcel->faces_end(); ++fit)
        _write_face(formatter, fit);
      formatter.write_faces_end();

      formatter.write_arrangement_end();
    }

  protected:

    /*! Write a vertex. */
    template <class Formatter>
    void _write_vertex(Formatter& formatter, Vertex_const_iterator vit)
    {
      // Map the current vertex to its index.
      const DVertex* v = &(*vit);

      m_v_index[v] = m_curr_v;
      ++m_curr_v;

      // Write the vertex.
      formatter.write_vertex_begin();
      formatter.write_vertex_index(static_cast<int>(v->parameter_space_in_x()));
      formatter.write_vertex_index(static_cast<int>(v->parameter_space_in_y()));

      if (! v->has_null_point())
      {
        // Write the associated point.
        formatter.write_vertex_index(1);
        formatter.write_point(v->point());

        // Write additional user-defined data.
        formatter.write_vertex_data(Vertex_const_handle(v));
      }
      else
      {
        // Mark that the vertex is not associated with a point.
        formatter.write_vertex_index(0);
      }

      formatter.write_vertex_end();
    }

    /*! Write an edge (a pair of halfedges). */
    template <class Formatter>
    void _write_edge(Formatter& formatter, Edge_const_iterator hit)
    {
      // Map the halfedge and its twin to their indices.
      const DHalfedge* he = &(*hit);
      const DHalfedge* he_t = he->opposite();

      m_he_index[&(*he)] = m_curr_he;
      ++m_curr_he;
      m_he_index[&(*he_t)] = m_curr_he;
      ++m_curr_he;

      // Write the edge.
      formatter.write_edge_begin();
      formatter.write_vertex_index(_index(he_t->vertex()));
      formatter.write_vertex_index(_index(he->vertex()));

      if (he->direction() == ARR_LEFT_TO_RIGHT)
        formatter.write_vertex_index(0);
      else
        formatter.write_vertex_index(1);

      if (! he->has_null_curve())
      {
        // Write the associated curve.
        formatter.write_vertex_index(1);
        formatter.write_x_monotone_curve(he->curve());

        // Write additional user-defined data.
        formatter.write_halfedge_data(Halfedge_const_handle(he));
        formatter.write_halfedge_data(Halfedge_const_handle(he_t));
      }
      else
      {
        // Mark that the edge is fictitious.
        formatter.write_vertex_index(0);
      }
      formatter.write_edge_end();
    }

    /*! Write a face. */
    template <class Formatter>
    void _write_face(Formatter& formatter, Face_const_iterator fit) const
    {
      const DFace* f = &(*fit);

      formatter.write_face_begin();

      // Write whether the face is unbounded and whether it is valid
      // (non-fictitious).
      if (f->is_unbounded())
        formatter.write_vertex_index(1);
      else
        formatter.write_vertex_index(0);

      if (! f->is_fictitious())
        formatter.write_vertex_index(1);
      else
        formatter.write_vertex_index(0);

      // Write the outer CCBs of the face.
      const std::size_t    n_occbs = f->number_of_outer_ccbs();
      Outer_ccb_iterator   oc_it;

      formatter.write_outer_ccbs_begin();
      formatter.write_size("number_of_outer_ccbs", n_occbs);
      for (oc_it = f->outer_ccbs_begin();
           oc_it != f->outer_ccbs_end(); ++oc_it)
      {
        const std::size_t              n = _circulator_size(*oc_it);

        formatter.write_size("halfedges_on_outer_ccb", n);
        _write_ccb(formatter, *oc_it);
      }
      formatter.write_inner_ccbs_end();

      // Write the inner CCBs of the face.
      const std::size_t    n_iccbs = f->number_of_inner_ccbs();
      Inner_ccb_iterator   ic_it;

      formatter.write_inner_ccbs_begin();
      formatter.write_size("number_of_inner_ccbs", n_iccbs);
      for (ic_it = f->inner_ccbs_begin(); ic_it != f->inner_ccbs_end(); ++ic_it)
      {
        const std::size_t n = _circulator_size(*ic_it);
        formatter.write_size("halfedges_on_inner_ccb", n);
        _write_ccb(formatter, *ic_it);
      }
      formatter.write_inner_ccbs_end();

      // Write the isolated vertices inside the face.
      std::size_t n_isolated = f->number_of_isolated_vertices();
      formatter.write_size("number_of_isolated_vertices", n_isolated);
      if (n_isolated) {
        formatter.write_isolated_vertices_begin();
        Isolated_vertex_iterator iso_vit;
        for (iso_vit = f->isolated_vertices_begin();
             iso_vit != f->isolated_vertices_end(); ++iso_vit)
          formatter.write_vertex_index(_index(&(*iso_vit)));
        formatter.write_isolated_vertices_end();
      }

      // Write additional user-defined data associated with the face.
      if (! f->is_fictitious())
        formatter.write_face_data(Face_const_handle(f));

      formatter.write_face_end();
    }

    /*! Write the edges along a given CCB. */
    template <class Formatter>
    void _write_ccb(Formatter& formatter, const DHalfedge* ccb) const
    {
      const DHalfedge* curr = ccb;

      formatter.write_ccb_halfedges_begin();
      do {
        formatter.write_halfedge_index(_index(curr));
        curr = curr->next();
      } while (curr != ccb);
      formatter.write_ccb_halfedges_end();
    }

    /*! Get the mapped index of a given vertex. */
    int _index(const DVertex* v) const
    {
      typename Vertex_index_map::const_iterator   pos = m_v_index.find(v);

      CGAL_assertion(pos != m_v_index.end());
      return (pos->second);
    }

    /*! Get the mapped index of a given halfegde. */
    int _index(const DHalfedge* he) const
    {
      typename Halfedge_index_map::const_iterator  pos = m_he_index.find(he);

      CGAL_assertion(pos != m_he_index.end());
      return (pos->second);
    }

    /*! Get the number of edges along a given CCB. */
    std::size_t _circulator_size(const DHalfedge* ccb) const
    {
      CGAL_assertion(ccb != NULL);

      std::size_t       n = 0;
      const DHalfedge*  curr = ccb;

      do {
        ++n;
        curr = curr->next();
      } while (curr != ccb);

      return (n);
    }
  };

} //namespace CGAL

#endif // CGAL_IO_ARRANGEMENT_2_WRITER_H
