// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <baruchzu@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

/*! \file
 * Definition of the internal _Arr_default_overlay_traits_base class template.
 */

#ifndef CGAL_ARR_DEFAULT_OVERLAY_TRAITS_BASE_H
#define CGAL_ARR_DEFAULT_OVERLAY_TRAITS_BASE_H

namespace CGAL {

/*!
 * \class
 * An overlay-traits class for computing the overlay of two arrangement that
 * are templated with the default DCEL classes, namely they store no extra
 * data with their DCEL features. The resulting arrangement is also assumed
 * to be templated with the default DCEL.
 */
template <class ArrangementA, class ArrangementB, class ArrangementR>
class _Arr_default_overlay_traits_base
{
public:

  typedef typename ArrangementA::Vertex_const_handle    Vertex_handle_A;
  typedef typename ArrangementA::Halfedge_const_handle  Halfedge_handle_A;
  typedef typename ArrangementA::Face_const_handle      Face_handle_A;

  typedef typename ArrangementB::Vertex_const_handle    Vertex_handle_B;
  typedef typename ArrangementB::Halfedge_const_handle  Halfedge_handle_B;
  typedef typename ArrangementB::Face_const_handle      Face_handle_B;
  
  typedef typename ArrangementR::Vertex_handle          Vertex_handle_R;
  typedef typename ArrangementR::Halfedge_handle        Halfedge_handle_R;
  typedef typename ArrangementR::Face_handle            Face_handle_R;

  /*! Destructor. */
  virtual ~_Arr_default_overlay_traits_base ()
  {}
  
  /*!
   * Create a vertex v that corresponds to the coinciding vertices v1 and v2.
   */
  virtual void create_vertex (Vertex_handle_A /* v1 */,
			      Vertex_handle_B /* v2 */,
			      Vertex_handle_R /* v */) const
  {}

  /*!
   * Create a vertex v that mathces v1, which lies of the edge e2.
   */
  virtual void create_vertex (Vertex_handle_A /* v1 */,
			      Halfedge_handle_B /* e2 */,
			      Vertex_handle_R /* v */) const
  {}

  /*!
   * Create a vertex v that mathces v1, contained in the face f2.
   */
  virtual void create_vertex (Vertex_handle_A /* v1 */,
			      Face_handle_B /* f2 */,
			      Vertex_handle_R /* v */) const
  {}

  /*!
   * Create a vertex v that mathces v2, which lies of the edge e1.
   */
  virtual void create_vertex (Halfedge_handle_A /* e1 */,
			      Vertex_handle_B /* v2 */,
			      Vertex_handle_R /* v */) const
  {}

  /*!
   * Create a vertex v that mathces v2, contained in the face f1.
   */
  virtual void create_vertex (Face_handle_A /* f1 */,
			      Vertex_handle_B /* v2 */,
			      Vertex_handle_R /* v */) const
  {}

  /*!
   * Create a vertex v that mathces the intersection of the edges e1 and e2.
   */
  virtual void create_vertex (Halfedge_handle_A /* e1 */,
			      Halfedge_handle_B /* e2 */,
			      Vertex_handle_R /* v */) const
  {}

  /*!
   * Create an edge e that matches the overlap between e1 and e2.
   */
  virtual void create_edge (Halfedge_handle_A /* e1 */,
			    Halfedge_handle_B /* e2 */,
			    Halfedge_handle_R /* e */) const
  {}

  /*!
   * Create an edge e that matches the edge e1, contained in the face f2.
   */
  virtual void create_edge (Halfedge_handle_A /* e1 */,
			    Face_handle_B /* f2 */,
			    Halfedge_handle_R /* e */) const
  {}

  /*!
   * Create an edge e that matches the edge e2, contained in the face f1.
   */
  virtual void create_edge (Face_handle_A /* f1 */,
			    Halfedge_handle_B /* e2 */,
			    Halfedge_handle_R /* e */) const
  {}

  /*!
   * Create a face f that matches the overlapping region between f1 and f2.
   */
  virtual void create_face (Face_handle_A /* f1 */,
			    Face_handle_B /* f2 */,
			    Face_handle_R /* f */) const
  {}

};

} //namespace CGAL

#endif
