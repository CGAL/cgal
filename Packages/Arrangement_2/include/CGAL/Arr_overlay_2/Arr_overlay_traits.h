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
// Author(s)     : Ron Wein <baruchzu@post.tau.ac.il>

#ifndef ARR_OVERLAY_TRAITS_H
#define ARR_OVERLAY_TRAITS_H

CGAL_BEGIN_NAMESPACE

/*!
 * \class
 * An overlay-traits class for computing the overlay of two arrangement that
 * are templated with the default DCEL classes, namely they store no extra
 * data with their DCEL features. The resulting arrangement is also assumed
 * to be templated with the default DCEL.
 */
template <class Arrangement1_, class Arrangement2_, class ResArr_>
class _Arr_default_overlay_traits
{
public:

  typedef typename Arrangement1_::Vertex_const_handle    Vertex_handle1;
  typedef typename Arrangement1_::Halfedge_const_handle  Halfedge_handle1;
  typedef typename Arrangement1_::Face_const_handle      Face_handle1;

  typedef typename Arrangement2_::Vertex_const_handle    Vertex_handle2;
  typedef typename Arrangement2_::Halfedge_const_handle  Halfedge_handle2;
  typedef typename Arrangement2_::Face_const_handle      Face_handle2;
  
  typedef typename ResArr_::Vertex_handle                Res_vertex_handle;
  typedef typename ResArr_::Halfedge_handle              Res_halfedge_handle;
  typedef typename ResArr_::Face_handle                  Res_face_handle;

  
  /*!
   * Create a vertex v that corresponds to the coinciding vertices v1 and v2.
   */
  virtual void create_vertex (Vertex_handle1 /* v1 */,
			      Vertex_handle2 /* v2 */,
			      Res_vertex_handle /* v */) const
  {}

  /*!
   * Create a vertex v that mathces v1, which lies of the edge e2.
   */
  virtual void create_vertex (Vertex_handle1 /* v1 */,
			      Halfedge_handle2 /* e2 */,
			      Res_vertex_handle /* v */) const
  {}

  /*!
   * Create a vertex v that mathces v1, contained in the face f2.
   */
  virtual void create_vertex (Vertex_handle1 /* v1 */,
			      Face_handle2 /* f2 */,
			      Res_vertex_handle /* v */) const
  {}

  /*!
   * Create a vertex v that mathces v2, which lies of the edge e1.
   */
  virtual void create_vertex (Halfedge_handle1 /* e1 */,
			      Vertex_handle2 /* v2 */,
			      Res_vertex_handle /* v */) const
  {}

  /*!
   * Create a vertex v that mathces v2, contained in the face f1.
   */
  virtual void create_vertex (Face_handle1 /* f1 */,
			      Vertex_handle2 /* v2 */,
			      Res_vertex_handle /* v */) const
  {}

  /*!
   * Create a vertex v that mathces the intersection of the edges e1 and e2.
   */
  virtual void create_vertex (Halfedge_handle1 /* e1 */,
			      Halfedge_handle2 /* e2 */,
			      Res_vertex_handle /* v */) const
  {}

  /*!
   * Create an edge e that matches the overlap between e1 and e2.
   */
  virtual void create_edge (Halfedge_handle1 /* e1 */,
			    Halfedge_handle2 /* e2 */,
			    Res_halfedge_handle /* e */) const
  {}

  /*!
   * Create an edge e that matches the edge e1, contained in the face f2.
   */
  virtual void create_edge (Halfedge_handle1 /* e1 */,
			    Face_handle2 /* f2 */,
			    Res_halfedge_handle /* e */) const
  {}

  /*!
   * Create an edge e that matches the edge e2, contained in the face f1.
   */
  virtual void create_edge (Face_handle1 /* f1 */,
			    Halfedge_handle2 /* e2 */,
			    Res_halfedge_handle /* e */) const
  {}

  /*!
   * Create a face f that matches the overlapping region between f1 and f2.
   */
  virtual void create_face (Face_handle1 /* f1 */,
			    Face_handle2 /* f2 */,
			    Res_face_handle /* f */) const
  {}

};

CGAL_END_NAMESPACE

#endif
