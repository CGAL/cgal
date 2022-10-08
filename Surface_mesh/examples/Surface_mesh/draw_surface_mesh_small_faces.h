// Copyright (c) 2018 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_SURFACE_MESH_SMALL_FACES_H
#define CGAL_DRAW_SURFACE_MESH_SMALL_FACES_H

#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Surface_mesh.h>
#include <CGAL/Random.h>

#include <cassert>

namespace draw_function_for_surface_mesh
{
template <class SM>
typename CGAL::Kernel_traits<typename SM::Point>::Kernel::Vector_3 get_face_normal(typename SM::Halfedge_index he,
                                          const SM &sm)
{
  typedef typename SM::Halfedge_index halfedge_descriptor;
  typedef typename SM::Point Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;

  typename Kernel::Vector_3 normal=CGAL::NULL_VECTOR;

  halfedge_descriptor end=he;
  unsigned int nb=0;
  do
  {
    CGAL::internal::newell_single_step_3(sm.point(sm.source(he)),
                                          sm.point(sm.target(he)), normal);
    ++nb;
    he=sm.next(he);
  }
  while (he!=end);
  assert(nb>0);
  return (typename Kernel::Construct_scaled_vector_3()(normal, 1.0/nb));
}

template <class SM>
typename CGAL::Kernel_traits<typename SM::Point>::Kernel::Vector_3 get_vertex_normal(typename SM::Halfedge_index he,
                    const SM &sm)
{
  typedef typename SM::Point Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
  typedef typename SM::Halfedge_index halfedge_descriptor;

  typename Kernel::Vector_3 normal=CGAL::NULL_VECTOR;

  halfedge_descriptor end=he;
  do
  {
    if (!sm.is_border(he))
    {
      typename Kernel::Vector_3 n=get_face_normal(he, sm);
      normal=typename Kernel::Construct_sum_of_vectors_3()(normal, n);
    }
    he=sm.next(sm.opposite(he));
  }
  while (he!=end);

  if (!typename Kernel::Equal_3()(normal, CGAL::NULL_VECTOR))
  { normal=(typename Kernel::Construct_scaled_vector_3()
            (normal, 1.0/CGAL::sqrt(normal.squared_length()))); }

  return normal;
}

template <class SM, typename BufferType = float>
void compute_face(typename SM::Face_index fh,
                  const SM &sm,
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                  typename CGAL::Kernel_traits<typename SM::Point>::Kernel::FT & m_min_size,
                  typename CGAL::Kernel_traits<typename SM::Point>::Kernel::FT & m_max_size)
{

  typedef typename SM::Point Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
  typedef typename SM::Face_index face_descriptor;
  typedef typename SM::Halfedge_index halfedge_descriptor;
  typedef typename Kernel::FT FT;

  // TODO: add custom functor, right?
  bool m_draw_small_faces = true;
  bool m_draw_big_faces =  true;

  unsigned int m_threshold = 85;

  // [Face creation]
  bool issmall=false;

  // Default color of faces
  CGAL::IO::Color c(75,160,255);

  // Compare the size of the face with the % m_threshold
  bool exist;
  typename SM::template Property_map<face_descriptor, FT> faces_size;
  boost::tie(faces_size, exist)=sm.template property_map<face_descriptor, FT>("f:size");
  assert(exist);

  // It it is smaller, color the face in red.
  if (get(faces_size, fh)<m_min_size+((m_max_size-m_min_size)/(100-m_threshold)))
  {
    c=CGAL::IO::Color(255,20,20);
    issmall=true;
  }

  if ((issmall && !m_draw_small_faces) || (!issmall && !m_draw_big_faces))
  { return; }

  // Add the color of the face, then all its points.
  graphic_buffer.face_begin(c);
  halfedge_descriptor hd=sm.halfedge(fh);
  do
  {
    graphic_buffer.add_point_in_face(sm.point(sm.source(hd)), get_vertex_normal(hd, sm));
    hd=sm.next(hd);
  }
  while(hd!=sm.halfedge(fh));
  graphic_buffer.face_end();
  /// [Face creation]
}

// Copy from draw_surface_mesh.h
template <class SM, typename BufferType = float>
void compute_edge(typename SM::Edge_index e,
                  const SM &sm,
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  /// [Edge creation]
  graphic_buffer.add_segment(sm.point(sm.source(sm.halfedge(e))),
              sm.point(sm.target(sm.halfedge(e))));
  /// [Edge creation]
}

template <class SM, typename BufferType = float>
void compute_vertex(typename SM::Vertex_index vh,
                    const SM &sm,
                    CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  /// [Vertex creation]
  graphic_buffer.add_point(sm.point(vh));
  /// [Vertex creation]
}

template <class SM, class DrawingFunctor, typename BufferType = float>
void compute_elements(const SM &sm,
                      CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                      const DrawingFunctor &m_drawing_functor)
{

  typedef typename SM::Point Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
  typedef typename SM::Vertex_index vertex_descriptor;
  typedef typename SM::Face_index face_descriptor;
  typedef typename SM::Edge_index edge_descriptor;
  typedef typename SM::Halfedge_index halfedge_descriptor;
  typedef typename Kernel::FT FT;

  FT m_min_size, m_max_size;

  if (sm.faces().begin()!=sm.faces().end())
  {
    bool exist;
    typename SM::template Property_map<face_descriptor, FT> faces_size;
    boost::tie(faces_size, exist)=sm.template property_map<face_descriptor, FT>("f:size");
    assert(exist);

    m_min_size=faces_size[*(sm.faces().begin())];
    m_max_size=m_min_size;
    FT cur_size;

    for (typename SM::Face_range::iterator f=sm.faces().begin(); f!=sm.faces().end(); ++f)
    {
      cur_size=faces_size[*f];
      if (cur_size<m_min_size) m_min_size=cur_size;
      if (cur_size>m_max_size) m_max_size=cur_size;
    }
  }

  for (typename SM::Face_range::iterator f=sm.faces().begin();
        f!=sm.faces().end(); ++f)
  {
    if (*f!=boost::graph_traits<SM>::null_face())
    { compute_face(*f, sm, graphic_buffer, m_min_size, m_max_size); }
  }

  for (typename SM::Edge_range::iterator e=sm.edges().begin();
        e!=sm.edges().end(); ++e)
  { compute_edge(*e, sm, graphic_buffer); }

  for (typename SM::Vertex_range::iterator v=sm.vertices().begin();
        v!=sm.vertices().end(); ++v)
  { compute_vertex(*v, sm, graphic_buffer); }
}

} // draw_function_for_surface_mesh

template <class SM, class DrawingFunctor, typename BufferType = float>
void add_in_graphic_buffer(const SM &sm, CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                               const DrawingFunctor &m_drawing_functor) {
  draw_function_for_surface_mesh::compute_elements(sm, graphic_buffer, m_drawing_functor);
}

template <class SM, typename BufferType = float>
void add_in_graphic_buffer(const SM &sm, CGAL::Graphic_buffer<BufferType> &graphic_buffer) {

  CGAL::Drawing_functor<SM,typename SM::Vertex_index,
                          typename SM::Edge_index,
                          typename SM::Face_index> drawing_functor;

  add_in_graphic_buffer(sm, graphic_buffer, drawing_functor);
}

template<class K, class DrawingFunctor, class BufferType = float>
void draw_surface_mesh_with_small_faces(CGAL::Surface_mesh<K>& amesh,
                                        const DrawingFunctor &drawing_functor)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(amesh, buffer, drawing_functor);
  draw_buffer(buffer);
}

template<class K, class BufferType = float>
void draw_surface_mesh_with_small_faces(CGAL::Surface_mesh<K>& amesh)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(amesh, buffer);
  draw_buffer(buffer);
}

#else // CGAL_USE_BASIC_VIEWER

template<class K, class DrawingFunctor, class BufferType = float>
void draw_surface_mesh_with_small_faces(CGAL::Surface_mesh<K>&)
{
  std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
}

template<class K>
void draw_surface_mesh_with_small_faces(CGAL::Surface_mesh<K>&, const DrawingFunctor &drawing_functor)
{
  std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
}

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_SURFACE_MESH_SMALL_FACES_H
