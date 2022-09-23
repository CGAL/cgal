// Copyright (c) 2022 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//            Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_LCC_H
#define CGAL_DRAW_LCC_H

#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <functional>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Random.h>

namespace CGAL {

  
  //   CGAL::IO::Color face_color(const DS &aVal, face_handle dh) const {
  //   CGAL::Random random((unsigned int)(aVal.darts().index(dh)));
  //   return get_random_color(random);
  // }

namespace draw_function_for_lcc
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3 Local_point;
typedef Local_kernel::Vector_3 Local_vector;

template <class LCC, class Local_kernel, int dim = LCC::ambient_dimension>
struct LCC_geom_utils;

template <class LCC, class Local_kernel>
struct LCC_geom_utils<LCC, Local_kernel, 3>
{
  static typename Local_kernel::Vector_3
  get_vertex_normal(const LCC &lcc, typename LCC::Dart_const_handle dh)
  {
    typename Local_kernel::Vector_3 n =
        internal::Geom_utils<typename LCC::Traits, Local_kernel>::
            get_local_vector(CGAL::compute_normal_of_cell_0<LCC>(lcc, dh));
    n = n / (CGAL::sqrt(n * n));
    return n;
  }
};

template <class LCC, class Local_kernel>
struct LCC_geom_utils<LCC, Local_kernel, 2>
{
  static typename Local_kernel::Vector_3
  get_vertex_normal(const LCC &, typename LCC::Dart_const_handle)
  {
    typename Local_kernel::Vector_3 n=CGAL::NULL_VECTOR;
    return n;
  }
};

template <typename BufferType=float, class LCC, class DrawingFunctorLCC>
void compute_face(typename LCC::Dart_const_handle dh,
                  typename LCC::Dart_const_handle voldh, const LCC *lcc,
                  const DrawingFunctorLCC &m_drawing_functor,
                  GraphicBuffer<BufferType> &graphic_buffer)
{
  if (!m_drawing_functor.are_faces_enabled() ||
      !m_drawing_functor.draw_face(*lcc, dh))
  { return; }

  // We fill only closed faces.
  typename LCC::Dart_const_handle cur=dh;
  do
  {
    if (!lcc->is_next_exist(cur))
    { return; } // open face=>not filled
    cur = lcc->next(cur);
  }
  while (cur!=dh);

  if (m_drawing_functor.colored_volume(*lcc, voldh))
  {
    CGAL::IO::Color c=m_drawing_functor.volume_color(*lcc, voldh);
    graphic_buffer.face_begin(c);
  }
  else if (m_drawing_functor.colored_face(*lcc, dh))
  {
    CGAL::IO::Color c=m_drawing_functor.face_color(*lcc, dh);
    graphic_buffer.face_begin(c);
  }
  else 
  { graphic_buffer.face_begin(); }

  cur=dh;
  do
  {
    graphic_buffer.add_point_in_face
      (lcc->point(cur),
       LCC_geom_utils<LCC, Local_kernel>::get_vertex_normal(*lcc, cur));
    cur=lcc->next(cur);
  }
  while (cur!=dh);

  graphic_buffer.face_end();
}

template <typename BufferType=float, class LCC, class DrawingFunctor>
void compute_edge(typename LCC::Dart_const_handle dh, const LCC *lcc,
                  const DrawingFunctor &m_drawing_functor,
                  GraphicBuffer<BufferType> &graphic_buffer)
{
  if (!m_drawing_functor.are_edges_enabled() ||
      !m_drawing_functor.draw_edge(*lcc, dh))
  { return; }

  const typename LCC::Point& p1=lcc->point(dh);
  typename LCC::Dart_const_handle d2=lcc->other_extremity(dh);
  if (d2!=nullptr)
  {
    if (m_drawing_functor.colored_edge(*lcc, dh))
    {
      graphic_buffer.add_segment(p1, lcc->point(d2),
                                 m_drawing_functor.edge_color(*lcc, dh));
    }
    else
    { graphic_buffer.add_segment(p1, lcc->point(d2)); }
  }
}

template <typename BufferType = float, class LCC, class DrawingFunctorLCC>
void compute_vertex(typename LCC::Dart_const_handle dh, const LCC *lcc,
                    const DrawingFunctorLCC &m_drawing_functor,
                    GraphicBuffer<BufferType> &graphic_buffer)
{
  if (!m_drawing_functor.are_vertices_enabled() ||
      !m_drawing_functor.draw_vertex(*lcc, dh))
  { return; }

  if (m_drawing_functor.colored_vertex(*lcc, dh))
  {
    graphic_buffer.add_point(lcc->point(dh),
                             m_drawing_functor.vertex_color(*lcc, dh));
  }
  else
  { graphic_buffer.add_point(lcc->point(dh)); }
}

template <typename BufferType = float, class LCC, class DrawingFunctor>
void compute_elements(GraphicBuffer<BufferType> &graphic_buffer, const LCC *lcc,
                      const DrawingFunctor &m_drawing_functor)
{
  if (lcc==nullptr)
  { return; }

  typename LCC::size_type markvolumes = lcc->get_new_mark();
  typename LCC::size_type markfaces = lcc->get_new_mark();
  typename LCC::size_type markedges = lcc->get_new_mark();
  typename LCC::size_type markvertices = lcc->get_new_mark();
  typename LCC::size_type oriented_mark = lcc->get_new_mark();

  lcc->orient(oriented_mark);

  for(typename LCC::Dart_range::const_iterator it=lcc->darts().begin(),
        itend=lcc->darts().end(); it!=itend; ++it)
  {
    if (!lcc->is_marked(it, markvolumes) &&
        m_drawing_functor.draw_volume(*lcc, it))
    {
      for(typename LCC::template Dart_of_cell_basic_range<3>::const_iterator
            itv=lcc->template darts_of_cell_basic<3>(it, markvolumes).begin(),
            itvend=lcc->template darts_of_cell_basic<3>(it, markvolumes).end();
          itv!=itvend; ++itv)
      {
        lcc->mark(itv, markvolumes);
        if (!lcc->is_marked(itv, markfaces) &&
            lcc->is_marked(itv, oriented_mark) &&
            m_drawing_functor.draw_face(*lcc, itv))
        {
          if (!m_drawing_functor.volume_wireframe(*lcc, itv) &&
              !m_drawing_functor.face_wireframe(*lcc, itv))
          { compute_face(itv, it, lcc, m_drawing_functor, graphic_buffer); }
          for(typename LCC::template Dart_of_cell_basic_range<2>::const_iterator
                itf=lcc->template darts_of_cell_basic<2>(itv, markfaces).begin(),
                itfend=lcc->template darts_of_cell_basic<2>(itv, markfaces).end();
              itf!=itfend; ++itf)
          {
            lcc->mark(itf, markfaces);
            if (!lcc->is_marked(itf, markedges) &&
                m_drawing_functor.draw_edge(*lcc, itf))
            {
              compute_edge(itf, lcc, m_drawing_functor, graphic_buffer);
              for(typename LCC::template Dart_of_cell_basic_range<1>::const_iterator
                    ite=lcc->template darts_of_cell_basic<1>(itf, markedges).begin(),
                    iteend=lcc->template darts_of_cell_basic<1>(itf, markedges).end();
                  ite!=iteend; ++ite)
              {
                lcc->mark(ite, markedges);
                if (!lcc->is_marked(ite, markvertices) &&
                    m_drawing_functor.draw_vertex(*lcc, ite))
                {
                  compute_vertex(ite, lcc, m_drawing_functor, graphic_buffer);
                  CGAL::mark_cell<LCC, 0>(*lcc, ite, markvertices);
                }
              }
            }
          }
        }
      }
    }
  }

  for (typename LCC::Dart_range::const_iterator it = lcc->darts().begin(),
         itend = lcc->darts().end(); it != itend; ++it)
  {
    lcc->unmark(it, markvertices);
    lcc->unmark(it, markedges);
    lcc->unmark(it, markfaces);
    lcc->unmark(it, markvolumes);
    lcc->unmark(it, oriented_mark);
  }

  lcc->free_mark(markvolumes);
  lcc->free_mark(markfaces);
  lcc->free_mark(markedges);
  lcc->free_mark(markvertices);
  lcc->free_mark(oriented_mark);
}

} // namespace draw_function

/**
 * @brief This function is responsible for filling the buffer to allow
 * visualization.
 *
 * @param graphic_buffer
 * @param m_drawing_functor
 * @param alcc
 */
template <typename BufferType = float, class LCC, class DrawingFunctor>
void add_in_graphic_buffer_lcc(GraphicBuffer<BufferType> &graphic_buffer,
                               const DrawingFunctor &m_drawing_functor,
                               const LCC *alcc = nullptr)
{
  if (alcc!=nullptr)
  {
    draw_function_for_lcc::compute_elements(graphic_buffer, alcc, m_drawing_functor);
  }
}

// Specialization of draw function.
#define CGAL_LCC_TYPE                                                          \
  CGAL::Linear_cell_complex_base<d_, ambient_dim, Traits_, Items_, Alloc_,     \
                                 Map, Refs, Storage_>

template<unsigned int d_, unsigned int ambient_dim, class Traits_,
         class Items_, class Alloc_,
         template <unsigned int, class, class, class, class> class Map,
         class Refs, class Storage_,
         class DrawingFunctor=Drawing_functor_with_volume<CGAL_LCC_TYPE,
                                                               typename CGAL_LCC_TYPE::Dart_const_handle,
                                                               typename CGAL_LCC_TYPE::Dart_const_handle,
                                                               typename CGAL_LCC_TYPE::Dart_const_handle,
                                                               typename CGAL_LCC_TYPE::Dart_const_handle>>
void draw(const CGAL_LCC_TYPE &alcc,
          const char *title = "LCC for CMap Basic Viewer",
          const DrawingFunctor &drawing_functor=DrawingFunctor())
{
  GraphicBuffer<float> buffer;
  add_in_graphic_buffer_lcc(buffer, drawing_functor, &alcc);
  draw_buffer(buffer);
}

#undef CGAL_LCC_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_LCC_H
