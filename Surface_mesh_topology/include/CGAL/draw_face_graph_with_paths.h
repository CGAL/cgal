// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
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
//
#ifndef CGAL_DRAW_FACE_GRAPH_WITH_PATHS_H
#define CGAL_DRAW_FACE_GRAPH_WITH_PATHS_H

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#include <iostream>
#include <initializer_list>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/assertions.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Random.h>

namespace CGAL {

  
namespace draw_function_for_lcc
{
  // We need to re-use the namespace draw_function_for_lcc because we want to specialize
  // the previous struct LCC_geom_utils
//   template <class LCC, class Local_kernel, int dim = LCC::ambient_dimension>
// struct LCC_geom_utils;

// Specialisation for face graph; otherwise use the LCC_geom_utils of LCC.
template<class Mesh, class Local_kernel>
struct LCC_geom_utils<CGAL::Face_graph_wrapper<Mesh>, Local_kernel, 3>
{
  static typename Get_traits<Mesh>::Vector
  get_face_normal(const CGAL::Face_graph_wrapper<Mesh>& mesh,
                  typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_handle dh)
  {
    typename Get_traits<Mesh>::Vector normal(CGAL::NULL_VECTOR);
    const typename Get_traits<Mesh>::Point*
        curr=&Get_traits<Mesh>::get_point(mesh.get_fg(), dh);
    typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_handle adart=dh;
    unsigned int nb=0;

    do
    {
      const typename Get_traits<Mesh>::Point*
          next=&Get_traits<Mesh>::get_point(mesh.get_fg(),
                                            mesh.other_extremity(adart));
      internal::newell_single_step_3_for_lcc(*curr, *next, normal);
      ++nb;
      curr=next;
      adart=mesh.next(adart);
    }
    while(adart!=dh);

    CGAL_assertion(nb>0);
    return typename Get_traits<Mesh>::Kernel::Construct_scaled_vector_3()
      (normal, 1.0/nb);
  }
  static typename Local_kernel::Vector_3
  get_vertex_normal(const CGAL::Face_graph_wrapper<Mesh>& mesh,
                    typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_handle dh)
  {
    typename Get_traits<Mesh>::Vector normal(CGAL::NULL_VECTOR);
    unsigned int nb = 0;

    for ( typename CGAL::Face_graph_wrapper<Mesh>::template Dart_of_cell_range<0>::
          const_iterator it=mesh.template darts_of_cell<0>(dh).begin(),
          itend=mesh.template darts_of_cell<0>(dh).end(); it!=itend; ++it )
    {
      normal=typename Get_traits<Mesh>::Kernel::Construct_sum_of_vectors_3()
          (normal, get_face_normal(mesh, it));
      ++nb;
    }

    if ( nb<2 ) return internal::Geom_utils
                  <typename Get_traits<Mesh>::Kernel, Local_kernel>::
                  get_local_vector(normal);

    return internal::Geom_utils
      <typename Get_traits<Mesh>::Kernel, Local_kernel>::
      get_local_vector(typename Get_traits<Mesh>::Kernel::
                       Construct_scaled_vector_3()(normal, 1.0/nb));
  }
};

} // namespace draw_function_for_lcc
  
namespace draw_function_for_face_graph_with_paths
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3 Local_point;
typedef Local_kernel::Vector_3 Local_vector;
  
// Destructor.
// lcc.free_mark(m_oriented_mark);



template <typename Mesh>
const typename CGAL::Get_traits<Mesh>::Point& get_point(typename Get_map<Mesh, Mesh>::type::Dart_const_handle dh, const Mesh &mesh)
{
  typedef typename CGAL::Get_traits<Mesh>::Kernel Kernel;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;
  typedef typename CGAL::Get_traits<Mesh>::Vector Vector;

  return CGAL::Get_traits<Mesh>::get_point(mesh, dh);
}

template <typename Mesh, typename BufferType = float>
void compute_face(typename Get_map<Mesh, Mesh>::type::Dart_const_handle dh,
                  const Mesh &mesh,
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                  const typename Get_map<Mesh, Mesh>::storage_type& lcc)
{
  typedef typename Get_map<Mesh, Mesh>::type      LCC;

  typedef typename LCC::size_type                 size_type;

  typedef typename LCC::Dart_const_handle         Dart_const_handle;
  typedef typename CGAL::Get_traits<Mesh>::Kernel Kernel;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;
  typedef typename CGAL::Get_traits<Mesh>::Vector Vector;

  // We fill only closed faces.
  Dart_const_handle cur=dh;
  Dart_const_handle min=dh;
  do
  {
    if (!lcc.is_next_exist(cur)) return; // open face=>not filled
    if (cur<min) min=cur;
    cur=lcc.next(cur);
  }
  while(cur!=dh);

  // CGAL::IO::Color c=m_fcolor.run(*lcc, dh);
  graphic_buffer.face_begin(); //c);

  cur=dh;
  do
  {
    graphic_buffer.add_point_in_face(draw_function_for_face_graph_with_paths::get_point(cur, mesh),
                                     draw_function_for_lcc::LCC_geom_utils<LCC, Local_kernel>::
                      get_vertex_normal(lcc, cur));
    cur=lcc.next(cur);
  }
  while(cur!=dh);

  graphic_buffer.face_end();
}


  // typename LCC::size_type m_amark; // If !=INVALID_MARK, show darts marked with this mark

template <typename Mesh, typename BufferType = float>
void compute_edge(typename Get_map<Mesh, Mesh>::type::Dart_const_handle dh,
                  const Mesh &mesh,
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                  const typename Get_map<Mesh, Mesh>::storage_type& lcc,
                  typename Get_map<Mesh, Mesh>::type::size_type m_amark,
                  bool m_draw_marked_darts = true)
{
  typedef typename Get_map<Mesh, Mesh>::type      LCC;
  typedef typename LCC::size_type                 size_type;

  typedef typename LCC::Dart_const_handle         Dart_const_handle;
  typedef typename CGAL::Get_traits<Mesh>::Kernel Kernel;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;
  typedef typename CGAL::Get_traits<Mesh>::Vector Vector;

  Point p1 = get_point(dh, mesh);
  Dart_const_handle d2 = lcc.other_extremity(dh);
  if (d2!=LCC::null_handle)
  {
    if (m_draw_marked_darts && m_amark!=LCC::INVALID_MARK &&
        (lcc.is_marked(dh, m_amark) || lcc.is_marked(lcc.opposite2(dh), m_amark)))
    { graphic_buffer.add_segment(p1, get_point(d2, mesh), CGAL::IO::Color(0, 0, 255)); }
    else
    { graphic_buffer.add_segment(p1, get_point(d2, mesh)); }
  }
}

template <typename Mesh, typename BufferType = float>
void compute_edge(typename Get_map<Mesh, Mesh>::type::Dart_const_handle dh,
                  const CGAL::IO::Color& color,
                  const Mesh &mesh,
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                  const typename Get_map<Mesh, Mesh>::storage_type& lcc)
{
  typedef typename Get_map<Mesh, Mesh>::type      LCC;
  typedef typename LCC::Dart_const_handle         Dart_const_handle;
  typedef typename CGAL::Get_traits<Mesh>::Kernel Kernel;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;
  typedef typename CGAL::Get_traits<Mesh>::Vector Vector;

  Point p1 = get_point(dh, mesh);
  Dart_const_handle d2 = lcc.other_extremity(dh);
  if (d2!=LCC::null_handle)
  { graphic_buffer.add_segment(p1, get_point(d2, mesh), color); }
}

template <typename Mesh, typename BufferType = float>
void compute_vertex(typename Get_map<Mesh, Mesh>::type::Dart_const_handle dh,
                    const Mesh &mesh,
                    CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  typedef typename CGAL::Get_traits<Mesh>::Kernel Kernel;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;
  typedef typename CGAL::Get_traits<Mesh>::Vector Vector;
  graphic_buffer.add_point(get_point(dh, mesh));
}

template <typename Mesh, typename BufferType = float>
void compute_path(std::size_t i,
                  typename Get_map<Mesh, Mesh>::type::size_type amark,
                  const Mesh &mesh,
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                  const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>* m_paths,
                  const typename Get_map<Mesh, Mesh>::storage_type& lcc)
{

  typedef typename CGAL::Get_traits<Mesh>::Kernel Kernel;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;
  typedef typename CGAL::Get_traits<Mesh>::Vector Vector;

  if ((*m_paths)[i].is_empty())
  { return; }

  CGAL::Random random(static_cast<unsigned int>(i));
  CGAL::IO::Color color = get_random_color(random);

  graphic_buffer.add_point(get_point((*m_paths)[i].get_ith_dart(0), mesh), color);
  for (std::size_t j=0; j<(*m_paths)[i].length(); ++j)
  {
    if ( !lcc.is_marked( (*m_paths)[i].get_ith_dart(j), amark) )
    {
      compute_edge<Mesh, BufferType>((*m_paths)[i].get_ith_dart(j), color, mesh, graphic_buffer, lcc);
      lcc.template mark_cell<1>((*m_paths)[i].get_ith_dart(j), amark);
    }
  }
}

template <class Mesh, class DrawingFunctor, typename BufferType = float>
void compute_elements(const Mesh &mesh,
                      CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                      const DrawingFunctor &m_drawing_functor,
                      const typename Get_map<Mesh, Mesh>::storage_type& lcc,
                      const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>* m_paths,
                      typename Get_map<Mesh, Mesh>::type::size_type amark /*= typename Get_map<Mesh, Mesh>::type::INVALID_MARK*/,
                      bool m_nofaces = false,
                      bool m_draw_marked_darts = true)
{
  typedef typename Get_map<Mesh, Mesh>::type      LCC;
  typedef typename LCC::size_type                 size_type;

  typedef typename CGAL::Get_traits<Mesh>::Kernel Kernel;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;
  typedef typename CGAL::Get_traits<Mesh>::Vector Vector;

  typedef typename LCC::Dart_const_handle         Dart_const_handle;

  typename LCC::Dart_range::const_iterator m_current_dart = lcc.darts().end();
  typename LCC::size_type m_oriented_mark = lcc.get_new_mark();
  std::size_t m_current_path = m_paths->size();
  typename LCC::size_type m_amark = amark==(std::numeric_limits<size_type>::max)()?
            LCC::INVALID_MARK:amark; // If !=INVALID_MARK, show darts marked with this mark

  lcc.orient(m_oriented_mark);


  typename LCC::size_type markfaces    = lcc.get_new_mark();
  typename LCC::size_type markedges    = lcc.get_new_mark();
  typename LCC::size_type markvertices = lcc.get_new_mark();

  if (m_current_dart!=lcc.darts().end())
  { // We want to draw only one dart
    Dart_const_handle selected_dart=m_current_dart; //lcc.dart_handle(m_current_dart);
    compute_edge<Mesh, BufferType>(selected_dart, CGAL::IO::Color(255,0,0), mesh, graphic_buffer, lcc);
    lcc.template mark_cell<1>(selected_dart, markedges);
    compute_vertex<Mesh, BufferType>(selected_dart, mesh, graphic_buffer);

    if ( !m_nofaces )
    { compute_face<Mesh, BufferType>(selected_dart, mesh, graphic_buffer, lcc); }

    for (typename LCC::Dart_range::const_iterator it=lcc.darts().begin(),
          itend=lcc.darts().end(); it!=itend; ++it )
    {
      if ( !lcc.is_marked(it, markedges) )
      {
        compute_edge<Mesh, BufferType>(it, mesh, graphic_buffer, lcc, m_amark, m_draw_marked_darts);
        lcc.template mark_cell<1>(it, markedges);
      }
    }
  }
  else
  {
    if (m_current_path==m_paths->size())
    {
      for (std::size_t i=0; i<m_paths->size(); ++i)
      { compute_path<Mesh, BufferType>(i, markedges, mesh, graphic_buffer, m_paths, lcc); }
    }
    else if (m_current_path!=m_paths->size()+1)
    { compute_path<Mesh, BufferType>(m_current_path, markedges, mesh, graphic_buffer, m_paths, lcc); }

    for (typename LCC::Dart_range::const_iterator it=lcc.darts().begin(),
          itend=lcc.darts().end(); it!=itend; ++it )
    {
      if (!m_nofaces && !lcc.is_marked(it, markfaces) &&
          !lcc.is_perforated(it) && lcc.is_marked(it, m_oriented_mark))
      {
        compute_face<Mesh, BufferType>(it, mesh, graphic_buffer, lcc);
        lcc.template mark_cell<2>(it, markfaces);
      }

      if ( !lcc.is_marked(it, markedges) )
      {
        compute_edge<Mesh, BufferType>(it, mesh, graphic_buffer, lcc, m_amark, m_draw_marked_darts);
        lcc.template mark_cell<1>(it, markedges);
      }

      if ( !lcc.is_marked(it, markvertices) )
      {
        compute_vertex<Mesh, BufferType>(it, mesh, graphic_buffer);
        lcc.template mark_cell<0>(it, markvertices);
      }
    }
  }

  lcc.free_mark(markfaces);
  lcc.free_mark(markedges);
  lcc.free_mark(markvertices);
}

} // namespace draw_function_for_face_graph_with_paths

template <typename BufferType = float, class Mesh, class DrawingFunctor>
void add_in_graphic_buffer(const Mesh &mesh,
                      CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                      const DrawingFunctor &m_drawing_functor,
                      // TODO: I think I need to use smart pointers with lcc, right?
                      const typename Get_map<Mesh, Mesh>::storage_type& lcc,
                      const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>* m_paths,
                      typename Get_map<Mesh, Mesh>::type::size_type amark= /*typename*/ Get_map<Mesh, Mesh>::type::INVALID_MARK,
                      bool m_nofaces = false) {
  draw_function_for_face_graph_with_paths::compute_elements(mesh, graphic_buffer, m_drawing_functor, lcc, m_paths, amark, m_nofaces);
}

template <typename BufferType = float, class Mesh>
void add_in_graphic_buffer(const Mesh &mesh,
                      CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                      // TODO: I think I need to use smart pointers with lcc, right?
                      const typename Get_map<Mesh, Mesh>::storage_type& lcc,
                      const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>* m_paths,
                      typename Get_map<Mesh, Mesh>::type::size_type amark= /*typename*/ Get_map<Mesh, Mesh>::type::INVALID_MARK,
                      bool m_nofaces = false) {

  // Default functor; user can add his own functor.
  Drawing_functor<Mesh,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_handle /*vh*/,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_handle /*eh*/,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_handle /*fh*/>
      drawing_functor;

  add_in_graphic_buffer(mesh, graphic_buffer, drawing_functor, lcc, m_paths, amark, m_nofaces);
}


template<typename Mesh, typename BufferType = float >
void draw(const Mesh& alcc,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& paths,
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          const char* title="Mesh Viewer With Path",
          bool nofill=false) {

  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(alcc, buffer, alcc, &paths, amark, nofill);
  draw_buffer(buffer);
}


template<typename Mesh, typename DrawingFunctor, typename BufferType = float>
void draw(const Mesh& alcc,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& paths,
          const DrawingFunctor& drawing_functor,
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          const char* title="Mesh Viewer With Path",
          bool nofill=false) {

  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(alcc, buffer, drawing_functor, alcc, &paths, amark, nofill);
  draw_buffer(buffer);
}


template<class Mesh, typename BufferType = float >
void draw(const Mesh& mesh,
          std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>> l,
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          const char* title="Mesh Viewer With Path",
          bool nofill=false)
{
  std::vector<Surface_mesh_topology::Path_on_surface<Mesh>> paths=l;

  CGAL::Graphic_buffer<BufferType> buffer;
  typename Get_map<Mesh, Mesh>::storage_type alcc(mesh);
  add_in_graphic_buffer(mesh, buffer, alcc, &paths, amark, nofill);
  draw_buffer(buffer);
}

/*
template<class Mesh, class DrawingFunctor>
void draw(const Mesh& alcc,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& paths,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          bool nofill=false,
          const DrawingFunctor& drawing_functor=DrawingFunctor())
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    CGAL::Qt::init_ogl_context(4,3);
    int argc=1;
    const char* argv[1]={"lccviewer"};
    QApplication app(argc,const_cast<char**>(argv));
    Face_graph_with_path_viewer<Mesh, DrawingFunctor> mainwindow(app.activeWindow(),
                                                                 alcc, &paths, amark,
                                                                 title, nofill,
                                                                 drawing_functor);
    mainwindow.show();
    app.exec();
  }
}

template<class Mesh>
void draw(const Mesh& alcc,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& paths,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          bool nofill=false)
{
  DefaultDrawingFunctorLCC f;
  draw<Mesh, DefaultDrawingFunctorLCC>(alcc, paths, title, amark, nofill, f);
}

template<class Mesh>
void draw(const Mesh& alcc,
          std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>> l,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          bool nofill=false)
{
  std::vector<Surface_mesh_topology::Path_on_surface<Mesh> > paths=l;
  draw(alcc, paths, title, amark, nofill);
}
*/

} // End namespace CGAL

#else  // CGAL_USE_BASIC_VIEWER

namespace CGAL
{

  template<class Mesh>
  void draw(const Mesh&,
            const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& ,
            const char* ="",
            typename Get_map<Mesh, Mesh>::type::size_type=
            (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
            bool=false)
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }

  template<class Mesh>
  void draw(const Mesh&,
            std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>>,
            const char* ="",
            typename Get_map<Mesh, Mesh>::type::size_type=
            (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
            bool=false)
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_FACE_GRAPH_WITH_PATHS_H
