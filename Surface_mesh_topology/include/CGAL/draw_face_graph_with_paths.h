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

#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>

#include <iostream>
#include <initializer_list>
#include <functional>

#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/assertions.h>
#include <CGAL/Random.h>

namespace CGAL {

// Specific graphics scene options
template <typename DS,
          typename vertex_handle,
          typename edge_handle,
          typename face_handle>
struct Graphics_scene_options_face_graph_with_paths :
    public CGAL::Graphics_scene_options<DS, vertex_handle, edge_handle, face_handle>
{
  Graphics_scene_options_face_graph_with_paths()
  {
    color_of_path=[](std::size_t i)->CGAL::IO::Color
    {
      CGAL::Random random(static_cast<unsigned int>(i));
      return get_random_color(random);
    };

    draw_path=[](std::size_t)->bool
    { return true; };
  }

  const CGAL::IO::Color& color_of_marked_edges() const
  { return m_color_of_marked_edges; }
  void color_of_marked_edges(const CGAL::IO::Color& c)
  { m_color_of_marked_edges=c; }

  std::function<CGAL::IO::Color(std::size_t)> color_of_path;
  std::function<bool(std::size_t)> draw_path;
protected:
  CGAL::IO::Color m_color_of_marked_edges=CGAL::IO::Color(0, 0, 255);
};

namespace draw_function_for_lcc
{
  // We need to reuse the namespace draw_function_for_lcc because we want to specialize
  // the previous struct LCC_geom_utils
//   template <class LCC, class Local_kernel, int dim = LCC::ambient_dimension>
// struct LCC_geom_utils;

// Specialization for face graph; otherwise use the LCC_geom_utils of LCC.
template<class Mesh, class Local_kernel>
struct LCC_geom_utils<CGAL::Face_graph_wrapper<Mesh>, Local_kernel, 3>
{
  static typename Get_traits<Mesh>::Vector
  get_face_normal(const CGAL::Face_graph_wrapper<Mesh>& mesh,
                  typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_descriptor dh)
  {
    typename Get_traits<Mesh>::Vector normal(CGAL::NULL_VECTOR);
    const typename Get_traits<Mesh>::Point*
        curr=&Get_traits<Mesh>::get_point(mesh.get_fg(), dh);
    typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_descriptor adart=dh;
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
                    typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_descriptor dh)
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

template <typename Mesh>
const typename CGAL::Get_traits<Mesh>::Point& get_point
(const Mesh &mesh, typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor dh)
{ return CGAL::Get_traits<Mesh>::get_point(mesh, dh); }

template <typename Mesh, class GSOptions>
void compute_face(const Mesh& mesh,
                  const typename Get_map<Mesh, Mesh>::storage_type& lcc,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor dh,
                  CGAL::Graphics_scene& graphics_scene,
                  GSOptions& gso)
{
  if(!gso.are_faces_enabled() || !gso.draw_face(mesh, dh))
  { return; }

  typedef typename Get_map<Mesh, Mesh>::type  LCC;
  typedef typename LCC::Dart_const_descriptor Dart_const_descriptor;

  // We fill only closed faces.
  Dart_const_descriptor cur=dh;
  Dart_const_descriptor min=dh;
  do
  {
    if (!lcc.is_next_exist(cur)) return; // open face=>not filled
    if (cur<min) min=cur;
    cur=lcc.next(cur);
  }
  while(cur!=dh);

  if(gso.colored_face(mesh, dh))
  { graphics_scene.face_begin(gso.face_color(mesh, dh)); }
  else
  { graphics_scene.face_begin(); }

  cur=dh;
  do
  {
    graphics_scene.add_point_in_face(draw_function_for_face_graph_with_paths::get_point(mesh, cur),
                                     draw_function_for_lcc::LCC_geom_utils<LCC, Local_kernel>::
                                     get_vertex_normal(lcc, cur));
    cur=lcc.next(cur);
  }
  while(cur!=dh);

  graphics_scene.face_end();
}

template <typename Mesh, class GSOptions>
void compute_edge(const Mesh &mesh,
                  const typename Get_map<Mesh, Mesh>::storage_type& lcc,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor dh,
                  CGAL::Graphics_scene& graphics_scene,
                  GSOptions& gso,
                  bool colored=false,
                  const CGAL::IO::Color& color=CGAL::IO::black())
{
  if(!gso.are_edges_enabled() || !gso.draw_edge(mesh, dh))
  { return; }

  typedef typename Get_map<Mesh, Mesh>::type      LCC;
  typedef typename LCC::Dart_const_descriptor     Dart_const_descriptor;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;

  Point p1=get_point(mesh, dh);
  Dart_const_descriptor d2=lcc.other_extremity(dh);
  if (d2!=LCC::null_descriptor)
  {
    if(colored)
    { graphics_scene.add_segment(p1, get_point(mesh, d2), color); }
    else
    { graphics_scene.add_segment(p1, get_point(mesh, d2)); }
  }
}

template <typename Mesh, class GSOptions>
void compute_edge(const Mesh &mesh,
                  const typename Get_map<Mesh, Mesh>::storage_type& lcc,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor dh,
                  CGAL::Graphics_scene& graphics_scene,
                  GSOptions& gso,
                  typename Get_map<Mesh, Mesh>::type::size_type amark)
{
  if(!gso.are_edges_enabled() || !gso.draw_edge(mesh, dh))
  { return; }

  if (amark!=Get_map<Mesh, Mesh>::type::INVALID_MARK &&
      (lcc.is_marked(dh, amark) || lcc.is_marked(lcc.opposite2(dh), amark)))
  { compute_edge(mesh, lcc, dh, graphics_scene, gso, true, gso.color_of_marked_edges()); }
  else
  {
    if(gso.colored_edge(mesh, dh))
    { compute_edge(mesh, lcc, dh, graphics_scene, gso, true,
                   gso.edge_color(mesh, dh)); }
    else
    { compute_edge(mesh, lcc, dh, graphics_scene, gso, false); }
  }
}

template <typename Mesh, class GSOptions>
void compute_vertex(const Mesh &mesh,
                    typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor dh,
                    CGAL::Graphics_scene& graphics_scene,
                    GSOptions& gso)
{
  if(gso.are_vertices_enabled() && gso.draw_vertex(mesh, dh))
  {
    if (gso.colored_vertex(mesh, dh))
    {
      graphics_scene.add_point(get_point(mesh, dh),
                               gso.vertex_color(mesh, dh));
    }
    else
    { graphics_scene.add_point(get_point(mesh, dh)); }
  }
}

template <typename Mesh, class GSOptions>
void compute_path(const Mesh &mesh,
                  const typename Get_map<Mesh, Mesh>::storage_type& lcc,
                  CGAL::Graphics_scene &graphics_scene,
                  const GSOptions &gso,
                  const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>& paths,
                  std::size_t i,
                  typename Get_map<Mesh, Mesh>::type::size_type amark)
{

  if (paths[i].is_empty() || !gso.draw_path(i))
  { return; }

  CGAL::IO::Color color=gso.color_of_path(i);
  graphics_scene.add_point(get_point(mesh, paths[i].get_ith_dart(0)), color);
  for (std::size_t j=0; j<paths[i].length(); ++j)
  {
    if (!lcc.is_marked(paths[i].get_ith_dart(j), amark))
    {
      compute_edge(mesh, lcc, paths[i].get_ith_dart(j), graphics_scene, gso, true, color);
      lcc.template mark_cell<1>(paths[i].get_ith_dart(j), amark);
    }
  }
}

template <class Mesh, class GSOptions>
void compute_elements(const Mesh &mesh,
                      CGAL::Graphics_scene &graphics_scene,
                      const GSOptions &gso,
                      const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>& paths,
                      typename Get_map<Mesh, Mesh>::type::size_type mark)
{
  typedef typename Get_map<Mesh, Mesh>::type      LCC;
  typedef typename LCC::size_type                 size_type;

  typename Get_map<Mesh, Mesh>::storage_type lcc(mesh);
  size_type oriented_mark = lcc.get_new_mark();
  std::size_t current_path = paths.size();
  typename LCC::size_type amark=mark==(std::numeric_limits<size_type>::max)()?
            LCC::INVALID_MARK:mark; // If !=INVALID_MARK, show darts marked with this mark

  lcc.orient(oriented_mark);

  size_type markfaces    = lcc.get_new_mark();
  size_type markedges    = lcc.get_new_mark();
  size_type markvertices = lcc.get_new_mark();

  if (current_path==paths.size())
  {
    for (std::size_t i=0; i<paths.size(); ++i)
    { compute_path(mesh, lcc, graphics_scene, gso, paths, i, markedges); }
  }
  else if (current_path!=paths.size()+1)
  { compute_path(mesh, lcc, graphics_scene, gso, paths, current_path, markedges); }

  for (typename LCC::Dart_range::const_iterator it=lcc.darts().begin(),
         itend=lcc.darts().end(); it!=itend; ++it )
  {
    if (gso.are_faces_enabled() && !lcc.is_marked(it, markfaces) &&
        !lcc.is_perforated(it) && lcc.is_marked(it, oriented_mark))
    {
      compute_face(mesh, lcc, it, graphics_scene, gso);
      lcc.template mark_cell<2>(it, markfaces);
    }

    if ( !lcc.is_marked(it, markedges) )
    {
      compute_edge(mesh, lcc, it, graphics_scene, gso, amark);
      lcc.template mark_cell<1>(it, markedges);
    }

    if ( !lcc.is_marked(it, markvertices) )
    {
      compute_vertex(mesh, it, graphics_scene, gso);
      lcc.template mark_cell<0>(it, markvertices);
    }
  }

  lcc.free_mark(markfaces);
  lcc.free_mark(markedges);
  lcc.free_mark(markvertices);
}

} // namespace draw_function_for_face_graph_with_paths

template <class Mesh, class GSOptions>
void add_to_graphics_scene(const Mesh& mesh,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& gso,
                           const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>& paths,
                           typename Get_map<Mesh, Mesh>::type::size_type amark=
                           Get_map<Mesh, Mesh>::type::INVALID_MARK)
{
  draw_function_for_face_graph_with_paths::compute_elements(mesh, graphics_scene,
                                                            gso, paths, amark);
}

template <class Mesh>
void add_to_graphics_scene(const Mesh& mesh,
                           CGAL::Graphics_scene& graphics_scene,
                           const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>& paths,
                           typename Get_map<Mesh, Mesh>::type::size_type amark=
                           Get_map<Mesh, Mesh>::type::INVALID_MARK)
{
  // Default graphics view options.
  Graphics_scene_options_face_graph_with_paths<Mesh,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor /*vh*/,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor /*eh*/,
                  typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor /*fh*/>
      gso;

  add_to_graphics_scene(mesh, graphics_scene, gso, paths, amark);
}

template <class Mesh, class GSOptions>
void add_to_graphics_scene(const Mesh& mesh,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& gso,
                           std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>> l,
                           typename Get_map<Mesh, Mesh>::type::size_type amark=
                           Get_map<Mesh, Mesh>::type::INVALID_MARK)
{
  std::vector<Surface_mesh_topology::Path_on_surface<Mesh>> paths=l;
  add_to_graphics_scene(mesh, graphics_scene, gso, paths, amark);
}

template <class Mesh>
void add_to_graphics_scene(const Mesh& mesh,
                           CGAL::Graphics_scene& graphics_scene,
                           std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>> l,
                           typename Get_map<Mesh, Mesh>::type::size_type amark=
                           Get_map<Mesh, Mesh>::type::INVALID_MARK)
{
  std::vector<Surface_mesh_topology::Path_on_surface<Mesh>> paths=l;
  add_to_graphics_scene(mesh, graphics_scene, paths, amark);
}

#ifdef CGAL_USE_BASIC_VIEWER

template<typename Mesh, typename GSOptions>
void draw(const Mesh& mesh,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& paths,
          const GSOptions& gso,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)())
{
  CGAL::Graphics_scene graphics_scene;
  add_to_graphics_scene(mesh, graphics_scene, paths, gso, amark);
  draw_graphics_scene(graphics_scene, title);
}

template<typename Mesh>
void draw(const Mesh& mesh,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& paths,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)())
{
  CGAL::Graphics_scene graphics_scene;
  add_to_graphics_scene(mesh, graphics_scene, paths, amark);
  draw_graphics_scene(graphics_scene, title);
}

template<class Mesh, typename GSOptions>
void draw(const Mesh& mesh,
          std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>> l,
          const GSOptions& gso,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)())
{
  std::vector<Surface_mesh_topology::Path_on_surface<Mesh>> paths=l;
  draw(mesh, paths, gso, title, amark);
}

template<class Mesh>
void draw(const Mesh& mesh,
          std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>> l,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)())
{
  std::vector<Surface_mesh_topology::Path_on_surface<Mesh>> paths=l;
  draw(mesh, paths, title, amark);
}

#endif // CGAL_USE_BASIC_VIEWER

} // End namespace CGAL

#endif // CGAL_DRAW_FACE_GRAPH_WITH_PATHS_H
