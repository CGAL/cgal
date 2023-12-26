// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jasmeet Singh    <jasmeet.singh.mec11@iitbhu.ac.in>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef DRAW_NEF_3_H
#define DRAW_NEF_3_H

#include <CGAL/license/Nef_3.h>
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/circulator.h>
#include <CGAL/Random.h>
#include <CGAL/assertions.h>

#include <unordered_map>

namespace CGAL {

namespace draw_function_for_nef_polyhedron
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Vector_3 Local_vector;

template <class Nef_Polyhedron>
Local_vector get_face_normal(typename Nef_Polyhedron::SHalfedge_const_handle she)
{
  typename Nef_Polyhedron::SHalfedge_around_facet_const_circulator he(she), end(he);
  Local_vector normal=CGAL::NULL_VECTOR;
  unsigned int nb=0;

  using GU=internal::Geom_utils<typename Kernel_traits
                                <typename Nef_Polyhedron::Point_3>::Kernel,
                                Local_kernel>;
  CGAL_For_all(he, end)
  {
    internal::newell_single_step_3(GU::get_local_point(he->next()->source()->
                                                       center_vertex()->point()),
                                   GU::get_local_point(he->source()->center_vertex()->
                                                       point()), normal);
    ++nb;
  }

  CGAL_assertion(nb > 0);
  return (typename Local_kernel::Construct_scaled_vector_3()(normal, 1.0/nb));
}

template <class Nef_Polyhedron>
Local_vector get_vertex_normal(typename Nef_Polyhedron::Vertex_const_handle vh)
{
  typename Nef_Polyhedron::SHalfedge_const_iterator it=vh->shalfedges_begin();
  typename Nef_Polyhedron::SHalfedge_const_handle end=it;
  Local_vector normal=CGAL::NULL_VECTOR;
  do
  {
    Local_vector n=get_face_normal<Nef_Polyhedron>(it);
    normal=typename Local_kernel::Construct_sum_of_vectors_3()(normal, n);
    it=it->snext();
  }
  while(it!=end);

  if (!typename Local_kernel::Equal_3()(normal, CGAL::NULL_VECTOR))
  {
    normal=(typename Local_kernel::Construct_scaled_vector_3()
            (normal, 1.0/CGAL::sqrt(normal.squared_length())));
  }

  return normal;
}

// Visitor class to iterate through shell objects
template <typename Nef_Polyhedron, typename GSOptions>
class Nef_Visitor {

  typedef typename Nef_Polyhedron::Halfedge_const_handle     Halfedge_const_handle;
  typedef typename Nef_Polyhedron::Halffacet_const_handle    Halffacet_const_handle;
  typedef typename Nef_Polyhedron::SHalfedge_const_handle    SHalfedge_const_handle;
  typedef typename Nef_Polyhedron::SHalfedge_around_facet_const_circulator   SHalfedge_around_facet_const_circulator;
  typedef typename Nef_Polyhedron::Vertex_const_handle       Vertex_const_handle;
  typedef typename Nef_Polyhedron::SFace_const_handle        SFace_const_handle;
  typedef typename Nef_Polyhedron::SHalfloop_const_handle    SHalfloop_const_handle;
  typedef typename Nef_Polyhedron::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;

public:
  Nef_Visitor(const Nef_Polyhedron&_nef,
              CGAL::Graphics_scene& _graphics_scene,
              const GSOptions&_gs_options) :
    n_faces(0), n_edges(0),
    graphics_scene(_graphics_scene),
    gs_options(_gs_options),
    nef(_nef)
  {}

  void visit(Vertex_const_handle vh)
  {
    if (!gs_options.are_vertices_enabled() ||
        !gs_options.draw_vertex(nef, vh))
    { return; }

    if(gs_options.colored_vertex(nef, vh))
    { graphics_scene.add_point(vh->point(),
                               gs_options.vertex_color(nef, vh)); }
    else
    { graphics_scene.add_point(vh->point()); }
  }

  void visit(Halffacet_const_handle opposite_facet)
  {
    if (!gs_options.are_faces_enabled() ||
        !gs_options.draw_face(nef, opposite_facet))
    { return; }

    Halffacet_const_handle f=opposite_facet->twin();
    if (facets_done.find(f)!=facets_done.end() ||
        facets_done.find(opposite_facet)!=facets_done.end())
    { return; }

    SHalfedge_const_handle se;
    Halffacet_cycle_const_iterator fc = f->facet_cycles_begin();

    se = SHalfedge_const_handle(fc); // non-zero if shalfedge is returned
    if(se==0)
    { return; } //return if not-shalfedge

    if(gs_options.colored_face(nef, f))
    { graphics_scene.face_begin(gs_options.face_color(nef, f)); }
    else
    { graphics_scene.face_begin(); }

    SHalfedge_around_facet_const_circulator hc_start(se);
    SHalfedge_around_facet_const_circulator hc_end(hc_start);
    Vertex_const_handle lastvh;
    CGAL_For_all(hc_start, hc_end)
    {
      Vertex_const_handle vh=hc_start->source()->center_vertex();
      lastvh=vh;
      graphics_scene.add_point_in_face(vh->point(),
                                       get_vertex_normal<Nef_Polyhedron>(vh));
    }

    // Now iterate through holes of the face
    ++fc;
    while(fc!=f->facet_cycles_end())
    {
      if(fc.is_shalfedge())
      {
        se = SHalfedge_const_handle(fc);
        hc_start=se;
        hc_end=hc_start;
        CGAL_For_all(hc_start, hc_end)
        {
          Vertex_const_handle vh=hc_start->source()->center_vertex();
          graphics_scene.add_point_in_face(vh->point(),
                                            get_vertex_normal<Nef_Polyhedron>(vh));
        }
        graphics_scene.add_point_in_face(hc_start->source()->center_vertex()->point(),
                                          get_vertex_normal<Nef_Polyhedron>
                                          (hc_start->source()->center_vertex()));
        graphics_scene.add_point_in_face(lastvh->point(),
                                          get_vertex_normal<Nef_Polyhedron>(lastvh));
      }
      ++fc;
    }

    graphics_scene.face_end();
    facets_done[f]=true;
    ++n_faces;
  }

  void visit(Halfedge_const_handle he)
  {
    if (!gs_options.are_edges_enabled() ||
        !gs_options.draw_edge(nef, he))
    { return; }

    Halfedge_const_handle twin=he->twin();
    if (edges_done.find(he)!=edges_done.end() ||
        edges_done.find(twin)!=edges_done.end())
    { return; } // Edge already added

    if(gs_options.colored_edge(nef, he))
    { graphics_scene.add_segment(he->source()->point(), he->target()->point(),
                                 gs_options.edge_color(nef, he)); }
    else
    { graphics_scene.add_segment(he->source()->point(), he->target()->point()); }
    edges_done[he]=true;
    ++n_edges;
  }

  void visit(SHalfedge_const_handle ) {}
  void visit(SHalfloop_const_handle ) {}
  void visit(SFace_const_handle ) {}
  int n_faces;
  int n_edges;
protected:
  std::unordered_map<Halffacet_const_handle, bool> facets_done;
  std::unordered_map<Halfedge_const_handle, bool> edges_done;
  CGAL::Graphics_scene& graphics_scene;
  const GSOptions& gs_options;
  const Nef_Polyhedron& nef;
};

template <class Nef_Polyhedron, class GSOptions>
void compute_elements(const Nef_Polyhedron &nef,
                      CGAL::Graphics_scene &graphics_scene,
                      const GSOptions &gs_options)
{

  typedef typename Nef_Polyhedron::Volume_const_iterator      Volume_const_iterator;
  typedef typename Nef_Polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;
  typedef typename Nef_Polyhedron::SFace_const_handle         SFace_const_handle;

  Volume_const_iterator c;
  Nef_Visitor<Nef_Polyhedron, GSOptions> V(nef, graphics_scene, gs_options);
  CGAL_forall_volumes(c, nef)
  {
    Shell_entry_const_iterator it;
    CGAL_forall_shells_of(it, c)
    { nef.visit_shell_objects(SFace_const_handle(it), V); }
  }

  graphics_scene.reverse_all_normals();
}

} // namespace draw_function_for_nef_polyhedron

#define CGAL_NEF3_TYPE Nef_polyhedron_3<Kernel_, Items_, Mark_>

// add_to_graphics_scene
template <typename Kernel_, typename Items_, typename Mark_,
          class GSOptions>
void add_to_graphics_scene(const CGAL_NEF3_TYPE &anef,
                           CGAL::Graphics_scene &graphics_scene,
                           const GSOptions &gs_options)
{
  draw_function_for_nef_polyhedron::compute_elements(anef,
                                                     graphics_scene,
                                                     gs_options);
}

template <typename Kernel_, typename Items_, typename Mark_>
void add_to_graphics_scene(const CGAL_NEF3_TYPE &anef,
                           CGAL::Graphics_scene &graphics_scene)
{
  // Default graphics view options.
  Graphics_scene_options<CGAL_NEF3_TYPE,
                  typename CGAL_NEF3_TYPE::Vertex_const_handle /*vh*/,
                  typename CGAL_NEF3_TYPE::Halfedge_const_handle /*eh*/,
                  typename CGAL_NEF3_TYPE::Halffacet_const_handle /*fh*/>
      gs_options;

  gs_options.colored_face=[](const CGAL_NEF3_TYPE&,
                             typename CGAL_NEF3_TYPE::Halffacet_const_handle) -> bool
  { return true; };

  gs_options.face_color=[](const CGAL_NEF3_TYPE&,
                           typename CGAL_NEF3_TYPE::Halffacet_const_handle fh) -> CGAL::IO::Color
  {
    if (fh==nullptr) // use to get the mono color
    { return CGAL::IO::Color(100, 125, 200); }

    CGAL::Random random((unsigned int)(std::size_t)(&(*fh)));
    return get_random_color(random);
  };

  add_to_graphics_scene(anef, graphics_scene, gs_options);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function
template <typename Kernel_, typename Items_, typename Mark_,
          class GSOptions>
void draw(const CGAL_NEF3_TYPE &anef,
          const GSOptions &gs_options,
          const char *title="Nef Polyhedron Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(anef, buffer, gs_options);
  draw_graphics_scene(buffer, title);
}

template <typename Kernel_, typename Items_, typename Mark_>
void draw(const CGAL_NEF3_TYPE &anef,
          const char *title="Nef Polyhedron Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(anef, buffer);
  draw_graphics_scene(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_NEF3_TYPE

} // End namespace CGAL

#endif // DRAW_NEF_3_H
