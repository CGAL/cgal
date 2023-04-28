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
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Graphic_storage.h>
#include <CGAL/Drawing_functor.h>
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
template <typename Nef_Polyhedron, typename DrawingFunctor, typename BufferType = float>
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
              CGAL::Graphic_storage<BufferType>& _graphic_buffer,
              const DrawingFunctor&_drawing_functor) :
    n_faces(0), n_edges(0),
    nef(_nef),
    graphic_buffer(_graphic_buffer),
    drawing_functor(_drawing_functor)
  {}

  void visit(Vertex_const_handle vh)
  {
    if (!drawing_functor.are_vertices_enabled() ||
        !drawing_functor.draw_vertex(nef, vh))
    { return; }

    if(drawing_functor.colored_vertex(nef, vh))
    { graphic_buffer.add_point(vh->point(),
                               drawing_functor.vertex_color(nef, vh)); }
    else
    { graphic_buffer.add_point(vh->point()); }
  }

  void visit(Halffacet_const_handle opposite_facet)
  {
    if (!drawing_functor.are_faces_enabled() ||
        !drawing_functor.draw_face(nef, opposite_facet))
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

    if(drawing_functor.colored_face(nef, f))
    { graphic_buffer.face_begin(drawing_functor.face_color(nef, f)); }
    else
    { graphic_buffer.face_begin(); }

    SHalfedge_around_facet_const_circulator hc_start(se);
    SHalfedge_around_facet_const_circulator hc_end(hc_start);
    Vertex_const_handle lastvh;
    CGAL_For_all(hc_start, hc_end)
    {
      Vertex_const_handle vh=hc_start->source()->center_vertex();
      lastvh=vh;
      graphic_buffer.add_point_in_face(vh->point(),
                                       get_vertex_normal<Nef_Polyhedron>(vh));
    }

    // Now iterate through holes of the face
    ++fc;
    while(fc!=f->facet_cycles_end())
    {
      se = SHalfedge_const_handle(fc);
      hc_start=se;
      hc_end=hc_start;
      CGAL_For_all(hc_start, hc_end)
      {
        Vertex_const_handle vh=hc_start->source()->center_vertex();
        graphic_buffer.add_point_in_face(vh->point(),
                                         get_vertex_normal<Nef_Polyhedron>(vh));
      }
      graphic_buffer.add_point_in_face(hc_start->source()->center_vertex()->point(),
                                       get_vertex_normal<Nef_Polyhedron>
                                       (hc_start->source()->center_vertex()));
      graphic_buffer.add_point_in_face(lastvh->point(),
                                       get_vertex_normal<Nef_Polyhedron>(lastvh));
      ++fc;
    }

    graphic_buffer.face_end();
    facets_done[f]=true;
    ++n_faces;
  }

  void visit(Halfedge_const_handle he)
  {
    if (!drawing_functor.are_edges_enabled() ||
        !drawing_functor.draw_edge(nef, he))
    { return; }

    Halfedge_const_handle twin=he->twin();
    if (edges_done.find(he)!=edges_done.end() ||
        edges_done.find(twin)!=edges_done.end())
    { return; } // Edge already added

    if(drawing_functor.colored_edge(nef, he))
    { graphic_buffer.add_segment(he->source()->point(), he->target()->point(),
                                 drawing_functor.edge_color(nef, he)); }
    else
    { graphic_buffer.add_segment(he->source()->point(), he->target()->point()); }
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
  CGAL::Graphic_storage<BufferType>& graphic_buffer;
  const DrawingFunctor& drawing_functor;
  const Nef_Polyhedron& nef;
};

template <typename BufferType=float, class Nef_Polyhedron, class DrawingFunctor>
void compute_elements(const Nef_Polyhedron &nef,
                      CGAL::Graphic_storage<BufferType> &graphic_buffer,
                      const DrawingFunctor &drawing_functor)
{

  typedef typename Nef_Polyhedron::Volume_const_iterator      Volume_const_iterator;
  typedef typename Nef_Polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;
  typedef typename Nef_Polyhedron::SFace_const_handle         SFace_const_handle;

  Volume_const_iterator c;
  Nef_Visitor<Nef_Polyhedron, DrawingFunctor, BufferType>
    V(nef, graphic_buffer, drawing_functor);
  CGAL_forall_volumes(c, nef)
  {
    Shell_entry_const_iterator it;
    CGAL_forall_shells_of(it, c)
    { nef.visit_shell_objects(SFace_const_handle(it), V); }
  }

  graphic_buffer.negate_all_normals();
}

} // namespace draw_function_for_nef_polyhedron

#define CGAL_NEF3_TYPE Nef_polyhedron_3<Kernel_, Items_, Mark_>

// add_in_graphic_buffer
template <typename Kernel_, typename Items_, typename Mark_,
          typename BufferType=float, class DrawingFunctor>
void add_in_graphic_buffer(const CGAL_NEF3_TYPE &anef,
                           CGAL::Graphic_storage<BufferType> &graphic_buffer,
                           const DrawingFunctor &drawing_functor)
{
  draw_function_for_nef_polyhedron::compute_elements(anef,
                                                     graphic_buffer,
                                                     drawing_functor);
}

template <typename Kernel_, typename Items_, typename Mark_,
          typename BufferType=float>
void add_in_graphic_buffer(const CGAL_NEF3_TYPE &anef,
                           CGAL::Graphic_storage<BufferType> &graphic_buffer)
{
  // Default functor; user can add his own functor.
  Drawing_functor<CGAL_NEF3_TYPE,
                  typename CGAL_NEF3_TYPE::Vertex_const_handle /*vh*/,
                  typename CGAL_NEF3_TYPE::Halfedge_const_handle /*eh*/,
                  typename CGAL_NEF3_TYPE::Halffacet_const_handle /*fh*/>
      drawing_functor;

  drawing_functor.colored_face = [] (const CGAL_NEF3_TYPE&,
                                     typename CGAL_NEF3_TYPE::Halffacet_const_handle) -> bool
  { return true; };

  drawing_functor.face_color = [] (const CGAL_NEF3_TYPE&,
                                   typename CGAL_NEF3_TYPE::Halffacet_const_handle fh) -> CGAL::IO::Color
  {
    if (fh==nullptr) // use to get the mono color
    { return CGAL::IO::Color(100, 125, 200); }

    CGAL::Random random((unsigned int)(std::size_t)(&(*fh)));
    return get_random_color(random);
  };

  add_in_graphic_buffer(anef, graphic_buffer, drawing_functor);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function
template <typename Kernel_, typename Items_, typename Mark_,
          typename BufferType=float, class DrawingFunctor>
void draw(const CGAL_NEF3_TYPE &anef,
          const DrawingFunctor &drawing_functor,
          const char *title="Nef Polyhedron Viewer")
{
  CGAL::Graphic_storage<BufferType> buffer;
  add_in_graphic_buffer(anef, buffer, drawing_functor);
  draw_graphic_storage(buffer, title);
}

template <typename Kernel_, typename Items_, typename Mark_,
          typename BufferType=float>
void draw(const CGAL_NEF3_TYPE &anef,
          const char *title="Nef Polyhedron Viewer")
{
  CGAL::Graphic_storage<BufferType> buffer;
  add_in_graphic_buffer(anef, buffer);
  draw_graphic_storage(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

} // End namespace CGAL

#endif // DRAW_NEF_3_H
