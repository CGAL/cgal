// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jasmeet Singh    <jasmeet.singh.mec11@iitbhu.ac.in>

#ifndef DRAW_NEF_3_H
#define DRAW_NEF_3_H

#include <CGAL/license/Nef_3.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/circulator.h>
#include <CGAL/Random.h>
#include <unordered_map>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorNefPolyhedron
{
  template<typename NefPolyhedron>
  static CGAL::Color run(const NefPolyhedron&,
                         typename NefPolyhedron::Halffacet_const_handle fh)
  {
    if (fh==boost::graph_traits<NefPolyhedron>::null_face()) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    CGAL::Random random((unsigned int)(std::size_t)(&(*fh)));
    return get_random_color(random);
  }
};

// Viewer class for Nef Polyhedron
template<class Nef_Polyhedron, class ColorFunctor>
class SimpleNefPolyhedronViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt                                   Base;
  typedef typename Nef_Polyhedron::Kernel                   Kernel;
  typedef typename Nef_Polyhedron::Vertex_const_iterator    Vertex_const_iterator;
  typedef typename Nef_Polyhedron::Vertex_const_handle      Vertex_const_handle;
  typedef typename Nef_Polyhedron::SNC_const_decorator      SNC_const_decorator;
  typedef typename Nef_Polyhedron::SNC_structure            SNC_structure;
  typedef typename SNC_structure::Infi_box                     Infi_box;

  typedef typename SNC_structure::Halffacet_cycle_const_iterator
      Halffacet_cycle_const_iterator;

  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;

  typedef typename Nef_Polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;
  typedef typename Nef_Polyhedron::SFace_const_handle SFace_const_handle;
  typedef typename Nef_Polyhedron::Volume_const_iterator Volume_const_iterator;

  typedef typename Nef_Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Nef_Polyhedron::Halffacet_const_handle Halffacet_const_handle;
  typedef typename Nef_Polyhedron::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Nef_Polyhedron::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename Kernel::Point_3 Point_3;

public:
  /// Construct the viewer
  /// @param anef the nef polyhedron to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not
  ///        computed: this can be useful for big objects)
  SimpleNefPolyhedronViewerQt(QWidget* parent,
                              const Nef_Polyhedron& anef,
                              const char* title="Basic Nef Polyhedron Viewer",
                              bool anofaces=false,
                              const ColorFunctor& fcolor=ColorFunctor()) :
  //First draw: no vertex; edges, faces; mon-color; inverse normal
    Base(parent, title, true, true, true, true, false),
    nef(anef),
    m_nofaces(anofaces),
    m_fcolor(fcolor)
  {
    compute_elements();
  }
protected:
//  void compute_face(Facet_const_handle fh)
//  {
//    CGAL::Color c=m_fcolor.run(poly, fh);
//    face_begin(c);
//    Halfedge_const_handle he=fh->facet_begin();
//    do
//    {
//      add_point_in_face(he->vertex()->point(),
//                        get_vertex_normal(he));
//      he=he->next();
//    }
//    while (he!=fh->facet_begin());
//    face_end();
//  }

//  void compute_edge(Halfedge_const_handle he)
//  {
//    add_segment(he->vertex()->point(),
//                he->opposite()->vertex()->point());
//    // We can use add_segment(p1, p2, c) with c a CGAL::Color to add a colored segment
//  }
  class Nef_Visitor {
  public:
    Nef_Visitor(SimpleNefPolyhedronViewerQt &v)
        : viewer(v), done(), n_faces(0) {}

    void visit(Vertex_const_handle vh) {
      viewer.add_point(vh->point());
    }

    void visit(Halffacet_const_handle opposite_facet)
    {
      Halffacet_const_handle f = opposite_facet->twin();

      if (done.find(f) != done.end() ||
          done.find(opposite_facet) != done.end()) {
        return;
      }

      SHalfedge_const_handle se;
      Halffacet_cycle_const_iterator fc;
      fc = f->facet_cycles_begin();

      se = SHalfedge_const_handle(fc); // non-zero if shalfedge is returned
      if(se == 0)
      { //return if not-shalfedge
        return;
      }

      viewer.face_begin(CGAL::Color(yellow()));

      SHalfedge_around_facet_const_circulator hc_start(se);
      SHalfedge_around_facet_const_circulator hc_end(hc_start);
      CGAL_For_all(hc_start, hc_end) {
        viewer.add_point_in_face(hc_start->source()->center_vertex()->point(),
                                 viewer.get_vertex_normal(hc_start));
      }
      viewer.face_end();
      done[f] = true;
      n_faces++;
    }

    void visit(Halfedge_const_handle ) {}
    void visit(SHalfedge_const_handle ) {}
    void visit(SHalfloop_const_handle ) {}
    void visit(SFace_const_handle ) {}
    int n_faces;
  protected:
    std::unordered_map<Halffacet_const_handle, bool> done;
    SimpleNefPolyhedronViewerQt& viewer;
  };

  void compute_vertex(Vertex_const_handle vh)
  {
    add_point(vh->point());
    // We can use add_point(p, c) with c a CGAL::Color to add a colored point
  }

  void compute_elements()
  {
    clear();

    Volume_const_iterator c;

    int ic = 0;
    Nef_Visitor V(*this);
    CGAL_forall_volumes(c, nef)
    {
      std::cout << "Volume " << ic++ << std::endl;
      int is = 0;
      Shell_entry_const_iterator it;
      CGAL_forall_shells_of(it, c)
      {
        nef.visit_shell_objects(SFace_const_handle(it), V);
      }
      std::cout << "Faces drawn: " << V.n_faces << std::endl;
    }
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    // Test key pressed:
    //    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    //    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }

    // Call: * compute_elements() if the model changed, followed by
    //       * redraw() if some viewing parameters changed that implies some
    //                  modifications of the buffers
    //                  (eg. type of normal, color/mono)
    //       * update() just to update the drawing

    // Call the base method to process others/classicals key
    Base::keyPressEvent(e);
  }

protected:
  Local_vector get_face_normal(SHalfedge_around_facet_const_circulator he)
  {
    Local_vector normal = CGAL::NULL_VECTOR;
    SHalfedge_around_facet_const_circulator end = he;
    unsigned int nb = 0;

    CGAL_For_all(he, end)
    {
      internal::newell_single_step_3(internal::Geom_utils<Kernel>::get_local_point(he->source()->center_vertex()->point()),
                                     internal::Geom_utils<Kernel>::get_local_point(he->next()->source()->center_vertex()->point()),
                                     normal);
      ++nb;
    }

    assert(nb > 0);
    return (typename Local_kernel::Construct_scaled_vector_3()(normal, 1.0 / nb));
  }

  Local_vector get_vertex_normal(SHalfedge_around_facet_const_circulator he)
  {
    Local_vector normal = CGAL::NULL_VECTOR;
    SHalfedge_around_facet_const_circulator end = he;

    CGAL_For_all(he, end)
    {
      Local_vector n = get_face_normal(he);
      normal = typename Local_kernel::Construct_sum_of_vectors_3()(normal, n);
    }

    if (!typename Local_kernel::Equal_3()(normal, CGAL::NULL_VECTOR))
    {
      normal = (typename Local_kernel::Construct_scaled_vector_3()(
          normal, 1.0 / CGAL::sqrt(normal.squared_length())));
    }

    return normal;
  }

protected:
  const Nef_Polyhedron &nef;
  bool m_nofaces;
  const ColorFunctor &m_fcolor;
};

template <class Nef_Polyhedron, class ColorFunctor>
void draw(const Nef_Polyhedron &anef, const char *title, bool nofill,
          const ColorFunctor &fcolor)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = false;
#endif

  if (!cgal_test_suite)
  {
    int argc = 1;
    const char *argv[2] = {"nef_polyhedron_viewer", "\0"};
    QApplication app(argc, const_cast<char **>(argv));
    SimpleNefPolyhedronViewerQt<Nef_Polyhedron, ColorFunctor> mainwindow(
        app.activeWindow(), anef, title, nofill, fcolor);
    mainwindow.show();
    app.exec();
  }
}

template <class Nef_Polyhedron>
void draw(const Nef_Polyhedron &anef, const char *title, bool nofill)
{
  DefaultColorFunctorNefPolyhedron c;
  draw(anef, title, nofill, c);
}

template <class Nef_Polyhedron>
void draw(const Nef_Polyhedron &anef, const char *title)
{
  draw(anef, title, false);
}

template <class Nef_Polyhedron> void draw(const Nef_Polyhedron &anef)
{
  draw(anef, "Basic Nef Polyhedron Viewer");
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DRAW_NEF_3_H
