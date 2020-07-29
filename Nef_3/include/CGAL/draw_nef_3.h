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
    if (fh == nullptr) // use to get the mono color
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

  typedef typename Nef_Polyhedron::Halffacet_cycle_const_iterator            Halffacet_cycle_const_iterator;
  typedef typename Nef_Polyhedron::SHalfedge_around_facet_const_circulator   SHalfedge_around_facet_const_circulator;

  typedef typename Nef_Polyhedron::Shell_entry_const_iterator   Shell_entry_const_iterator;
  typedef typename Nef_Polyhedron::SHalfedge_const_iterator     SHalfedge_const_iterator;
  typedef typename Nef_Polyhedron::Volume_const_iterator        Volume_const_iterator;

  typedef typename Nef_Polyhedron::Vertex_const_handle       Vertex_const_handle;
  typedef typename Nef_Polyhedron::SFace_const_handle        SFace_const_handle;
  typedef typename Nef_Polyhedron::Halfedge_const_handle     Halfedge_const_handle;
  typedef typename Nef_Polyhedron::Halffacet_const_handle    Halffacet_const_handle;
  typedef typename Nef_Polyhedron::SHalfedge_const_handle    SHalfedge_const_handle;
  typedef typename Nef_Polyhedron::SHalfloop_const_handle    SHalfloop_const_handle;

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
  //First draw: vertex; edges, faces; mon-color; inverse normal
    Base(parent, title, false, true, true, true, false),
    nef(anef),
    m_nofaces(anofaces),
    m_fcolor(fcolor)
  {
    compute_elements();
  }
protected:
  // Visitor class to iterate through shell objects
  class Nef_Visitor {
  public:
    Nef_Visitor(SimpleNefPolyhedronViewerQt &v)
      : n_faces(0), n_edges(0), viewer(v) {}

    void visit(Vertex_const_handle vh) {
      viewer.add_point(vh->point());
    }

    void visit(Halffacet_const_handle opposite_facet)
    {
      Halffacet_const_handle f = opposite_facet->twin();

      if (facets_done.find(f) != facets_done.end() ||
          facets_done.find(opposite_facet) != facets_done.end()) {
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

      CGAL::Color c = viewer.run_color(f);
      viewer.face_begin(c);

      SHalfedge_around_facet_const_circulator hc_start(se);
      SHalfedge_around_facet_const_circulator hc_end(hc_start);
      CGAL_For_all(hc_start, hc_end) {
        Vertex_const_handle vh = hc_start->source()->center_vertex();
        viewer.add_point_in_face(vh->point(),
                                 viewer.get_vertex_normal(vh));
      }
      viewer.face_end();
      facets_done[f] = true;
      n_faces++;
    }

    void visit(Halfedge_const_handle he)
    {
      Halfedge_const_handle twin = he->twin();
      if (edges_done.find(he) != edges_done.end() ||
          edges_done.find(twin) != edges_done.end())
      {
        // Edge already added
        return;
      }

      viewer.add_segment(he->source()->point(), he->target()->point());
      edges_done[he] = true;
      n_edges++;
    }

    void visit(SHalfedge_const_handle ) {}
    void visit(SHalfloop_const_handle ) {}
    void visit(SFace_const_handle ) {}
    int n_faces;
    int n_edges;
  protected:
    std::unordered_map<Halffacet_const_handle, bool> facets_done;
    std::unordered_map<Halfedge_const_handle, bool> edges_done;
    SimpleNefPolyhedronViewerQt& viewer;
  };

  void compute_elements()
  {
    clear();

    Volume_const_iterator c;

    Nef_Visitor V(*this);
    CGAL_forall_volumes(c, nef)
    {
      Shell_entry_const_iterator it;
      CGAL_forall_shells_of(it, c)
      {
        nef.visit_shell_objects(SFace_const_handle(it), V);
      }
    }

    negate_all_normals();
  }

  CGAL::Color run_color(Halffacet_const_handle fh)
  {
    return m_fcolor.run(nef, fh);
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
  Local_vector get_face_normal(SHalfedge_const_handle she)
  {
    SHalfedge_around_facet_const_circulator he(she);
    Local_vector normal = CGAL::NULL_VECTOR;
    SHalfedge_around_facet_const_circulator end = he;
    unsigned int nb = 0;

    CGAL_For_all(he, end)
    {
      internal::newell_single_step_3(this->get_local_point
                                     (he->next()->source()->center_vertex()->point()),
                                     this->get_local_point(he->source()->center_vertex()->
                                                           point()), normal);
      ++nb;
    }

    assert(nb > 0);
    return (typename Local_kernel::Construct_scaled_vector_3()(normal, 1.0 / nb));
  }

  Local_vector get_vertex_normal(Vertex_const_handle vh)
  {
    Local_vector normal = CGAL::NULL_VECTOR;

    SHalfedge_const_iterator it = vh->shalfedges_begin();
    SHalfedge_const_handle end = it;
    do {
      Local_vector n = get_face_normal(it);
      normal = typename Local_kernel::Construct_sum_of_vectors_3()(normal, n);
      it = it->snext();
    } while( it != end );

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

#define CGAL_NEF3_TYPE Nef_polyhedron_3<Kernel_, Items_, Mark_>

template <typename Kernel_, typename Items_, typename Mark_>
void draw(const CGAL_NEF3_TYPE &anef,
          const char *title = "Nef Polyhedron Viewer",
          bool nofill = false) {
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    int argc = 1;
    const char *argv[2] = {"nef_polyhedron_viewer", "\0"};
    QApplication app(argc, const_cast<char **>(argv));
    DefaultColorFunctorNefPolyhedron fcolor;
    SimpleNefPolyhedronViewerQt<CGAL_NEF3_TYPE,
                                DefaultColorFunctorNefPolyhedron>
        mainwindow(app.activeWindow(), anef, title, nofill, fcolor);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DRAW_NEF_3_H
