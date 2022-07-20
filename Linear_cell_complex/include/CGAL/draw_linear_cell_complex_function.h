// Copyright (c) 2022 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_LCC_H
#define CGAL_DRAW_LCC_H

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/GraphicBuffer.h>

#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Qt/init_ogl_context.h>

#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Random.h>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultDrawingFunctorLCC {
  /// @return true iff the volume containing dh is drawn.
  template <typename LCC>
  bool draw_volume(const LCC &, typename LCC::Dart_const_handle) const {
    return true;
  }
  /// @return true iff the face containing dh is drawn.
  template <typename LCC>
  bool draw_face(const LCC &, typename LCC::Dart_const_handle) const {
    return true;
  }
  /// @return true iff the edge containing dh is drawn.
  template <typename LCC>
  bool draw_edge(const LCC &, typename LCC::Dart_const_handle) const {
    return true;
  }
  /// @return true iff the vertex containing dh is drawn.
  template <typename LCC>
  bool draw_vertex(const LCC &, typename LCC::Dart_const_handle) const {
    return true;
  }

  /// @return true iff the volume containing dh is drawn in wireframe.
  template <typename LCC>
  bool volume_wireframe(const LCC &, typename LCC::Dart_const_handle) const {
    return false;
  }
  /// @return true iff the face containing dh is drawn in wireframe.
  template <typename LCC>
  bool face_wireframe(const LCC &, typename LCC::Dart_const_handle) const {
    return false;
  }

  /// @return true iff the volume containing dh is colored.
  template <typename LCC>
  bool colored_volume(const LCC &, typename LCC::Dart_const_handle) const {
    return true;
  }
  /// @return true iff the face containing dh is colored.
  ///  if we have also colored_volume(alcc, dh), the volume color is
  ///  ignored and only the face color is considered.
  template <typename LCC>
  bool colored_face(const LCC &, typename LCC::Dart_const_handle) const {
    return false;
  }
  /// @return true iff the edge containing dh is colored.
  template <typename LCC>
  bool colored_edge(const LCC &, typename LCC::Dart_const_handle) const {
    return false;
  }
  /// @return true iff the vertex containing dh is colored.
  template <typename LCC>
  bool colored_vertex(const LCC &, typename LCC::Dart_const_handle) const {
    return false;
  }

  /// @return the color of the volume containing dh
  ///  used only if colored_volume(alcc, dh) and !colored_face(alcc, dh)
  template <typename LCC>
  CGAL::IO::Color volume_color(const LCC &alcc,
                               typename LCC::Dart_const_handle dh) const {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the face containing dh
  ///  used only if colored_face(alcc, dh)
  template <typename LCC>
  CGAL::IO::Color face_color(const LCC &alcc,
                             typename LCC::Dart_const_handle dh) const {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the edge containing dh
  ///  used only if colored_edge(alcc, dh)
  template <typename LCC>
  CGAL::IO::Color edge_color(const LCC &alcc,
                             typename LCC::Dart_const_handle dh) const {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the vertex containing dh
  ///  used only if colored_vertex(alcc, dh)
  template <typename LCC>
  CGAL::IO::Color vertex_color(const LCC &alcc,
                               typename LCC::Dart_const_handle dh) const {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
};

template <class LCC, class Local_kernel, int dim = LCC::ambient_dimension>
struct LCC_geom_utils;

template <class LCC, class Local_kernel>
struct LCC_geom_utils<LCC, Local_kernel, 3> {
  static typename Local_kernel::Vector_3
  get_vertex_normal(const LCC &lcc, typename LCC::Dart_const_handle dh) {
    typename Local_kernel::Vector_3 n =
        internal::Geom_utils<typename LCC::Traits, Local_kernel>::
            get_local_vector(CGAL::compute_normal_of_cell_0<LCC>(lcc, dh));
    n = n / (CGAL::sqrt(n * n));
    return n;
  }
};

template <class LCC, class Local_kernel>
struct LCC_geom_utils<LCC, Local_kernel, 2> {
  static typename Local_kernel::Vector_3
  get_vertex_normal(const LCC &, typename LCC::Dart_const_handle) {
    typename Local_kernel::Vector_3 n = CGAL::NULL_VECTOR;
    return n;
  }
};

// #define LCC  CGAL_LCC_TYPE                                

// TODO?


template <typename BufferType = float, class LCC, class Local_kernel, class DrawingFunctorLCC>
void compute_face(typename LCC::Dart_const_handle dh, typename LCC::Dart_const_handle voldh, const LCC *lcc,
                  bool m_nofaces, bool m_random_face_color,
                  const DrawingFunctorLCC &m_drawing_functor,
                  GraphicBuffer<BufferType> &graphic_buffer) {

  typedef typename LCC::Dart_const_handle Dart_const_handle;
  typedef typename LCC::Traits Kernel;
  typedef typename Kernel::Point Point;
  typedef typename Kernel::Vector Vector;

  if (m_nofaces || !m_drawing_functor.draw_face(*lcc, dh))
    return;

  // We fill only closed faces.
  Dart_const_handle cur = dh;
  Dart_const_handle min = dh;
  do {
    if (!lcc->is_next_exist(cur))
      return; // open face=>not filled
    if (cur < min)
      min = cur;
    cur = lcc->next(cur);
  } while (cur != dh);

  if (m_random_face_color) {
    CGAL::Random random((unsigned int)(lcc->darts().index(dh)));
    CGAL::IO::Color c = get_random_color(random);
    graphic_buffer.face_begin(c);
  } else if (m_drawing_functor.colored_face(*lcc, dh)) {
    CGAL::IO::Color c = m_drawing_functor.face_color(*lcc, dh);
    graphic_buffer.face_begin(c);
  } else if (m_drawing_functor.colored_volume(*lcc, voldh)) {
    CGAL::IO::Color c = m_drawing_functor.volume_color(*lcc, voldh);
    graphic_buffer.face_begin(c);
  } else {
    graphic_buffer.face_begin();
  }

  cur = dh;
  do {
    graphic_buffer.add_point_in_face(
        lcc->point(cur),
        LCC_geom_utils<LCC, Local_kernel>::get_vertex_normal(*lcc, cur));
    cur = lcc->next(cur);
  } while (cur != dh);

  graphic_buffer.face_end();
}

template <typename BufferType = float, class LCC, class DrawingFunctorLCC>
void compute_edge(typename LCC::Dart_const_handle dh, const LCC *lcc,
                  const DrawingFunctorLCC &m_drawing_functor,
                  GraphicBuffer<BufferType> &graphic_buffer) {
  typedef typename LCC::Dart_const_handle Dart_const_handle;
  typedef typename LCC::Traits Kernel;
  typedef typename Kernel::Point Point;

  if (!m_drawing_functor.draw_edge(*lcc, dh))
    return;

  Point p1 = lcc->point(dh);
  Dart_const_handle d2 = lcc->other_extremity(dh);
  if (d2 != nullptr) {
    if (m_drawing_functor.colored_edge(*lcc, dh)) {
      graphic_buffer.add_segment(p1, lcc->point(d2),
                                 m_drawing_functor.edge_color(*lcc, dh));
    } else {
      graphic_buffer.add_segment(p1, lcc->point(d2));
    }
  }
}

template <typename BufferType = float, class LCC, class DrawingFunctorLCC>
void compute_vertex(typename LCC::Dart_const_handle dh, const LCC *lcc,
                    const DrawingFunctorLCC &m_drawing_functor,
                    GraphicBuffer<BufferType> &graphic_buffer) {
  if (!m_drawing_functor.draw_vertex(*lcc, dh))
    return;

  if (m_drawing_functor.colored_vertex(*lcc, dh)) {
    graphic_buffer.add_point(lcc->point(dh),
                             m_drawing_functor.vertex_color(*lcc, dh));
  } else {
    graphic_buffer.add_point(lcc->point(dh));
  }
}

template <typename BufferType = float, class LCC, class DrawingFunctorLCC>
void compute_elements(const LCC *lcc,
                      const DrawingFunctorLCC &m_drawing_functor,
                      bool m_nofaces, bool m_random_face_color,
                      GraphicBuffer<BufferType> &graphic_buffer) {

  if (lcc == nullptr)
    return;

  typename LCC::size_type markvolumes = lcc->get_new_mark();
  typename LCC::size_type markfaces = lcc->get_new_mark();
  typename LCC::size_type markedges = lcc->get_new_mark();
  typename LCC::size_type markvertices = lcc->get_new_mark();
  typename LCC::size_type oriented_mark = lcc->get_new_mark();

  lcc->orient(oriented_mark);

  for (typename LCC::Dart_range::const_iterator it = lcc->darts().begin(),
                                                itend = lcc->darts().end();
       it != itend; ++it) {
    if (!lcc->is_marked(it, markvolumes) &&
        m_drawing_functor.draw_volume(*lcc, it)) {
      for (typename LCC::template Dart_of_cell_basic_range<3>::const_iterator
               itv = lcc->template darts_of_cell_basic<3>(it, markvolumes)
                         .begin(),
               itvend =
                   lcc->template darts_of_cell_basic<3>(it, markvolumes).end();
           itv != itvend; ++itv) {
        lcc->mark(itv, markvolumes); // To be sure that all darts of the basic
                                     // iterator will be marked
        if (!lcc->is_marked(itv, markfaces) &&
            lcc->is_marked(itv, oriented_mark) &&
            m_drawing_functor.draw_face(*lcc, itv)) {
          if (!m_drawing_functor.volume_wireframe(*lcc, itv) &&
              !m_drawing_functor.face_wireframe(*lcc, itv)) {
            compute_face(itv, it, lcc, m_nofaces, m_random_face_color,
                         m_drawing_functor, graphic_buffer);
          }
          for (typename LCC::template Dart_of_cell_basic_range<
                   2>::const_iterator
                   itf = lcc->template darts_of_cell_basic<2>(itv, markfaces)
                             .begin(),
                   itfend = lcc->template darts_of_cell_basic<2>(itv, markfaces)
                                .end();
               itf != itfend; ++itf) {
            if (!m_drawing_functor.volume_wireframe(*lcc, itv) &&
                !m_drawing_functor.face_wireframe(*lcc, itv)) {
              lcc->mark(itf, markfaces);
            } // To be sure that all darts of the basic iterator will be marked
            if (!lcc->is_marked(itf, markedges) &&
                m_drawing_functor.draw_edge(*lcc, itf)) {
              compute_edge(itf, lcc, m_drawing_functor, graphic_buffer);
              for (typename LCC::template Dart_of_cell_basic_range<
                       1>::const_iterator
                       ite =
                           lcc->template darts_of_cell_basic<1>(itf, markedges)
                               .begin(),
                       iteend =
                           lcc->template darts_of_cell_basic<1>(itf, markedges)
                               .end();
                   ite != iteend; ++ite) {
                lcc->mark(ite, markedges); // To be sure that all darts of the
                                           // basic iterator will be marked
                if (!lcc->is_marked(ite, markvertices) &&
                    m_drawing_functor.draw_vertex(*lcc, ite)) {
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
                                                itend = lcc->darts().end();
       it != itend; ++it) {
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

// This function is responsible for filling the buffer to allow visualization.
template <typename BufferType = float, class LCC, class DrawingFunctorLCC>
void add_in_graphic_buffer_lcc(GraphicBuffer<BufferType> &graphic_buffer,
                              const DrawingFunctorLCC &m_drawing_functor,
                               const LCC *alcc = nullptr, bool nofill = false,
                               bool m_random_face_color = false
                               ) {

  if (alcc != nullptr) {
    compute_elements(graphic_buffer, alcc, m_drawing_functor);
  }
}

// TODO: Move to Basic_viewer_qt.h
template <typename BufferType = float>
void draw_buffer(GraphicBuffer<BufferType> &graphic_buffer) {

#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite) {
    CGAL::Qt::init_ogl_context(4, 3);
    int argc = 1;
    const char *argv[2] = {"lccviewer", nullptr};
    QApplication app(argc, const_cast<char **>(argv));

    Basic_viewer_qt<float> basic_viewer(app.activeWindow(), graphic_buffer);

    basic_viewer.show();
    app.exec();
  }
}

// Specialization of draw function.
#define CGAL_LCC_TYPE                                                          \
  CGAL::Linear_cell_complex_base<d_, ambient_dim, Traits_, Items_, Alloc_,     \
                                 Map, Refs, Storage_>

template <unsigned int d_, unsigned int ambient_dim, class Traits_,
          class Items_, class Alloc_,
          template <unsigned int, class, class, class, class> class Map,
          class Refs, class Storage_,
          class DrawingFunctorLCC = DefaultDrawingFunctorLCC>
void draw(const CGAL_LCC_TYPE &alcc,
          const char *title = "LCC for CMap Basic Viewer", bool nofill = false,
          const DrawingFunctorLCC &drawing_functor = DrawingFunctorLCC()) {
  GraphicBuffer<float> buffer;
  add_in_graphic_buffer_lcc(buffer, drawing_functor, &alcc, nofill, false);
  draw_buffer(buffer);
}

// Todo a function taking a const DrawingFunctorLCC& drawing_functor as
// parameter
#undef CGAL_LCC_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_LCC_H
