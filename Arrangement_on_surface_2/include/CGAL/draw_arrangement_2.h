#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <map>

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Arrangement_2.h>

namespace CGAL {

// Viewer class for Polygon_2
template <typename Arrangement_2>
class Arr_2_viewer_qt : public Basic_viewer_qt {
  typedef Basic_viewer_qt                       Base;
  typedef Arrangement_2                         Arr;
  typedef typename Arr::Point_2                 Point;
  typedef typename Arr::Face_const_handle       Face_const_handle;

public:
  /// Construct the viewer.
  /// @param arr the arrangement to view
  /// @param title the title of the window
  Arr_2_viewer_qt(QWidget* parent, const Arr& arr,
                  const char* title = "2D Arrangement Basic Viewer") :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    m_arr(arr)
  { add_elements(); }

protected:
  //!
  void add_ccb(typename Arr::Ccb_halfedge_const_circulator circ, bool has_area = true) {
    typename Arr::Ccb_halfedge_const_circulator curr = circ;
    do {
      typename Arr::Halfedge_const_handle e = curr;
      add_segment(e->source()->point(), e->target()->point());
      add_point(e->source()->point());
      if (has_area) add_point_in_face(e->source()->point());
    } while (++curr != circ);
  }

  //!
  void add_face(Face_const_handle face, CGAL::IO::Color c) {
    if (! face->is_unbounded()) {
      face_begin(c);
      add_ccb(face->outer_ccb());
      for (auto iv = face->isolated_vertices_begin();
           iv != face->isolated_vertices_end(); ++iv)
        add_point(iv->point());
      face_end();
    }

    std::map<Face_const_handle, bool> visited;
    for (auto it = face->inner_ccbs_begin(); it != face->inner_ccbs_end(); ++it) {
      auto new_face = (*it)->twin()->face();
      if (new_face == face) {
        add_ccb(*it, false);
        continue;
      }
      if (visited.find(new_face) != visited.end()) continue;
      visited[new_face] = true;
      add_face(new_face, c);
    }
  }

  //!
  void add_elements() {
    clear();
    if (m_arr.is_empty()) return;
    CGAL::IO::Color c(75,160,255);
    for (auto it = m_arr.unbounded_faces_begin(); it != m_arr.unbounded_faces_end(); ++it)
      add_face(it, c);
  }

  //!
  virtual void keyPressEvent(QKeyEvent* e) {
    // Test key pressed:
    //    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    //    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }

    // Call: * add_elements() if the model changed, followed by
    //       * redraw() if some viewing parameters changed that implies some
    //                  modifications of the buffers
    //                  (eg. type of normal, color/mono)
    //       * update() just to update the drawing

    // Call the base method to process others/classicals key
    Base::keyPressEvent(e);
  }

protected:
  const Arrangement_2& m_arr;
};

template <typename GeometryTraits_2, typename Dcel>
void draw(const CGAL::Arrangement_2<GeometryTraits_2, Dcel>& arr,
          const char* title = "2D Arrangement Basic Viewer") {
  typedef GeometryTraits_2              Gt;
  typedef CGAL::Arrangement_2<Gt, Dcel> Arr;

  CGAL::Qt::init_ogl_context(4,3);

  int argc = 1;
  const char* argv[2] = {"t2_viewer", nullptr};
  QApplication app(argc, const_cast<char**>(argv));
  Arr_2_viewer_qt<Arr>mainwindow(app.activeWindow(), arr, title);
  mainwindow.show();
  app.exec();
}

}

#endif
