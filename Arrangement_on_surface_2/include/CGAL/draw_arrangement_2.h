#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <unordered_map>

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_conic_traits_2.h>

namespace CGAL {

// Viewer class for Polygon_2
template <typename GeometryTraits_2, typename Dcel>
class Arr_2_basic_viewer_qt : public Basic_viewer_qt {
  typedef GeometryTraits_2                      Gt;
  typedef CGAL::Arrangement_2<Gt, Dcel>         Arr;
  typedef Basic_viewer_qt                       Base;
  typedef typename Arr::Point_2                 Point;
  typedef typename Arr::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arr::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arr::Face_const_handle       Face_const_handle;
  typedef typename Arr::Ccb_halfedge_const_circulator
                                                Ccb_halfedge_const_circulator;

public:
  /// Construct the viewer.
  /// @param arr the arrangement to view
  /// @param title the title of the window
  Arr_2_basic_viewer_qt(QWidget* parent, const Arr& arr,
                        const char* title = "2D Arrangement Basic Viewer") :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    m_arr(arr)
  {}

  //!
  void add_elements() {
    clear();
    if (m_arr.is_empty()) return;
    CGAL::IO::Color c(75,160,255);
    for (auto it = m_arr.unbounded_faces_begin();
         it != m_arr.unbounded_faces_end(); ++it)
      add_face(it, c);

    // Add edges that do not separe faces.
    for (auto it = m_arr.edges_begin(); it != m_arr.edges_end(); ++it)
      if (it->face() == it->twin()->face()) add_edge(it);

    // Add all points
    for (auto it = m_arr.vertices_begin(); it != m_arr.vertices_end(); ++it)
      add_vertex(it);
  }

protected:
  //!
  virtual void add_ccb(Ccb_halfedge_const_circulator circ) {
    auto curr = circ;
    do {
      this->add_point_in_face(curr->source()->point());
      add_edge(curr);
    } while (++curr != circ);
  }

  //!
  virtual void add_edge(Halfedge_const_handle e)
  { this->add_segment(e->source()->point(), e->target()->point()); }

  //!
  virtual void add_vertex(Vertex_const_handle v)
  { this->add_point(v->point()); }

  //!
  void add_face(Face_const_handle face, CGAL::IO::Color color) {
    if (! face->is_unbounded()) {
      face_begin(color);
      add_ccb(face->outer_ccb());
      for (auto iv = face->isolated_vertices_begin();
           iv != face->isolated_vertices_end(); ++iv)
        add_point(iv->point());
      face_end();
    }

    for (auto it = face->inner_ccbs_begin(); it != face->inner_ccbs_end(); ++it)
    {
      std::unordered_map<Face_const_handle, bool> visited;
      auto curr = *it;
      do {
        auto new_face = curr->twin()->face();
        if (new_face == face) continue;
        if (visited.find(new_face) == visited.end()) {
          visited[new_face] = true;
          add_face(new_face, color);
        }
      } while (++curr != *it);
    }
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
  //! The arrangement to draw
  const Arr& m_arr;
};

//! Basic viewer of a 2D arrangement.
template <typename GeometryTraits_2, typename Dcel>
class Arr_2_viewer_qt : public Arr_2_basic_viewer_qt<GeometryTraits_2, Dcel> {
public:
  typedef GeometryTraits_2                      Gt;
  typedef CGAL::Arrangement_2<Gt, Dcel>         Arr;
  typedef Arr_2_basic_viewer_qt<Gt, Dcel>       Base;
  typedef typename Arr::Point_2                 Point;
  typedef typename Arr::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arr::Face_const_handle       Face_const_handle;
  typedef typename Arr::Ccb_halfedge_const_circulator
                                                Ccb_halfedge_const_circulator;

  /// Construct the viewer.
  /// @param arr the arrangement to view
  /// @param title the title of the window
  Arr_2_viewer_qt(QWidget* parent, const Arr& arr,
                  const char* title = "2D Arrangement Basic Viewer") :
    Base(parent, arr, title)
  {}
};

//!
template <typename RatKernel, typename AlgKernel, typename NtTraits,
          typename Dcel>
class Arr_2_viewer_qt<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>,
                      Dcel> :
    public Arr_2_basic_viewer_qt<Arr_conic_traits_2<RatKernel, AlgKernel,
                                                    NtTraits>, Dcel> {
public:
  typedef Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>    Gt;
  typedef CGAL::Arrangement_2<Gt, Dcel>                         Arr;
  typedef Arr_2_basic_viewer_qt<Gt, Dcel>                       Base;
  typedef typename Arr::Point_2                 Point;
  typedef typename Arr::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arr::Face_const_handle       Face_const_handle;
  typedef typename Arr::Ccb_halfedge_const_circulator
                                                Ccb_halfedge_const_circulator;

  /// Construct the viewer.
  /// @param arr the arrangement to view
  /// @param title the title of the window
  Arr_2_viewer_qt(QWidget* parent, const Arr& arr,
                  const char* title = "2D Arrangement Basic Viewer") :
    Base(parent, arr, title)
  {}

  //!
  virtual void add_ccb(Ccb_halfedge_const_circulator circ) {
    const auto* traits = this->m_arr.geometry_traits();
    auto approx = traits->approximate_2_object();
    auto cmp_xy = traits->compare_xy_2_object();
    auto cmp_y = traits->compare_y_at_x_right_2_object();

    // Find the first halfedge directed from left to right
    auto curr = circ;
    do if (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT) break;
    while (++curr != circ);
    Halfedge_const_handle ext = curr;

    // Find the halfedge incident to the lexicographically smallest vertex,
    //  such that there is no other halfedge underneath.
    do {
      // Discard edges not directed from left to right:
      if (curr->direction() != CGAL::ARR_LEFT_TO_RIGHT) continue;


      auto res = cmp_xy(curr->source()->point(), ext->source()->point());

      // Discard the edges inciden to a point strictly larger than the point
      // incident to the stored extreme halfedge:
      if (res == LARGER) continue;

      // Store the edge inciden to a point strictly smaller:
      if (res == SMALLER) {
        ext = curr;
        continue;
      }

      // The incident points are equal; compare the halfedges themselves:
      if (cmp_y(curr->curve(), ext->curve(), curr->source()->point()) ==
          SMALLER)
        ext = curr;
    } while (++curr != circ);

    // Iterate, starting from the lexicographically smallest vertex
    curr = ext;
    do {
      while (curr->face() == curr->twin()->face())
        curr = curr->twin()->next();
      std::vector<typename Gt::Approximate_point_2> polyline;
      approx(curr->curve(), 10, std::back_inserter(polyline));
      auto it = polyline.begin();
      auto prev = it++;
      for (; it != polyline.end(); prev = it++) {
        this->add_segment(*prev, *it);
        this->add_point_in_face(*prev);
      }
      curr = curr->next();
    } while (curr != ext);
  }

  //!
  virtual void add_edge(Halfedge_const_handle e) {
    const auto* traits = this->m_arr.geometry_traits();
    auto approx = traits->approximate_2_object();
    std::vector<typename Gt::Approximate_point_2> polyline;
    approx(e->curve(), 10, std::back_inserter(polyline));

    auto it = polyline.begin();
    auto prev = it++;
    for (; it != polyline.end(); prev = it++) this->add_segment(*prev, *it);
  }
};

//!
template <typename GeometryTraits_2, typename Dcel>
void draw(const CGAL::Arrangement_2<GeometryTraits_2, Dcel>& arr,
          const char* title = "2D Arrangement Basic Viewer") {
  typedef GeometryTraits_2              Gt;
  typedef CGAL::Arrangement_2<Gt, Dcel> Arr;

  CGAL::Qt::init_ogl_context(4,3);

  int argc = 1;
  const char* argv[2] = {"t2_viewer", nullptr};
  QApplication app(argc, const_cast<char**>(argv));
  Arr_2_viewer_qt<Gt, Dcel> mainwindow(app.activeWindow(), arr, title);
  mainwindow.add_elements();
  mainwindow.show();
  app.exec();
}

}

#endif
