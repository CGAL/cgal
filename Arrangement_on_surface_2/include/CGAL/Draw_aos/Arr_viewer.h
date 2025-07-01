
#ifndef ARR_VIEWER_H
#define ARR_VIEWER_H
#include <CGAL/Draw_aos/helpers.h>
#include <array>
#include <boost/iterator/function_output_iterator.hpp>
#include "CGAL/Arr_trapezoid_ric_point_location.h"
#include "CGAL/Arrangement_on_surface_2.h"
#include "CGAL/Bbox_2.h"
#include "CGAL/Draw_aos/Arr_bounded_renderer.h"
#include "CGAL/Draw_aos/Arr_render_context.h"
#include "CGAL/Graphics_scene.h"
#include "CGAL/Graphics_scene_options.h"
#include "CGAL/Qt/camera.h"
#include "CGAL/unordered_flat_map.h"
#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtOpenGLWidgets/QtOpenGLWidgets>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QMouseEvent>
#include <QtGui/QKeyEvent>
#include <CGAL/Basic_viewer.h>
#include <cstddef>
#include <boost/range/iterator_range.hpp>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

namespace CGAL {

class Arr_viewer : public Qt::Basic_viewer
{
  using Basic_viewer = Qt::Basic_viewer;
  using Vertex_const_handle = Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
  using Face_const_handle = Arrangement::Face_const_handle;
  using Graphics_scene_options =
      Graphics_scene_options<Arrangement, Vertex_const_handle, Halfedge_const_handle, Face_const_handle>;
  using Point_location = Arr_trapezoid_ric_point_location<Arrangement>;

private:
  // Function to check if the camera's state has changed
  bool is_camera_changed() {
    QMatrix4x4 proj_mat, mv_mat;
    this->camera_->computeProjectionMatrix();
    this->camera_->computeModelViewMatrix();
    this->camera_->getProjectionMatrix(proj_mat.data());
    this->camera_->getModelViewMatrix(mv_mat.data());
    if(proj_mat == m_last_proj_matrix && mv_mat == m_last_modelview_matrix) {
      return false;
    }
    m_last_proj_matrix = proj_mat;
    m_last_modelview_matrix = mv_mat;
    return true;
  }

  // Bbox_2 initial_bbox() const {
  //   Bbox_2 bbox;
  //   for(const auto& vh : m_arr.vertex_handles()) {
  //     bbox += vh->point().bbox();
  //   }
  //   if(bbox.x_span() == 0 || bbox.y_span() == 0) {
  //     // make a default bbox around the degenrate rect
  //     bbox = Bbox_2(bbox.xmin() - 1, bbox.ymin() - 1, bbox.xmax() + 1, bbox.ymax() + 1);
  //   }
  //   return bbox;
  // }

  /**
   * @brief Computes the bounding box of the view from orthogonal camera.
   *
   * Currently we reverse the entire projection to find out the view bounding box.
   * @return Bbox_2
   */
  Bbox_2 view_bbox_from_camera() const {
    QMatrix4x4 modelview_matrix, projection_matrix;
    camera_->getModelViewMatrix(modelview_matrix.data());
    camera_->getProjectionMatrix(projection_matrix.data());
    QMatrix4x4 inverse_mvp = (projection_matrix * modelview_matrix).inverted();

    // Define 4 corners of the near plane in NDC (-1 to 1 in x and y)
    QVector4D clip_space_corners[] = {QVector4D(-1.0, -1.0, 0.0, 1.0), QVector4D(-1.0, 1.0, 0.0, 1.0),
                                      QVector4D(1.0, -1.0, 0.0, 1.0), QVector4D(1.0, 1.0, 0.0, 1.0)};

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();

    for(const QVector4D& corner : clip_space_corners) {
      QVector4D world = inverse_mvp * corner;
      if(world.w() != 0.0)
        world /= world.w();
      double x = world.x();
      double y = world.y();

      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }

    return Bbox_2(xmin, ymin, xmax, ymax);
  }

  double get_approx_error(const Bbox_2& bbox) const {
    std::array<GLint, 4> viewport;
    camera_->getViewport(viewport.data());
    double width = static_cast<double>(viewport[2]);
    // return bbox.x_span() / std::min(600.0, width);
    // We are testing linear traits, lets set it to inf
    return 10000;
  }

public:
  Arr_viewer(QWidget* parent,
             const Arrangement& arr,
             Graphics_scene_options options,
             const char* title = "Arrangement Viewer")
      : Basic_viewer(parent, m_scene, title)
      , m_scene_options(options)
      , m_arr(arr)
      , m_pl(arr) {}

  void render_arr(const Arrangement& arr, const Point_location& pl, const Bbox_2& bbox) {
    Arr_render_context ctx(arr, pl, get_approx_error(bbox));
    Arr_bounded_renderer renderer(ctx, bbox);
    const auto& cache = renderer.render();

    // add faces
    for(const auto& [fh, face_tris] : cache.face_cache()) {
      const auto& points = face_tris.points;
      const auto& tris = face_tris.triangles;
      bool draw_face = m_scene_options.colored_face(arr, fh);
      for(const auto& t : tris) {
        if(draw_face) {
          m_scene.face_begin(m_scene_options.face_color(arr, fh));
        } else {
          m_scene.face_begin();
        }
        for(const auto idx : t) {
          m_scene.add_point_in_face(points[idx]);
        }
        m_scene.face_end();
      }
    }

    // add edges
    for(const auto& [he, polyline] : cache.halfedge_cache()) {
      if(polyline.size() < 2) {
        continue; // skip degenerate edges
      }

      bool draw_colored_edge = m_scene_options.colored_edge(arr, he);
      auto color = draw_colored_edge ? m_scene_options.edge_color(arr, he) : CGAL::IO::Color();
      for(size_t i = 0; i < polyline.size() - 1; ++i) {
        const auto& cur_pt = polyline[i];
        const auto& next_pt = polyline[i + 1];
        auto mid_pt = CGAL::midpoint(cur_pt, next_pt);
        if(mid_pt.x() <= bbox.xmin() || mid_pt.x() > bbox.xmax() || mid_pt.y() <= bbox.ymin() ||
           mid_pt.y() > bbox.ymax())
        {
          continue;
        }
        if(draw_colored_edge) {
          m_scene.add_segment(cur_pt, next_pt, color);
        } else {
          m_scene.add_segment(cur_pt, next_pt);
        }
      }
    }
    // add vertices
    for(const auto& [vh, pt] : cache.vertex_cache()) {
      if(m_scene_options.colored_vertex(arr, vh)) {
        m_scene.add_point(pt, m_scene_options.vertex_color(arr, vh));
      } else {
        m_scene.add_point(pt);
      }
    }
  }

  void rerender(Bbox_2 bbox) {
    m_scene.clear();
    render_arr(m_arr, m_pl, bbox);
    Basic_viewer::redraw();
  }

  virtual void draw() override {
    if(is_camera_changed()) {
      Bbox_2 bbox = view_bbox_from_camera();
      // shrink the bbox by 10% for testing
      double dx = (bbox.xmax() - bbox.xmin()) * 0.1;
      double dy = (bbox.ymax() - bbox.ymin()) * 0.1;
      bbox = Bbox_2(bbox.xmin() + dx, bbox.ymin() + dy, bbox.xmax() - dx, bbox.ymax() - dy);
      std::cout << "Camera changed, recomputing arrangement bounding box: " << bbox << std::endl;
      rerender(bbox);
    }
    Basic_viewer::draw();
  }

  virtual ~Arr_viewer() {}

private:
  Graphics_scene m_scene;
  Graphics_scene_options m_scene_options;
  Arrangement m_arr;
  Point_location m_pl;
  QMatrix4x4 m_last_proj_matrix;
  QMatrix4x4 m_last_modelview_matrix;
};

void draw_viewer(const Arrangement& arr) {
  Qt::init_ogl_context(4, 3);
  int argc;
  QApplication app(argc, nullptr);
  Graphics_scene_options<Arrangement, Arrangement::Vertex_const_handle, Arrangement::Halfedge_const_handle,
                         Arrangement::Face_const_handle>
      gso;
  gso.enable_faces();
  gso.enable_edges();
  gso.enable_vertices();
  gso.face_color = [](const Arrangement&, const Arrangement::Face_const_handle& fh) -> CGAL::IO::Color {
    CGAL::Random random((size_t(&*fh)));
    return get_random_color(random);
  };
  gso.colored_face = [](const Arrangement&, const Arrangement::Face_const_handle&) { return true; };
  gso.vertex_color = [](const Arrangement&, const Arrangement::Vertex_const_handle& vh) -> CGAL::IO::Color {
    CGAL::Random random((size_t(&*vh)));
    return get_random_color(random);
  };
  // Arr_viewer viewer(app.activeWindow(), move_degenerate_features(arr), gso, "Arrangement Viewer");
  Arr_viewer viewer(app.activeWindow(), arr, gso, "Arrangement Viewer");
  std::cout << "Preprocess complete" << std::endl;
  viewer.show();
  app.exec();
}

} // namespace CGAL
#endif // ARR_VIEWER_H