// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_QT_BASIC_VIEWER_H
#define CGAL_QT_BASIC_VIEWER_H
#ifndef CGAL_USE_BASIC_VIEWER_QT
#define CGAL_USE_BASIC_VIEWER_QT 1
#endif // CGAL_USE_BASIC_VIEWER_QT
#if defined(CGAL_USE_BASIC_VIEWER_QT)

// TODO #include <CGAL/license/GraphicsView.h>
#include <CGAL/Graphics_scene.h>
#include <cstring>
#include <iostream>
#include <string>
#include <tuple>

#ifdef __GNUC__
#if __GNUC__ >= 9
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif
#endif

#include <QApplication>
#include <QKeyEvent>

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>

#ifdef __GNUC__
#if __GNUC__ >= 9
#pragma GCC diagnostic pop
#endif
#endif

#include <cfloat>
#include <cstdlib>
#include <vector>

#include <CGAL/Basic_shaders.h>
#include <CGAL/Qt/CreateOpenGLContext.h>
#include <CGAL/Qt/constraint.h>
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/assertions.h>

#define CGAL_BASIC_VIEWER_INIT_SIZE_X 500
#define CGAL_BASIC_VIEWER_INIT_SIZE_Y 450

namespace CGAL {
namespace Qt {
//------------------------------------------------------------------------------
class Basic_viewer : public CGAL::QGLViewer
{
public:
  using BufferType = float;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
  typedef Local_kernel::Point_3 Local_point;
  typedef Local_kernel::Vector_3 Local_vector;
  using GS = Graphics_scene;

  // Constructor/Destructor
  Basic_viewer(QWidget* parent,
               const Graphics_scene& buf,
               const char* title = "",
               bool draw_vertices = false,
               bool draw_edges = true,
               bool draw_faces = true,
               bool use_default_color = false,
               bool inverse_normal = false,
               bool draw_rays = true,
               bool draw_lines = true,
               bool draw_text = true,
               bool no_2D_mode = false)
      : CGAL::QGLViewer(parent)
      , m_scene(buf)
      , m_draw_vertices(draw_vertices)
      , m_draw_edges(draw_edges)
      , m_draw_rays(draw_rays)
      , m_draw_lines(draw_lines)
      , m_draw_faces(draw_faces)
      , m_draw_text(draw_text)
      , m_draw_normals(false)
      , m_draw_cylinder_edge(false)
      , m_draw_sphere_vertex(false)
      , m_draw_mesh_triangles(false)
      , m_flat_shading(true)
      , m_use_default_color(use_default_color)
      , m_use_default_color_normal(false)
      , m_display_face_normal(false)
      , m_inverse_normal(inverse_normal)
      , m_no_2D_mode(no_2D_mode)
      , m_geometry_feature_enabled(true)
      , m_prev_scene_empty(true)
      , m_default_color_normal(220, 60, 20)
      , m_ambient_color(0.6f, 0.5f, 0.5f, 0.5f)
      , m_are_buffers_initialized(false) {
    // Define 'Control+Q' as the new exit shortcut (default was 'Escape')
    setShortcut(qglviewer::EXIT_VIEWER, ::Qt::CTRL, ::Qt::Key_Q);

    // Add custom key description (see keyPressEvent).
    setKeyDescription(::Qt::Key_C, "Switch clipping plane display mode");
    setKeyDescription(::Qt::AltModifier, ::Qt::Key_C, "Toggle clipping plane rendering on/off");
    setKeyDescription(::Qt::Key_E, "Toggles edges display");
    setKeyDescription(::Qt::Key_M, "Toggles mono color");
    setKeyDescription(::Qt::Key_N, "Inverse direction of normals");
    setKeyDescription(::Qt::Key_S, "Switch between flat/Gouraud shading display");
    setKeyDescription(::Qt::Key_T, "Toggles text display");
    setKeyDescription(::Qt::Key_U, "Move camera direction upside down");
    setKeyDescription(::Qt::Key_V, "Toggles vertices display");
    setKeyDescription(::Qt::Key_W, "Toggles faces display");
    setKeyDescription(::Qt::Key_Plus, "Increase size of edges");
    setKeyDescription(::Qt::Key_Minus, "Decrease size of edges");
    setKeyDescription(::Qt::ControlModifier, ::Qt::Key_Plus, "Increase size of vertices");
    setKeyDescription(::Qt::ControlModifier, ::Qt::Key_Minus, "Decrease size of vertices");
    setKeyDescription(::Qt::Key_PageDown, "Increase light (all colors, use shift/alt/ctrl for one rgb component)");
    setKeyDescription(::Qt::Key_PageUp, "Decrease light (all colors, use shift/alt/ctrl for one rgb component)");
    setKeyDescription(::Qt::Key_O, "Toggles 2D mode only");
    setKeyDescription(::Qt::ControlModifier, ::Qt::Key_V, "Toggle vertices display as sphere");
    setKeyDescription(::Qt::ControlModifier, ::Qt::Key_E, "Toggle edges display as cylinder");
    setKeyDescription(::Qt::ControlModifier, ::Qt::Key_N, "Toggle normals display");
    setKeyDescription(::Qt::ShiftModifier, ::Qt::Key_N, "Toggle face/vertex normals for normal display");
    setKeyDescription(::Qt::ControlModifier, ::Qt::Key_M, "Toggle normals mono color");
    setKeyDescription(::Qt::ControlModifier, ::Qt::Key_T, "Toggle triangles display");
    setKeyDescription(::Qt::Key_F2, "Take a screenshot");

    // Add custom mouse description
    setMouseBindingDescription(::Qt::Key_C, ::Qt::ControlModifier, ::Qt::LeftButton,
                               "Rotate the clipping plane when enabled");
    setMouseBindingDescription(::Qt::Key_C, ::Qt::ControlModifier, ::Qt::RightButton,
                               "Translate the clipping plane when enabled");
    setMouseBindingDescription(::Qt::Key_C, ::Qt::ControlModifier, ::Qt::MiddleButton,
                               "Control the clipping plane transparency when enabled");

    setMouseBinding(::Qt::ControlModifier, ::Qt::LeftButton, qglviewer::FRAME, qglviewer::NO_MOUSE_ACTION);
    setMouseBinding(::Qt::ControlModifier, ::Qt::RightButton, qglviewer::FRAME, qglviewer::NO_MOUSE_ACTION);
    setMouseBinding(::Qt::ControlModifier, ::Qt::MiddleButton, qglviewer::FRAME, qglviewer::NO_MOUSE_ACTION);
    setWheelBinding(::Qt::ControlModifier, qglviewer::FRAME, qglviewer::NO_MOUSE_ACTION);

    setMouseBinding(::Qt::Key_C, ::Qt::ControlModifier, ::Qt::LeftButton, qglviewer::FRAME, qglviewer::ROTATE);
    setMouseBinding(::Qt::Key_C, ::Qt::ControlModifier, ::Qt::RightButton, qglviewer::FRAME, qglviewer::TRANSLATE);
    setMouseBinding(::Qt::Key_C, ::Qt::ControlModifier, ::Qt::MiddleButton, qglviewer::FRAME, qglviewer::ZOOM);
    setWheelBinding(::Qt::Key_C, ::Qt::ControlModifier, qglviewer::FRAME, qglviewer::ZOOM);

    if(title[0] == 0)
      setWindowTitle("CGAL Basic Viewer");
    else
      setWindowTitle(title);

    resize(CGAL_BASIC_VIEWER_INIT_SIZE_X, CGAL_BASIC_VIEWER_INIT_SIZE_Y);
  }

  ~Basic_viewer() {
    makeCurrent();
    for(unsigned int i = 0; i < NB_GL_BUFFERS; ++i)
      buffers[i].destroy();

    for(unsigned int i = 0; i < NB_VAO_BUFFERS; ++i)
      vao[i].destroy();

    Q_FOREACH(QOpenGLShader* shader, rendering_program_p_l.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_line.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_face.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_clipping_plane.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_sphere.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_cylinder.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_normal.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_triangle.shaders())
      delete shader;
    delete m_frame_plane;
  }

  /*****SETTERS*****/

  inline void size_vertices(float s) { m_size_vertices = s; }
  inline void size_edges(float s) { m_size_edges = s; }
  inline void size_rays(float s) { m_size_rays = s; }
  inline void size_lines(float s) { m_size_lines = s; }

  inline void draw_vertices(bool b) { m_draw_vertices = b; }
  inline void draw_edges(bool b) { m_draw_edges = b; }
  inline void draw_rays(bool b) { m_draw_rays = b; }
  inline void draw_lines(bool b) { m_draw_lines = b; }
  inline void draw_faces(bool b) { m_draw_faces = b; }
  inline void draw_text(bool b) { m_draw_text = b; }
  inline void draw_mesh_triangles(bool b) { m_draw_mesh_triangles = b; }

  inline void use_default_color(bool b) { m_use_default_color = b; }
  inline void use_default_color_normal(bool b) { m_use_default_color_normal = b; }
  inline void flat_shading(bool b) { m_flat_shading = b; }
  inline void reverse_normal(bool b) { m_inverse_normal = b; }
  inline void default_color_normals(const CGAL::IO::Color& c) { m_default_color_normal = c; }
  inline void normal_height_factor(float h) { m_height_factor_normals = h; }

  inline void toggle_draw_vertices() { m_draw_vertices = !m_draw_vertices; }
  inline void toggle_draw_edges() { m_draw_edges = !m_draw_edges; }
  inline void toggle_draw_rays() { m_draw_rays = !m_draw_rays; }
  inline void toggle_draw_lines() { m_draw_lines = !m_draw_lines; }
  inline void toggle_draw_faces() { m_draw_faces = !m_draw_faces; }
  inline void toggle_draw_text() { m_draw_text = !m_draw_text; }
  inline void toggle_use_default_color() { m_use_default_color = !m_use_default_color; }
  inline void toggle_use_normal_default_color() { m_use_default_color_normal = !m_use_default_color_normal; }
  inline void toggle_flat_shading() { m_flat_shading = !m_flat_shading; }
  inline void toggle_reverse_normal() { m_inverse_normal = !m_inverse_normal; }

  /*****************/
  /*****GETTERS*****/
  inline bool draw_vertices() const { return m_draw_vertices; }
  inline bool draw_edges() const { return m_draw_edges; }
  inline bool draw_rays() const { return m_draw_rays; }
  inline bool draw_lines() const { return m_draw_lines; }
  inline bool draw_faces() const { return m_draw_faces; }
  inline bool draw_text() const { return m_draw_text; }
  inline bool use_default_color() const { return m_use_default_color; }
  inline bool reverse_normal() const { return m_inverse_normal; }
  inline bool flat_shading() const { return m_flat_shading; }

  Local_kernel::Plane_3 clipping_plane() const {
    CGAL::qglviewer::Vec n = m_frame_plane->inverseTransformOf(CGAL::qglviewer::Vec(0.f, 0.f, 1.f));
    const CGAL::qglviewer::Vec& pos = m_frame_plane->position();
    return Local_kernel::Plane_3(n[0], n[1], n[2], -n * pos);
  }

  bool clipping_plane_enabled() const { return (m_use_clipping_plane != CLIPPING_PLANE_OFF); }

  const Graphics_scene& graphics_scene() const { return m_scene; }

  /*****************/

  virtual void redraw() {
    initialize_buffers();
    update();
    if(m_prev_scene_empty) {
      initialize_vertices_and_edges_size();
    }
    m_prev_scene_empty = (m_scene.empty());
  }

  void reverse_all_normals() {
    m_inverse_normal = !m_inverse_normal;
    m_scene.reverse_all_normals();
  }

  // Returns true if the data structure lies on a plane
  bool is_two_dimensional() const { return !m_no_2D_mode && m_scene.is_two_dimensional(); }

  virtual void draw() {
    glEnable(GL_DEPTH_TEST);

    QMatrix4x4 clipping_mMatrix;
    clipping_mMatrix.setToIdentity();
    for(int i = 0; i < 16; i++) {
      clipping_mMatrix.data()[i] = m_frame_plane->matrix()[i];
    }
    QVector4D clipPlane = clipping_mMatrix * QVector4D(0.0, 0.0, 1.0, 0.0);
    QVector4D plane_point = clipping_mMatrix * QVector4D(0, 0, 0, 1);
    if(!m_are_buffers_initialized) {
      initialize_buffers();
    }

    QVector3D color;
    attrib_buffers(this);

    if(m_draw_vertices) {
      if(m_draw_sphere_vertex && m_geometry_feature_enabled) {
        auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
          rendering_program_sphere.bind();
          if(m_use_default_color) {
            auto vertex_color = m_scene.get_default_color_point();
            color = QVector3D((double)vertex_color.red() / (double)255, (double)vertex_color.green() / (double)255,
                              (double)vertex_color.blue() / (double)255);
            rendering_program_sphere.setUniformValue("u_DefaultColor", color);
            rendering_program_sphere.setUniformValue("u_UseDefaultColor", static_cast<GLint>(1));
          } else {
            rendering_program_sphere.setUniformValue("u_UseDefaultColor", static_cast<GLint>(0));
          }
          rendering_program_sphere.setUniformValue("u_Radius",
                                                   static_cast<GLfloat>(sceneRadius() * m_size_vertices * 0.002));
          rendering_program_sphere.setUniformValue("u_ClipPlane", clipPlane);
          rendering_program_sphere.setUniformValue("u_PointPlane", plane_point);
          rendering_program_sphere.setUniformValue("u_RenderingMode", rendering_mode);

          vao[VAO_POINTS].bind();
          glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_POINTS)));
        };

        enum {
          DRAW_ALL = -1,    // draw all
          DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
          DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
        };

        if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) {
          renderer(DRAW_INSIDE_ONLY);
        } else {
          renderer(DRAW_ALL);
        }

        rendering_program_sphere.release();
      } else {
        auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
          rendering_program_p_l.bind();

          if(m_use_default_color) {
            auto vertex_color = m_scene.get_default_color_point();
            color = QVector3D((double)vertex_color.red() / (double)255, (double)vertex_color.green() / (double)255,
                              (double)vertex_color.blue() / (double)255);
            rendering_program_p_l.setUniformValue("u_DefaultColor", color);
            rendering_program_p_l.setUniformValue("u_UseDefaultColor", static_cast<GLint>(1));
          } else {
            rendering_program_p_l.setUniformValue("u_UseDefaultColor", static_cast<GLint>(0));
          }
          rendering_program_p_l.setUniformValue("u_PointSize", GLfloat(m_size_vertices));
          rendering_program_p_l.setUniformValue("u_IsOrthographic", GLint(is_two_dimensional()));

          rendering_program_p_l.setUniformValue("u_ClipPlane", clipPlane);
          rendering_program_p_l.setUniformValue("u_PointPlane", plane_point);
          rendering_program_p_l.setUniformValue("u_RenderingMode", rendering_mode);

          vao[VAO_POINTS].bind();
          glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_POINTS)));
        };

        enum {
          DRAW_ALL = -1,    // draw all
          DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
          DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
        };

        if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) {
          renderer(DRAW_INSIDE_ONLY);
        } else {
          renderer(DRAW_ALL);
        }

        rendering_program_p_l.release();
      }
    }

    if(m_draw_edges && !m_draw_mesh_triangles) {
      if(m_draw_cylinder_edge && m_geometry_feature_enabled) {
        auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
          rendering_program_cylinder.bind();

          if(m_use_default_color) {
            auto edge_color = m_scene.get_default_color_segment();
            color = QVector3D((double)edge_color.red() / (double)255, (double)edge_color.green() / (double)255,
                              (double)edge_color.blue() / (double)255);
            rendering_program_cylinder.setUniformValue("u_DefaultColor", color);
            rendering_program_cylinder.setUniformValue("u_UseDefaultColor", static_cast<GLint>(1));
          } else {
            rendering_program_cylinder.setUniformValue("u_UseDefaultColor", static_cast<GLint>(0));
          }
          rendering_program_cylinder.setUniformValue("u_Radius",
                                                     static_cast<GLfloat>(sceneRadius() * m_size_edges * 0.001));
          rendering_program_cylinder.setUniformValue("u_ClipPlane", clipPlane);
          rendering_program_cylinder.setUniformValue("u_PointPlane", plane_point);
          rendering_program_cylinder.setUniformValue("u_RenderingMode", rendering_mode);

          vao[VAO_SEGMENTS].bind();
          glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_SEGMENTS)));
        };

        enum {
          DRAW_ALL = -1,    // draw all
          DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
          DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
        };

        if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) {
          renderer(DRAW_INSIDE_ONLY);
        } else {
          renderer(DRAW_ALL);
        }

        rendering_program_cylinder.release();
      } else {
        auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
          QVector2D viewport = {CGAL_BASIC_VIEWER_INIT_SIZE_X, CGAL_BASIC_VIEWER_INIT_SIZE_Y};

          rendering_program_line.bind();

          if(m_use_default_color) {
            auto edge_color = m_scene.get_default_color_segment();
            color = QVector3D((double)edge_color.red() / (double)255, (double)edge_color.green() / (double)255,
                              (double)edge_color.blue() / (double)255);
            rendering_program_line.setUniformValue("u_DefaultColor", color);
            rendering_program_line.setUniformValue("u_UseDefaultColor", static_cast<GLint>(1));
          } else {
            rendering_program_line.setUniformValue("u_UseDefaultColor", static_cast<GLint>(0));
          }
          rendering_program_line.setUniformValue("u_PointSize", static_cast<GLfloat>(m_size_edges));
          rendering_program_line.setUniformValue("u_IsOrthographic", static_cast<GLint>(is_two_dimensional()));

          rendering_program_line.setUniformValue("u_Viewport", viewport);
          rendering_program_line.setUniformValue("u_ClipPlane", clipPlane);
          rendering_program_line.setUniformValue("u_PointPlane", plane_point);
          rendering_program_line.setUniformValue("u_RenderingMode", rendering_mode);

          vao[VAO_SEGMENTS].bind();
          glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_SEGMENTS)));
        };

        enum {
          DRAW_ALL = -1,    // draw all
          DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
          DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
        };

        if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) {
          renderer(DRAW_INSIDE_ONLY);
        } else {
          renderer(DRAW_ALL);
        }

        rendering_program_line.release();
      }
    }

    if(m_draw_rays) {
      auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
        rendering_program_p_l.bind();

        if(m_use_default_color) {
          auto ray_color = m_scene.get_default_color_ray();
          color = QVector3D((double)ray_color.red() / (double)255, (double)ray_color.green() / (double)255,
                            (double)ray_color.blue() / (double)255);
          rendering_program_p_l.setUniformValue("u_DefaultColor", color);
          rendering_program_p_l.setUniformValue("u_UseDefaultColor", static_cast<GLint>(1));
        } else {
          rendering_program_p_l.setUniformValue("u_UseDefaultColor", static_cast<GLint>(0));
        }
        rendering_program_p_l.setUniformValue("u_PointSize", GLfloat(m_size_rays));
        rendering_program_p_l.setUniformValue("u_IsOrthographic", GLint(is_two_dimensional()));

        rendering_program_p_l.setUniformValue("u_ClipPlane", clipPlane);
        rendering_program_p_l.setUniformValue("u_PointPlane", plane_point);
        rendering_program_p_l.setUniformValue("u_RenderingMode", rendering_mode);

        vao[VAO_RAYS].bind();
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_RAYS)));
      };

      enum {
        DRAW_ALL = -1,    // draw all
        DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
        DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
      };

      if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) {
        renderer(DRAW_INSIDE_ONLY);
      } else {
        renderer(DRAW_ALL);
      }

      rendering_program_p_l.release();
    }

    if(m_draw_lines) {
      auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
        rendering_program_p_l.bind();

        if(m_use_default_color) {
          auto line_color = m_scene.get_default_color_line();
          color = QVector3D((double)line_color.red() / (double)255, (double)line_color.green() / (double)255,
                            (double)line_color.blue() / (double)255);
          rendering_program_p_l.setUniformValue("u_DefaultColor", color);
          rendering_program_p_l.setUniformValue("u_UseDefaultColor", static_cast<GLint>(1));
        } else {
          rendering_program_p_l.setUniformValue("u_UseDefaultColor", static_cast<GLint>(0));
        }
        rendering_program_p_l.setUniformValue("u_PointSize", GLfloat(m_size_lines));
        rendering_program_p_l.setUniformValue("u_IsOrthographic", GLint(is_two_dimensional()));

        rendering_program_p_l.setUniformValue("u_ClipPlane", clipPlane);
        rendering_program_p_l.setUniformValue("u_PointPlane", plane_point);
        rendering_program_p_l.setUniformValue("u_RenderingMode", rendering_mode);

        vao[VAO_LINES].bind();
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_LINES)));
      };

      enum {
        DRAW_ALL = -1,    // draw all
        DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
        DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
      };

      if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) {
        renderer(DRAW_INSIDE_ONLY);
      } else {
        renderer(DRAW_ALL);
      }

      rendering_program_p_l.release();
    }

    // Fix Z-fighting by drawing faces at a depth
    GLfloat offset_factor;
    GLfloat offset_units;
    if(is_two_dimensional()) {
      glGetFloatv(GL_POLYGON_OFFSET_FACTOR, &offset_factor);
      glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
      glPolygonOffset(0.1f, 0.9f);
    }

    if(m_draw_faces) {

      // reference: https://stackoverflow.com/questions/37780345/opengl-how-to-create-order-independent-transparency
      // rendering_mode == -1: draw all as solid;
      // rendering_mode == 0: draw solid only;
      // rendering_mode == 1: draw transparent only;
      auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(2.0, 2.0);
        glDepthFunc(GL_LESS);

        rendering_program_face.bind();
        if(m_use_default_color) {
          auto face_color = m_scene.get_default_color_face();
          color = QVector3D((double)face_color.red() / (double)255, (double)face_color.green() / (double)255,
                            (double)face_color.blue() / (double)255);
          rendering_program_face.setUniformValue("u_DefaultColor", color);
          rendering_program_face.setUniformValue("u_UseDefaultColor", static_cast<GLint>(1));
        } else {
          rendering_program_face.setUniformValue("u_UseDefaultColor", static_cast<GLint>(0));
        }
        rendering_program_face.setUniformValue("u_RenderingMode", rendering_mode);
        rendering_program_face.setUniformValue("u_RenderingTransparency", clipping_plane_rendering_transparency);
        rendering_program_face.setUniformValue("u_ClipPlane", clipPlane);
        rendering_program_face.setUniformValue("u_PointPlane", plane_point);

        vao[VAO_FACES].bind();
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_FACES)));
        glDisable(GL_POLYGON_OFFSET_FILL);
      };

      auto renderer_clipping_plane = [this](bool clipping_plane_rendering) {
        if(!isOpenGL_4_3())
          return;
        if(!clipping_plane_rendering)
          return;
        // render clipping plane here
        rendering_program_clipping_plane.bind();
        vao[VAO_CLIPPING_PLANE].bind();
        glLineWidth(0.1f);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>((m_array_for_clipping_plane.size() / 3)));
        glLineWidth(1.0f);
        vao[VAO_CLIPPING_PLANE].release();
        rendering_program_clipping_plane.release();
      };

      enum {
        DRAW_SOLID_ALL = -1,  // draw all mesh in solid mode
        DRAW_SOLID_HALF,      // draw only the mesh inside the clipping plane as solid
        DRAW_TRANSPARENT_HALF // draw only the mesh outside the clipping plane as transparent
      };

      if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF) {
        // The z-buffer will prevent transparent objects from being displayed behind other transparent objects.
        // Before rendering all transparent objects, disable z-testing first.

        // 1. draw solid first
        renderer(DRAW_SOLID_HALF);

        // 2. draw transparent layer second with back face culling to avoid messy triangles
        glDepthMask(false); // disable z-testing
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CW);
        renderer(DRAW_TRANSPARENT_HALF);

        // 3. draw solid again without culling and blend to make sure the solid mesh is visible
        glDepthMask(true); // enable z-testing
        glDisable(GL_CULL_FACE);
        glDisable(GL_BLEND);
        renderer(DRAW_SOLID_HALF);

        // 4. render clipping plane here
        renderer_clipping_plane(clipping_plane_rendering);
      } else if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_WIRE_HALF ||
                m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY)
      {
        // 1. draw solid HALF
        renderer(DRAW_SOLID_HALF);

        // 2. render clipping plane here
        renderer_clipping_plane(clipping_plane_rendering);
      } else {
        // 1. draw solid FOR ALL
        renderer(DRAW_SOLID_ALL);
      }

      if(is_two_dimensional())
        glPolygonOffset(offset_factor, offset_units);

      rendering_program_face.release();
    }

    if(m_draw_normals) {
      auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
        rendering_program_normal.bind();
        if(m_use_default_color_normal) {
          color = QVector3D((double)m_default_color_normal.red() / (double)255,
                            (double)m_default_color_normal.green() / (double)255,
                            (double)m_default_color_normal.blue() / (double)255);
          rendering_program_normal.setUniformValue("u_DefaultColor", color);
          rendering_program_normal.setUniformValue("u_UseDefaultColor", static_cast<GLint>(1));
        } else {
          rendering_program_normal.setUniformValue("u_UseDefaultColor", static_cast<GLint>(0));
        }
        rendering_program_normal.setUniformValue("u_Factor", static_cast<GLfloat>(m_height_factor_normals));
        rendering_program_normal.setUniformValue("u_SceneRadius", static_cast<GLfloat>(sceneRadius()));
        if(m_display_face_normal) {
          rendering_program_normal.setUniformValue("u_DisplayFaceNormal", static_cast<GLint>(1));
        } else {
          rendering_program_normal.setUniformValue("u_DisplayFaceNormal", static_cast<GLint>(0));
        }

        rendering_program_normal.setUniformValue("u_ClipPlane", clipPlane);
        rendering_program_normal.setUniformValue("u_PointPlane", plane_point);
        rendering_program_normal.setUniformValue("u_RenderingMode", rendering_mode);

        vao[VAO_FACES].bind();
        glLineWidth(m_size_normals);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_FACES)));
      };

      enum {
        DRAW_ALL = -1,    // draw all
        DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
        DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
      };

      if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) {
        renderer(DRAW_INSIDE_ONLY);
      } else {
        renderer(DRAW_ALL);
      }

      rendering_program_normal.release();
    }

    if(m_draw_mesh_triangles) {

      auto renderer = [this, &clipPlane, &plane_point](float rendering_mode) {
        rendering_program_triangle.bind();
        rendering_program_triangle.setUniformValue("u_RenderingMode", rendering_mode);
        rendering_program_triangle.setUniformValue("u_ClipPlane", clipPlane);
        rendering_program_triangle.setUniformValue("u_PointPlane", plane_point);

        vao[VAO_FACES].bind();
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(m_scene.number_of_elements(GS::POS_FACES)));
      };

      enum {
        DRAW_ALL = -1,    // draw all
        DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
        DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
      };

      if(m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY) {
        renderer(DRAW_INSIDE_ONLY);
      } else {
        renderer(DRAW_ALL);
      }

      rendering_program_triangle.release();
    }

    if(m_draw_text) {
      glDisable(GL_LIGHTING);
      for(std::size_t i = 0; i < m_scene.m_texts_size(); ++i) {
        auto& m_texts_vec = m_scene.get_m_texts();
        CGAL::qglviewer::Vec screenPos = camera()->projectedCoordinatesOf(CGAL::qglviewer::Vec(
            std::get<0>(m_texts_vec[i]).x(), std::get<0>(m_texts_vec[i]).y(), std::get<0>(m_texts_vec[i]).z()));

        drawText((int)screenPos[0], (int)screenPos[1], QString(std::get<1>(m_texts_vec[i]).c_str()));
      }
      glEnable(GL_LIGHTING);
    }

    // Multiply matrix to get in the frame coordinate system.
    // glMultMatrixd(manipulatedFrame()->matrix()); // Linker error
    // Scale down the drawings
    // glScalef(0.3f, 0.3f, 0.3f); // Linker error
    // Draw an axis using the QGLViewer static function
    // drawAxis();
  }

protected:
  void compile_shaders() {
    rendering_program_face.removeAllShaders();
    rendering_program_p_l.removeAllShaders();
    rendering_program_line.removeAllShaders();
    rendering_program_clipping_plane.removeAllShaders();
    rendering_program_sphere.removeAllShaders();
    rendering_program_cylinder.removeAllShaders();
    rendering_program_normal.removeAllShaders();
    rendering_program_triangle.removeAllShaders();

    // Create the buffers
    for(unsigned int i = 0; i < NB_GL_BUFFERS; ++i) {
      if(!buffers[i].isCreated() && !buffers[i].create()) {
        std::cerr << "VBO Creation number " << i << " FAILED" << std::endl;
      }
    }

    for(unsigned int i = 0; i < NB_VAO_BUFFERS; ++i) {
      if(!vao[i].isCreated() && !vao[i].create()) {
        std::cerr << "VAO Creation number " << i << " FAILED" << std::endl;
      }
    }

    // Vertices and segments shader

    // const char* source_ = isOpenGL_4_3()
    //     ? VERTEX_SOURCE_P_L
    //     : VERTEX_SOURCE_P_L_COMP;

    const char* source_ = isOpenGL_4_3() ? VERTEX_SOURCE_P_L : VERTEX_SOURCE_P_L_COMP;

    QOpenGLShader* vertex_shader_p_l = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_p_l->compileSourceCode(source_)) {
      std::cerr << "Compiling vertex source FAILED" << std::endl;
    }

    source_ = isOpenGL_4_3() ? FRAGMENT_SOURCE_P_L : FRAGMENT_SOURCE_P_L_COMP;

    QOpenGLShader* fragment_shader_p_l = new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader_p_l->compileSourceCode(source_)) {
      std::cerr << "Compiling fragment source FAILED" << std::endl;
    }

    if(!rendering_program_p_l.addShader(vertex_shader_p_l)) {
      std::cerr << "adding vertex shader FAILED" << std::endl;
    }
    if(!rendering_program_p_l.addShader(fragment_shader_p_l)) {
      std::cerr << "adding fragment shader FAILED" << std::endl;
    }
    if(!rendering_program_p_l.link()) {
      std::cerr << "linking Program FAILED" << std::endl;
    }

    // Faces shader

    source_ = isOpenGL_4_3() ? VERTEX_SOURCE_COLOR : VERTEX_SOURCE_COLOR_COMP;

    QOpenGLShader* vertex_shader_face = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_face->compileSourceCode(source_)) {
      std::cerr << "Compiling vertex source FAILED" << std::endl;
    }

    source_ = isOpenGL_4_3() ? FRAGMENT_SOURCE_COLOR : FRAGMENT_SOURCE_COLOR_COMP;

    QOpenGLShader* fragment_shader_face = new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader_face->compileSourceCode(source_)) {
      std::cerr << "Compiling fragment source FAILED" << std::endl;
    }

    if(!rendering_program_face.addShader(vertex_shader_face)) {
      std::cerr << "adding vertex shader FAILED" << std::endl;
    }
    if(!rendering_program_face.addShader(fragment_shader_face)) {
      std::cerr << "adding fragment shader FAILED" << std::endl;
    }
    if(!rendering_program_face.link()) {
      std::cerr << "linking Program FAILED" << std::endl;
    }

    if(isOpenGL_4_3()) {
      // clipping plane shader
      source_ = VERTEX_SOURCE_CLIPPING_PLANE;

      QOpenGLShader* vertex_shader_clipping_plane = new QOpenGLShader(QOpenGLShader::Vertex);
      if(!vertex_shader_clipping_plane->compileSourceCode(source_)) {
        std::cerr << "Compiling vertex source for clipping plane FAILED" << std::endl;
      }

      source_ = FRAGMENT_SOURCE_CLIPPING_PLANE;

      QOpenGLShader* fragment_shader_clipping_plane = new QOpenGLShader(QOpenGLShader::Fragment);
      if(!fragment_shader_clipping_plane->compileSourceCode(source_)) {
        std::cerr << "Compiling fragment source for clipping plane FAILED" << std::endl;
      }

      if(!rendering_program_clipping_plane.addShader(vertex_shader_clipping_plane)) {
        std::cerr << "Adding vertex shader for clipping plane FAILED" << std::endl;
      }
      if(!rendering_program_clipping_plane.addShader(fragment_shader_clipping_plane)) {
        std::cerr << "Adding fragment shader for clipping plane FAILED" << std::endl;
      }
      if(!rendering_program_clipping_plane.link()) {
        std::cerr << "Linking Program for clipping plane FAILED" << std::endl;
      }
    }

    // source_ = isOpenGL_4_3()
    //         ? VERTEX_SOURCE_CLIPPING_PLANE
    //         : vertex_source_clipping_plane_comp;

    // QOpenGLShader *vertex_shader_clipping_plane = new QOpenGLShader(QOpenGLShader::Vertex);
    // if (!vertex_shader_clipping_plane->compileSourceCode(source_))
    // { std::cerr << "Compiling vertex source for clipping plane FAILED" << std::endl; }

    // source_ = isOpenGL_4_3()
    //         ? FRAGMENT_SOURCE_CLIPPING_PLANE
    //         : fragment_source_clipping_plane_comp;

    // QOpenGLShader *fragment_shader_clipping_plane = new QOpenGLShader(QOpenGLShader::Fragment);
    // if (!fragment_shader_clipping_plane->compileSourceCode(source_))
    // { std::cerr << "Compiling fragment source for clipping plane FAILED" << std::endl; }

    // if (!rendering_program_clipping_plane.addShader(vertex_shader_clipping_plane))
    // { std::cerr << "Adding vertex shader for clipping plane FAILED" << std::endl;}
    // if (!rendering_program_clipping_plane.addShader(fragment_shader_clipping_plane))
    // { std::cerr << "Adding fragment shader for clipping plane FAILED" << std::endl; }
    // if (!rendering_program_clipping_plane.link())
    // { std::cerr << "Linking Program for clipping plane FAILED" << std::endl; }

    // Sphere shader
    if(isOpenGL_4_3()) {
      source_ = VERTEX_SOURCE_SHAPE;

      QOpenGLShader* vertex_shader_sphere = new QOpenGLShader(QOpenGLShader::Vertex);
      if(!vertex_shader_sphere->compileSourceCode(source_)) {
        std::cerr << "Compiling vertex source for sphere FAILED" << std::endl;
      }

      source_ = GEOMETRY_SOURCE_SPHERE;

      QOpenGLShader* geometry_shader_sphere = new QOpenGLShader(QOpenGLShader::Geometry);
      if(!geometry_shader_sphere->compileSourceCode(source_)) {
        std::cerr << "Compiling geometry source for sphere FAILED" << std::endl;
      }

      source_ = FRAGMENT_SOURCE_P_L;

      QOpenGLShader* fragment_shader_sphere = new QOpenGLShader(QOpenGLShader::Fragment);
      if(!fragment_shader_sphere->compileSourceCode(source_)) {
        std::cerr << "Compiling fragment source for sphere FAILED" << std::endl;
      }

      if(!rendering_program_sphere.addShader(vertex_shader_sphere)) {
        std::cerr << "Adding vertex shader for sphere FAILED" << std::endl;
      }
      if(!rendering_program_sphere.addShader(geometry_shader_sphere)) {
        std::cerr << "Adding geometry shader for sphere FAILED" << std::endl;
      }
      if(!rendering_program_sphere.addShader(fragment_shader_sphere)) {
        std::cerr << "Adding fragment shader for clipping plane FAILED" << std::endl;
      }
      if(!rendering_program_sphere.link()) {
        std::cerr << "Linking Program for sphere FAILED" << std::endl;
      }
    }

    // Cylinder shader
    if(isOpenGL_4_3()) {
      // clipping plane shader
      source_ = VERTEX_SOURCE_SHAPE;

      QOpenGLShader* vertex_shader_cylinder = new QOpenGLShader(QOpenGLShader::Vertex);
      if(!vertex_shader_cylinder->compileSourceCode(source_)) {
        std::cerr << "Compiling vertex source for cylinder FAILED" << std::endl;
      }

      source_ = GEOMETRY_SOURCE_CYLINDER;

      QOpenGLShader* geometry_shader_cylinder = new QOpenGLShader(QOpenGLShader::Geometry);
      if(!geometry_shader_cylinder->compileSourceCode(source_)) {
        std::cerr << "Compiling geometry source for cylinder FAILED" << std::endl;
      }

      source_ = FRAGMENT_SOURCE_P_L;

      QOpenGLShader* fragment_shader_cylinder = new QOpenGLShader(QOpenGLShader::Fragment);
      if(!fragment_shader_cylinder->compileSourceCode(source_)) {
        std::cerr << "Compiling fragment source for cylinder FAILED" << std::endl;
      }

      if(!rendering_program_cylinder.addShader(vertex_shader_cylinder)) {
        std::cerr << "Adding vertex shader for cylinder FAILED" << std::endl;
      }
      if(!rendering_program_cylinder.addShader(geometry_shader_cylinder)) {
        std::cerr << "Adding geometry shader for cylinder FAILED" << std::endl;
      }
      if(!rendering_program_cylinder.addShader(fragment_shader_cylinder)) {
        std::cerr << "Adding fragment shader for clipping plane FAILED" << std::endl;
      }
      if(!rendering_program_cylinder.link()) {
        std::cerr << "Linking Program for cylinder FAILED" << std::endl;
      }
    }

    // Normal shader
    if(isOpenGL_4_3()) {
      source_ = VERTEX_SOURCE_NORMAL;

      QOpenGLShader* vertex_shader_normal = new QOpenGLShader(QOpenGLShader::Vertex);
      if(!vertex_shader_normal->compileSourceCode(source_)) {
        std::cerr << "Compiling vertex source for normal FAILED" << std::endl;
      }

      source_ = GEOMETRY_SOURCE_NORMAL;

      QOpenGLShader* geometry_shader_normal = new QOpenGLShader(QOpenGLShader::Geometry);
      if(!geometry_shader_normal->compileSourceCode(source_)) {
        std::cerr << "Compiling geometry source for normal FAILED" << std::endl;
      }

      source_ = FRAGMENT_SOURCE_P_L;

      QOpenGLShader* fragment_shader_normal = new QOpenGLShader(QOpenGLShader::Fragment);
      if(!fragment_shader_normal->compileSourceCode(source_)) {
        std::cerr << "Compiling fragment source for normal FAILED" << std::endl;
      }

      if(!rendering_program_normal.addShader(vertex_shader_normal)) {
        std::cerr << "Adding vertex shader for normal FAILED" << std::endl;
      }
      if(!rendering_program_normal.addShader(geometry_shader_normal)) {
        std::cerr << "Adding geometry shader for normal FAILED" << std::endl;
      }
      if(!rendering_program_normal.addShader(fragment_shader_normal)) {
        std::cerr << "Adding fragment shader for clipping plane FAILED" << std::endl;
      }
      if(!rendering_program_normal.link()) {
        std::cerr << "Linking Program for normal FAILED" << std::endl;
      }
    }

    // Normal shader
    if(isOpenGL_4_3()) {
      source_ = VERTEX_SOURCE_TRIANGLE;

      QOpenGLShader* vertex_shader_triangle = new QOpenGLShader(QOpenGLShader::Vertex);
      if(!vertex_shader_triangle->compileSourceCode(source_)) {
        std::cerr << "Compiling vertex source for triangle FAILED" << std::endl;
      }

      source_ = GEOMETRY_SOURCE_TRIANGLE;

      QOpenGLShader* geometry_shader_triangle = new QOpenGLShader(QOpenGLShader::Geometry);
      if(!geometry_shader_triangle->compileSourceCode(source_)) {
        std::cerr << "Compiling geometry source for triangle FAILED" << std::endl;
      }

      source_ = FRAGMENT_SOURCE_P_L;

      QOpenGLShader* fragment_shader_triangle = new QOpenGLShader(QOpenGLShader::Fragment);
      if(!fragment_shader_triangle->compileSourceCode(source_)) {
        std::cerr << "Compiling fragment source for triangle FAILED" << std::endl;
      }

      if(!rendering_program_triangle.addShader(vertex_shader_triangle)) {
        std::cerr << "Adding vertex shader for triangle FAILED" << std::endl;
      }
      if(!rendering_program_triangle.addShader(geometry_shader_triangle)) {
        std::cerr << "Adding geometry shader for triangle FAILED" << std::endl;
      }
      if(!rendering_program_triangle.addShader(fragment_shader_triangle)) {
        std::cerr << "Adding fragment shader for clipping plane FAILED" << std::endl;
      }
      if(!rendering_program_triangle.link()) {
        std::cerr << "Linking Program for triangle FAILED" << std::endl;
      }
    }

    // Line shader
    if(isOpenGL_4_3()) {
      source_ = VERTEX_SOURCE_LINE_WIDTH;

      QOpenGLShader* vertex_shader_line = new QOpenGLShader(QOpenGLShader::Vertex);
      if(!vertex_shader_line->compileSourceCode(source_)) {
        std::cerr << "Compiling vertex source for line FAILED" << std::endl;
      }

      source_ = GEOMETRY_SOURCE_LINE_WIDTH;

      QOpenGLShader* geometry_shader_line = new QOpenGLShader(QOpenGLShader::Geometry);
      if(!geometry_shader_line->compileSourceCode(source_)) {
        std::cerr << "Compiling geometry source for line FAILED" << std::endl;
      }

      source_ = FRAGMENT_SOURCE_P_L;

      QOpenGLShader* fragment_shader_line = new QOpenGLShader(QOpenGLShader::Fragment);
      if(!fragment_shader_line->compileSourceCode(source_)) {
        std::cerr << "Compiling fragment source for line FAILED" << std::endl;
      }

      if(!rendering_program_line.addShader(vertex_shader_line)) {
        std::cerr << "Adding vertex shader for line FAILED" << std::endl;
      }
      if(!rendering_program_line.addShader(geometry_shader_line)) {
        std::cerr << "Adding geometry shader for line FAILED" << std::endl;
      }
      if(!rendering_program_line.addShader(fragment_shader_line)) {
        std::cerr << "Adding fragment shader for line FAILED" << std::endl;
      }
      if(!rendering_program_line.link()) {
        std::cerr << "Linking Program for line FAILED" << std::endl;
      }
    }
  }

  void initialize_buffers() {
    set_camera_mode();
    rendering_program_p_l.bind();

    unsigned int bufn = 0;
    std::vector<float> positions, normals, colors;
    // 1) POINT SHADER

    vao[VAO_POINTS].bind();
    positions = m_scene.get_array_of_index(GS::POS_POINTS);
    colors = m_scene.get_array_of_index(GS::COLOR_POINTS);

    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(positions.data(), static_cast<int>(positions.size() * sizeof(float)));
    rendering_program_p_l.enableAttributeArray("a_Pos");
    rendering_program_p_l.setAttributeBuffer("a_Pos", GL_FLOAT, 0, 3);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(colors.data(), static_cast<int>(colors.size() * sizeof(float)));
    rendering_program_p_l.enableAttributeArray("a_Color");
    rendering_program_p_l.setAttributeBuffer("a_Color", GL_FLOAT, 0, 3);

    // 2) SEGMENT SHADER

    vao[VAO_SEGMENTS].bind();
    positions = m_scene.get_array_of_index(GS::POS_SEGMENTS);
    colors = m_scene.get_array_of_index(GS::COLOR_SEGMENTS);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(positions.data(), static_cast<int>(positions.size() * sizeof(float)));
    rendering_program_p_l.enableAttributeArray("a_Pos");
    rendering_program_p_l.setAttributeBuffer("a_Pos", GL_FLOAT, 0, 3);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(colors.data(), static_cast<int>(colors.size() * sizeof(float)));
    rendering_program_p_l.enableAttributeArray("a_Color");
    rendering_program_p_l.setAttributeBuffer("a_Color", GL_FLOAT, 0, 3);

    // 3) RAYS SHADER

    vao[VAO_RAYS].bind();
    positions = m_scene.get_array_of_index(GS::POS_RAYS);
    colors = m_scene.get_array_of_index(GS::COLOR_RAYS);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(positions.data(), static_cast<int>(positions.size() * sizeof(float)));
    rendering_program_p_l.enableAttributeArray("a_Pos");
    rendering_program_p_l.setAttributeBuffer("a_Pos", GL_FLOAT, 0, 3);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(colors.data(), static_cast<int>(colors.size() * sizeof(float)));
    rendering_program_p_l.enableAttributeArray("a_Color");
    rendering_program_p_l.setAttributeBuffer("a_Color", GL_FLOAT, 0, 3);

    // 4) LINES SHADER

    vao[VAO_LINES].bind();
    positions = m_scene.get_array_of_index(GS::POS_LINES);
    colors = m_scene.get_array_of_index(GS::COLOR_LINES);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(positions.data(), static_cast<int>(positions.size() * sizeof(float)));
    rendering_program_p_l.enableAttributeArray("a_Pos");
    rendering_program_p_l.setAttributeBuffer("a_Pos", GL_FLOAT, 0, 3);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(colors.data(), static_cast<int>(colors.size() * sizeof(float)));
    rendering_program_p_l.enableAttributeArray("a_Color");
    rendering_program_p_l.setAttributeBuffer("a_Color", GL_FLOAT, 0, 3);

    // 5) FACE SHADER

    vao[VAO_FACES].bind();
    positions = m_scene.get_array_of_index(GS::POS_FACES);
    normals = m_scene.get_array_of_index(m_flat_shading ? GS::FLAT_NORMAL_FACES : GS::SMOOTH_NORMAL_FACES);
    colors = m_scene.get_array_of_index(GS::COLOR_FACES);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(positions.data(), static_cast<int>(positions.size() * sizeof(float)));
    rendering_program_face.enableAttributeArray("a_Pos");
    rendering_program_face.setAttributeBuffer("a_Pos", GL_FLOAT, 0, 3);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(normals.data(), static_cast<int>(normals.size() * sizeof(float)));
    rendering_program_face.enableAttributeArray("a_Normal");
    rendering_program_face.setAttributeBuffer("a_Normal", GL_FLOAT, 0, 3);

    ++bufn;
    CGAL_assertion(bufn < NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(colors.data(), static_cast<int>(colors.size() * sizeof(float)));
    rendering_program_face.enableAttributeArray("a_Color");
    rendering_program_face.setAttributeBuffer("a_Color", GL_FLOAT, 0, 3);

    // 6) clipping plane shader
    if(isOpenGL_4_3()) {
      generate_clipping_plane();

      rendering_program_clipping_plane.bind();

      vao[VAO_CLIPPING_PLANE].bind();
      ++bufn;
      CGAL_assertion(bufn < NB_GL_BUFFERS);
      buffers[bufn].bind();
      buffers[bufn].allocate(m_array_for_clipping_plane.data(),
                             static_cast<int>(m_array_for_clipping_plane.size() * sizeof(BufferType)));
      rendering_program_clipping_plane.enableAttributeArray("a_Pos");
      rendering_program_clipping_plane.setAttributeBuffer("a_Pos", GL_FLOAT, 0, 3);

      buffers[bufn].release();

      rendering_program_clipping_plane.release();
    }

    m_are_buffers_initialized = true;
  }

  void attrib_buffers(CGAL::QGLViewer* viewer) {
    QMatrix4x4 mvpMatrix;
    QMatrix4x4 mvMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(unsigned int i = 0; i < 16; i++) {
      mvpMatrix.data()[i] = (float)mat[i];
    }
    viewer->camera()->getModelViewMatrix(mat);
    for(unsigned int i = 0; i < 16; i++) {
      mvMatrix.data()[i] = (float)mat[i];
    }
    // define material
    QVector4D diffuse(0.9f, 0.9f, 0.9f, 0.9f);

    QVector4D specular(0.0f, 0.0f, 0.0f, 1.0f);

    CGAL::Bbox_3 bb;
    if(bb == m_scene.bounding_box()) // Case of "empty" bounding box
    {
      bb = Local_point(CGAL::ORIGIN).bbox();
      bb = bb + Local_point(1, 1, 1).bbox(); // To avoid a warning from Qglviewer
    } else {
      bb = m_scene.bounding_box();
    }

    QVector4D position((bb.xmax() - bb.xmin()) / 2, (bb.ymax() - bb.ymin()) / 2, bb.zmax(), 0.0);
    GLfloat shininess = 1.0f;

    rendering_program_face.bind();
    int mvpLocation = rendering_program_face.uniformLocation("u_Mvp");
    int mvLocation = rendering_program_face.uniformLocation("u_Mv");
    int lightLocation[5];
    lightLocation[0] = rendering_program_face.uniformLocation("u_LightPos");
    lightLocation[1] = rendering_program_face.uniformLocation("u_LightDiff");
    lightLocation[2] = rendering_program_face.uniformLocation("u_LightSpec");
    lightLocation[3] = rendering_program_face.uniformLocation("u_LightAmb");
    lightLocation[4] = rendering_program_face.uniformLocation("u_SpecPower");

    rendering_program_face.setUniformValue(lightLocation[0], position);
    rendering_program_face.setUniformValue(lightLocation[1], diffuse);
    rendering_program_face.setUniformValue(lightLocation[2], specular);
    rendering_program_face.setUniformValue(lightLocation[3], m_ambient_color);
    rendering_program_face.setUniformValue(lightLocation[4], shininess);
    rendering_program_face.setUniformValue(mvpLocation, mvpMatrix);
    rendering_program_face.setUniformValue(mvLocation, mvMatrix);
    rendering_program_face.release();

    rendering_program_p_l.bind();
    mvpLocation = rendering_program_p_l.uniformLocation("u_Mvp");
    rendering_program_p_l.setUniformValue(mvpLocation, mvpMatrix);
    rendering_program_p_l.release();

    // cylinder edge feature
    rendering_program_cylinder.bind();
    mvpLocation = rendering_program_cylinder.uniformLocation("u_Mvp");
    rendering_program_cylinder.setUniformValue(mvpLocation, mvpMatrix);
    rendering_program_cylinder.release();

    // sphere vertex feature
    rendering_program_sphere.bind();
    mvpLocation = rendering_program_sphere.uniformLocation("u_Mvp");
    rendering_program_sphere.setUniformValue(mvpLocation, mvpMatrix);
    rendering_program_sphere.release();

    if(isOpenGL_4_3()) {
      QMatrix4x4 clipping_mMatrix;
      clipping_mMatrix.setToIdentity();
      for(int i = 0; i < 16; i++) {
        clipping_mMatrix.data()[i] = m_frame_plane->matrix()[i];
      }

      rendering_program_clipping_plane.bind();
      int vpLocation = rendering_program_clipping_plane.uniformLocation("u_Vp");
      int mLocation = rendering_program_clipping_plane.uniformLocation("u_M");
      rendering_program_clipping_plane.setUniformValue(vpLocation, mvpMatrix);
      rendering_program_clipping_plane.setUniformValue(mLocation, clipping_mMatrix);
      rendering_program_clipping_plane.release();
    }

    if(isOpenGL_4_3()) {
      rendering_program_normal.bind();

      QMatrix4x4 projection;
      double mat[16];
      viewer->camera()->getProjectionMatrix(mat);
      for(unsigned int i = 0; i < 16; i++) {
        projection.data()[i] = (float)mat[i];
      }

      int mvLocation = rendering_program_normal.uniformLocation("u_Mv");
      int pLocation = rendering_program_normal.uniformLocation("u_Projection");
      rendering_program_normal.setUniformValue(mvLocation, mvMatrix);
      rendering_program_normal.setUniformValue(pLocation, projection);
      rendering_program_normal.release();
    }

    if(isOpenGL_4_3()) {
      rendering_program_triangle.bind();

      int mvpLocation = rendering_program_triangle.uniformLocation("u_Mvp");
      rendering_program_triangle.setUniformValue(mvpLocation, mvpMatrix);
      rendering_program_triangle.release();
    }

    if(isOpenGL_4_3()) {
      rendering_program_line.bind();

      int mvpLocation = rendering_program_line.uniformLocation("u_Mvp");
      rendering_program_line.setUniformValue(mvpLocation, mvpMatrix);
      rendering_program_line.release();
    }
  }

  void set_camera_mode() {
    if(is_two_dimensional()) {
      camera()->setType(CGAL::qglviewer::Camera::ORTHOGRAPHIC);
      //      Camera Constraint:
      constraint.setRotationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::AXIS);
      constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);

      double cx = 0., cy = 0., cz = 0.;
      if(m_scene.has_zero_x()) {
        cx = 1.;
      } else if(m_scene.has_zero_y()) {
        cy = 1.;
      } else {
        cz = 1.;
      }

      camera()->setViewDirection(CGAL::qglviewer::Vec(-cx, -cy, -cz));
      constraint.setRotationConstraintDirection(CGAL::qglviewer::Vec(cx, cy, cz));
      camera()->frame()->setConstraint(&constraint);
    } else {
      camera()->setType(CGAL::qglviewer::Camera::PERSPECTIVE);
      camera()->frame()->setConstraint(nullptr);
    }
  }

  virtual void init() {
    set_camera_mode();
    initializeOpenGLFunctions();

    // Light default parameters
    glLineWidth(m_size_edges);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.f, 1.f);
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glDisable(GL_BLEND);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_LINE_SMOOTH);
    glDisable(GL_POLYGON_SMOOTH_HINT);
    glBlendFunc(GL_ONE, GL_ZERO);
    glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);

    compile_shaders();

    CGAL::Bbox_3 bb;
    if(bb == m_scene.bounding_box()) // Case of "empty" bounding box
    {
      bb = Local_point(CGAL::ORIGIN).bbox();
      bb = bb + Local_point(1, 1, 1).bbox(); // To avoid a warning from Qglviewer
    } else {
      bb = m_scene.bounding_box();
    }
    this->camera()->setSceneBoundingBox(CGAL::qglviewer::Vec(bb.xmin(), bb.ymin(), bb.zmin()),
                                        CGAL::qglviewer::Vec(bb.xmax(), bb.ymax(), bb.zmax()));

    m_frame_plane = new CGAL::qglviewer::ManipulatedFrame;

    // Check for geometry shader availability
    int max_geometry_output_vertices = 0;
    glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES, &max_geometry_output_vertices);
    int max_geometry_output_components = 0;
    glGetIntegerv(GL_MAX_GEOMETRY_TOTAL_OUTPUT_COMPONENTS, &max_geometry_output_components);

    if(max_geometry_output_vertices < 128 || max_geometry_output_components < 1024) {
      std::cout << "Cylinder edge and sphere vertex feature disabled! (max_geometry_output_vertices="
                << max_geometry_output_vertices << ", max_geometry_output_components=" << max_geometry_output_components
                << ")\n";
      m_geometry_feature_enabled = false;
    }

    /// This code cannot be done in the constructor, because the Graphics_scene
    /// is not yet created (cf. for example LCC demo).
    if(m_inverse_normal) {
      reverse_all_normals();
    }

    initialize_vertices_and_edges_size();

    this->showEntireScene();
  }

  void initialize_vertices_and_edges_size() {
    if(!m_scene.empty()) {
      auto& bbox = m_scene.bounding_box();
      double d = CGAL::sqrt(CGAL::squared_distance(Local_point(bbox.xmin(), bbox.ymin(), bbox.zmin()),
                                                   Local_point(bbox.xmax(), bbox.ymax(), bbox.zmax())));
      // std::cout<<"Length of the diagonal: "<<d<<std::endl;
      m_size_vertices = 1.5 * d;
      m_size_edges = d;
      m_size_rays = m_size_edges;
      m_size_lines = m_size_edges;
      m_size_normals = d / 3;
      m_height_factor_normals = 0.02;
    }
  }

  void generate_clipping_plane() {
    qreal size = ((m_scene.bounding_box().xmax() - m_scene.bounding_box().xmin()) +
                  (m_scene.bounding_box().ymax() - m_scene.bounding_box().ymin()) +
                  (m_scene.bounding_box().zmax() - m_scene.bounding_box().zmin()));
    const unsigned int nb_subdivisions = 30;

    auto& array = m_array_for_clipping_plane;
    array.clear();
    for(unsigned int i = 0; i <= nb_subdivisions; ++i) {
      const float pos = float(size * (2.0 * i / nb_subdivisions - 1.0));
      array.push_back(pos);
      array.push_back(float(-size));
      array.push_back(0.f);

      array.push_back(pos);
      array.push_back(float(+size));
      array.push_back(0.f);

      array.push_back(float(-size));
      array.push_back(pos);
      array.push_back(0.f);

      array.push_back(float(size));
      array.push_back(pos);
      array.push_back(0.f);
    }

    // Normal
    array.push_back(0.f);
    array.push_back(0.f);
    array.push_back(0.f);
    array.push_back(0.f);
    array.push_back(0.f);
    array.push_back(1.f);
  }

  virtual void mouseDoubleClickEvent(QMouseEvent* e) {
    if((e->modifiers() == ::Qt::ShiftModifier) && (e->button() == ::Qt::LeftButton)) {
      if(manipulatedFrame()) {
        camera()->frame()->alignWithFrame(manipulatedFrame(), true);
      }
    } else {
      CGAL::QGLViewer::mouseDoubleClickEvent(e);
    }
  }

  virtual void keyPressEvent(QKeyEvent* e) {
    if(!on_key_pressed || !on_key_pressed(e, this)) {
      const ::Qt::KeyboardModifiers modifiers = e->modifiers();
      if((e->key() == ::Qt::Key_C) && (modifiers == ::Qt::NoButton)) {
        if(!isOpenGL_4_3())
          return;
        if(!is_two_dimensional()) {
          // toggle clipping plane
          m_use_clipping_plane = (m_use_clipping_plane + 1) % CLIPPING_PLANE_END_INDEX;
          if(m_use_clipping_plane == CLIPPING_PLANE_OFF) {
            setManipulatedFrame(nullptr);
          } else {
            setManipulatedFrame(m_frame_plane);
          }

          switch(m_use_clipping_plane) {
          case CLIPPING_PLANE_OFF:
            displayMessage(QString("Draw clipping = false"));
            break;
          case CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF:
            clipping_plane_rendering = true;
            displayMessage(QString("Draw clipping = solid half & transparent half"));
            break;
          case CLIPPING_PLANE_SOLID_HALF_WIRE_HALF:
            displayMessage(QString("Draw clipping = solid half & wireframe half"));
            break;
          case CLIPPING_PLANE_SOLID_HALF_ONLY:
            displayMessage(QString("Draw clipping = solid half only"));
            break;
          default:
            break;
          }
          update();
        }
      }

      else if((e->key() == ::Qt::Key_C) && (modifiers == ::Qt::AltModifier))
      {
        if(!isOpenGL_4_3())
          return;
        if(m_use_clipping_plane != CLIPPING_PLANE_OFF) {
          clipping_plane_rendering = !clipping_plane_rendering;
          displayMessage(QString("Draw clipping plane=%1.").arg(clipping_plane_rendering ? "true" : "false"));
          update();
        }
      } else if((e->key() == ::Qt::Key_E) && (modifiers == ::Qt::NoButton)) {
        m_draw_edges = !m_draw_edges;
        displayMessage(QString("Draw edges=%1.").arg(m_draw_edges ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_M) && (modifiers == ::Qt::NoButton)) {
        m_use_default_color = !m_use_default_color;
        displayMessage(QString("Mono color=%1.").arg(m_use_default_color ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_N) && (modifiers == ::Qt::NoButton)) {
        reverse_all_normals();
        displayMessage(QString("Inverse normal=%1.").arg(m_inverse_normal ? "true" : "false"));
        redraw();
      } else if((e->key() == ::Qt::Key_S) && (modifiers == ::Qt::NoButton)) {
        m_flat_shading = !m_flat_shading;
        if(m_flat_shading)
          displayMessage("Flat shading.");
        else
          displayMessage("Gouraud shading.");
        redraw();
      } else if((e->key() == ::Qt::Key_T) && (modifiers == ::Qt::NoButton)) {
        m_draw_text = !m_draw_text;
        displayMessage(QString("Draw text=%1.").arg(m_draw_text ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_U) && (modifiers == ::Qt::NoButton)) {
        if(is_two_dimensional()) {
          displayMessage(QString("Move camera direction upside down."));
          /* CGAL::qglviewer::Vec cur=camera()->viewDirection(); // TODO !
          double cx=cur.x, cy=cur.y, cz=cur.z;
          if (has_zero_x())      { cx=-cx; }
          else if (has_zero_y()) { cy=-cy; }
          else                   { cz=-cz; }
          double cx=0., cy=0., cz=0.;
          if (has_zero_x())      { cx=(cur.x<0?-1.:1); }
          else if (has_zero_y()) { cy=(cur.y<0?-1.:1); }
          else                   { cz=(cur.z<0?-1.:1); }*/

          camera()->setUpVector(-camera()->upVector());
          // camera()->frame()->setConstraint(NULL);
          //  camera()->setViewDirection(CGAL::qglviewer::Vec(-cx,-cy,-cz));
          // constraint.setRotationConstraintDirection(CGAL::qglviewer::Vec(cx, cy, cz));
          // camera()->frame()->setConstraint(&constraint);
          // update();
          redraw();
        }
      } else if((e->key() == ::Qt::Key_V) && (modifiers == ::Qt::NoButton)) {
        m_draw_vertices = !m_draw_vertices;
        displayMessage(QString("Draw vertices=%1.").arg(m_draw_vertices ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_W) && (modifiers == ::Qt::NoButton)) {
        m_draw_faces = !m_draw_faces;
        displayMessage(QString("Draw faces=%1.").arg(m_draw_faces ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_Plus) && (!modifiers.testFlag(::Qt::ControlModifier))) // No ctrl
      {
        m_size_edges += .5;
        displayMessage(QString("Size of edges=%1.").arg(m_size_edges));
        update();
      } else if((e->key() == ::Qt::Key_Minus) && (!modifiers.testFlag(::Qt::ControlModifier))) // No ctrl
      {
        if(m_size_edges > .5)
          m_size_edges -= .5;
        displayMessage(QString("Size of edges=%1.").arg(m_size_edges));
        update();
      } else if((e->key() == ::Qt::Key_Plus) && (modifiers.testFlag(::Qt::ControlModifier))) {
        m_size_vertices += .5;
        displayMessage(QString("Size of points=%1.").arg(m_size_vertices));
        update();
      } else if((e->key() == ::Qt::Key_Minus) && (modifiers.testFlag(::Qt::ControlModifier))) {
        if(m_size_vertices > .5)
          m_size_vertices -= .5;
        displayMessage(QString("Size of points=%1.").arg(m_size_vertices));
        update();
      } else if((e->key() == ::Qt::Key_PageUp) && (modifiers == ::Qt::NoButton)) {
        m_ambient_color.setX(m_ambient_color.x() + .1);
        if(m_ambient_color.x() > 1.)
          m_ambient_color.setX(1.);
        m_ambient_color.setY(m_ambient_color.x() + .1);
        if(m_ambient_color.y() > 1.)
          m_ambient_color.setY(1.);
        m_ambient_color.setZ(m_ambient_color.x() + .1);
        if(m_ambient_color.z() > 1.)
          m_ambient_color.setZ(1.);
        displayMessage(QString("Light color=(%1 %2 %3).")
                           .arg(m_ambient_color.x())
                           .arg(m_ambient_color.y())
                           .arg(m_ambient_color.z()));
        update();
      } else if((e->key() == ::Qt::Key_PageDown) && (modifiers == ::Qt::NoButton)) {
        m_ambient_color.setX(m_ambient_color.x() - .1);
        if(m_ambient_color.x() < 0.)
          m_ambient_color.setX(0.);
        m_ambient_color.setY(m_ambient_color.y() - .1);
        if(m_ambient_color.y() < 0.)
          m_ambient_color.setY(0.);
        m_ambient_color.setZ(m_ambient_color.z() - .1);
        if(m_ambient_color.z() < 0.)
          m_ambient_color.setZ(0.);
        displayMessage(QString("Light color=(%1 %2 %3).")
                           .arg(m_ambient_color.x())
                           .arg(m_ambient_color.y())
                           .arg(m_ambient_color.z()));
        update();
      } else if((e->key() == ::Qt::Key_PageUp) && (modifiers == ::Qt::ShiftModifier)) {
        m_ambient_color.setX(m_ambient_color.x() + .1);
        if(m_ambient_color.x() > 1.)
          m_ambient_color.setX(1.);
        displayMessage(QString("Light color=(%1 %2 %3).")
                           .arg(m_ambient_color.x())
                           .arg(m_ambient_color.y())
                           .arg(m_ambient_color.z()));
        update();
      } else if((e->key() == ::Qt::Key_PageUp) && (modifiers == ::Qt::AltModifier)) {
        m_ambient_color.setY(m_ambient_color.y() + .1);
        if(m_ambient_color.y() > 1.)
          m_ambient_color.setY(1.);
        displayMessage(QString("Light color=(%1 %2 %3).")
                           .arg(m_ambient_color.x())
                           .arg(m_ambient_color.y())
                           .arg(m_ambient_color.z()));
        update();
      } else if((e->key() == ::Qt::Key_PageUp) && (modifiers == ::Qt::ControlModifier)) {
        m_ambient_color.setZ(m_ambient_color.z() + .1);
        if(m_ambient_color.z() > 1.)
          m_ambient_color.setZ(1.);
        displayMessage(QString("Light color=(%1 %2 %3).")
                           .arg(m_ambient_color.x())
                           .arg(m_ambient_color.y())
                           .arg(m_ambient_color.z()));
        update();
      } else if((e->key() == ::Qt::Key_PageDown) && (modifiers == ::Qt::ShiftModifier)) {
        m_ambient_color.setX(m_ambient_color.x() - .1);
        if(m_ambient_color.x() < 0.)
          m_ambient_color.setX(0.);
        displayMessage(QString("Light color=(%1 %2 %3).")
                           .arg(m_ambient_color.x())
                           .arg(m_ambient_color.y())
                           .arg(m_ambient_color.z()));
        update();
      } else if((e->key() == ::Qt::Key_PageDown) && (modifiers == ::Qt::AltModifier)) {
        m_ambient_color.setY(m_ambient_color.y() - .1);
        if(m_ambient_color.y() < 0.)
          m_ambient_color.setY(0.);
        displayMessage(QString("Light color=(%1 %2 %3).")
                           .arg(m_ambient_color.x())
                           .arg(m_ambient_color.y())
                           .arg(m_ambient_color.z()));
        update();
      } else if((e->key() == ::Qt::Key_PageDown) && (modifiers == ::Qt::ControlModifier)) {
        m_ambient_color.setZ(m_ambient_color.z() - .1);
        if(m_ambient_color.z() < 0.)
          m_ambient_color.setZ(0.);
        displayMessage(QString("Light color=(%1 %2 %3).")
                           .arg(m_ambient_color.x())
                           .arg(m_ambient_color.y())
                           .arg(m_ambient_color.z()));
        update();
      } else if((e->key() == ::Qt::Key_O) && (modifiers == ::Qt::NoButton)) {
        bool old_2D = is_two_dimensional();
        m_no_2D_mode = !m_no_2D_mode;
        if(old_2D != is_two_dimensional()) {
          if(is_two_dimensional()) {
            displayMessage(QString("Viewer is in 2D mode."));
          } else {
            displayMessage(QString("Viewer is in 3D mode."));
          }
          set_camera_mode();
          update();
        }
      } else if((e->key() == ::Qt::Key_V) && (modifiers == ::Qt::ControlModifier)) {
        m_draw_sphere_vertex = !m_draw_sphere_vertex;
        displayMessage(QString("Draw sphere vertex=%1.").arg(m_draw_sphere_vertex ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_E) && (modifiers == ::Qt::ControlModifier)) {
        m_draw_cylinder_edge = !m_draw_cylinder_edge;
        displayMessage(QString("Draw cylinder edge=%1.").arg(m_draw_cylinder_edge ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_N) && (modifiers == ::Qt::ControlModifier)) {
        m_draw_normals = !m_draw_normals;
        displayMessage(QString("Draw normals=%1.").arg(m_draw_normals ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_T) && (modifiers == ::Qt::ControlModifier)) {
        m_draw_mesh_triangles = !m_draw_mesh_triangles;
        displayMessage(QString("Draw triangles=%1.").arg(m_draw_mesh_triangles ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_M) && (modifiers == ::Qt::ControlModifier)) {
        m_use_default_color_normal = !m_use_default_color_normal;
        displayMessage(QString("Normal mono color=%1.").arg(m_use_default_color_normal ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_N) && (modifiers == ::Qt::ShiftModifier)) {
        m_display_face_normal = !m_display_face_normal;
        displayMessage(QString("Display face normal=%1.").arg(m_display_face_normal ? "true" : "false"));
        update();
      } else if((e->key() == ::Qt::Key_F2)) {
        capture_screenshot(QString("./screenshot.png"));
        displayMessage(QString("Screenshot saved in ./screenshot"));
      } else {
        CGAL::QGLViewer::keyPressEvent(e);
      } // By default call QGLViewer key press
    }
  }

  virtual QString helpString() const { return helpString("CGAL Basic Viewer"); }

  virtual QString helpString(const char* title) const {
    QString text(QString("<h2>") + QString(title) + QString("</h2>"));
    text += "Use the mouse to move the camera around the object. ";
    text += "You can respectively revolve around, zoom and translate with "
            "the three mouse buttons. ";
    text += "Left and middle buttons pressed together rotate around the "
            "camera view direction axis<br><br>";
    text += "Pressing <b>Alt</b> and one of the function keys "
            "(<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
    text += "Simply press the function key again to restore it. "
            "Several keyFrames define a ";
    text += "camera path. Paths are saved when you quit the application "
            "and restored at next start.<br><br>";
    text += "Press <b>F</b> to display the frame rate, <b>A</b> for the "
            "world axis, ";
    text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> "
            "to save a snapshot. ";
    text += "See the <b>Keyboard</b> tab in this window for a complete "
            "shortcut list.<br><br>";
    text += "Double clicks automates single click actions: A left button "
            "double click aligns the closer axis with the camera (if close enough). ";
    text += "A middle button double click fits the zoom of the camera and "
            "the right button re-centers the scene.<br><br>";
    text += "A left button double click while holding right button pressed "
            "defines the camera <i>Revolve Around Point</i>. ";
    text += "See the <b>Mouse</b> tab and the documentation web pages for "
            "details.<br><br>";
    text += "Press <b>Escape</b> to exit the viewer.";
    return text;
  }

  void capture_screenshot(const QString& file_path) {
    QScreen* screen;
    screen = QApplication::primaryScreen();

    // auto geom = screen->geometry();
    auto qpx_pixmap = screen->grabWindow(this->winId());
    qpx_pixmap.save(file_path);
  }

public:
  std::function<bool(QKeyEvent*, CGAL::Qt::Basic_viewer*)> on_key_pressed;

protected:
  const Graphics_scene& m_scene;

  bool m_draw_vertices;
  bool m_draw_edges;
  bool m_draw_rays;
  bool m_draw_lines;
  bool m_draw_faces;
  bool m_draw_text;
  bool m_draw_normals;
  bool m_draw_cylinder_edge;
  bool m_draw_sphere_vertex;
  bool m_draw_mesh_triangles;
  bool m_flat_shading;
  bool m_use_default_color;
  bool m_use_default_color_normal;
  bool m_display_face_normal;
  bool m_inverse_normal;
  bool m_no_2D_mode;
  bool m_geometry_feature_enabled;
  bool m_prev_scene_empty;

  enum {
    CLIPPING_PLANE_OFF = 0,
    CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF,
    CLIPPING_PLANE_SOLID_HALF_WIRE_HALF,
    CLIPPING_PLANE_SOLID_HALF_ONLY,
    CLIPPING_PLANE_END_INDEX
  };

  int m_use_clipping_plane = CLIPPING_PLANE_OFF;
  CGAL::qglviewer::ManipulatedFrame* m_frame_plane = nullptr;

  // Buffer for clipping plane is not stored in the scene because it is not
  // filled by users but by the basic viewer.
  std::vector<BufferType> m_array_for_clipping_plane;

  double m_size_vertices = 1.;
  double m_size_edges = 1.;
  double m_size_rays = 1.;
  double m_size_lines = 1.;
  double m_size_normals = .2;
  double m_height_factor_normals = .02;

  CGAL::IO::Color m_default_color_normal;
  QVector4D m_ambient_color;

  bool m_are_buffers_initialized;

  // CGAL::qglviewer::LocalConstraint constraint;
  CGAL::qglviewer::WorldConstraint constraint;

  static const unsigned int NB_GL_BUFFERS = (GS::END_POS - GS::BEGIN_POS) + (GS::END_COLOR - GS::BEGIN_COLOR) +
                                            3; // +2 for normals (mono and color), +1 for clipping plane

  QOpenGLBuffer buffers[NB_GL_BUFFERS]; // +1 for the buffer of clipping plane

  // The following enum gives the indices of the different vao.
  enum { VAO_POINTS = 0, VAO_SEGMENTS, VAO_RAYS, VAO_LINES, VAO_FACES, VAO_CLIPPING_PLANE, NB_VAO_BUFFERS };
  QOpenGLVertexArrayObject vao[NB_VAO_BUFFERS];

  QOpenGLShaderProgram rendering_program_face;
  QOpenGLShaderProgram rendering_program_p_l;
  QOpenGLShaderProgram rendering_program_line;
  QOpenGLShaderProgram rendering_program_clipping_plane;
  QOpenGLShaderProgram rendering_program_sphere;
  QOpenGLShaderProgram rendering_program_cylinder;
  QOpenGLShaderProgram rendering_program_normal;
  QOpenGLShaderProgram rendering_program_triangle;

  // variables for clipping plane
  bool clipping_plane_rendering = true; // will be toggled when alt+c is pressed, which is used for indicating whether
                                        // or not to render the clipping plane ;
  float clipping_plane_rendering_transparency = 0.5f; // to what extent the transparent part should be rendered;
};

//------------------------------------------------------------------------------
class QApplication_and_basic_viewer
{
public:
  QApplication_and_basic_viewer(const CGAL::Graphics_scene& buffer, const char* title = "CGAL Basic Viewer")
      : m_application(nullptr)
      , m_basic_viewer(nullptr)
      , m_argc(1) {
    m_argv[0] = new char[strlen(title) + 1];
    memcpy(m_argv[0], title, strlen(title) + 1);
    m_argv[1] = nullptr;

#if defined(CGAL_TEST_SUITE)
    bool cgal_test_suite = true;
#else
    bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

    if(cgal_test_suite) {
      return;
    }

    Qt::init_ogl_context(4, 3);
    m_application = new QApplication(m_argc, const_cast<char**>(m_argv));
    m_basic_viewer = new Basic_viewer(m_application->activeWindow(), buffer, title);
  }

  ~QApplication_and_basic_viewer() {
    delete[] m_argv[0];
    delete m_basic_viewer;
    delete m_application;
  }

  operator bool() const { return m_application != nullptr; }

  void run() {
    if(m_application != nullptr) {
      m_basic_viewer->show();
      m_application->exec();
    }
  }

  Basic_viewer& basic_viewer() {
    CGAL_assertion(m_basic_viewer != nullptr);
    return *m_basic_viewer;
  }

protected:
  QApplication* m_application;
  Basic_viewer* m_basic_viewer;
  char* m_argv[2];
  int m_argc;
};

} // End namespace Qt

// A shortcut to use directly CGAL::Basic_viewer instead of CGAL::Qt::Basic_viewer.
// Can be changed later if we have several viewers.
using Qt::Basic_viewer;

inline void draw_graphics_scene(const Graphics_scene& graphics_scene, const char* title = "CGAL Basic Viewer (Qt)") {
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if(!cgal_test_suite) {
    Qt::init_ogl_context(4, 3);

    int argc = 1;
    const char* argv[2] = {title, nullptr};
    QApplication app(argc, const_cast<char**>(argv));
    Basic_viewer basic_viewer(app.activeWindow(), graphics_scene, title);

    basic_viewer.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_QT_BASIC_VIEWER_H
