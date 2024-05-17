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

#ifndef CGAL_BASIC_VIEWER_H
#define CGAL_BASIC_VIEWER_H

// TODO #include <CGAL/license/GraphicsView.h>
#include <cstring>
#include <iostream>
#include <tuple>
#include <string>
#include <CGAL/Graphics_scene.h>

#ifdef CGAL_USE_BASIC_VIEWER

#ifdef __GNUC__
#if  __GNUC__ >= 9
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif
#endif

#include <QApplication>
#include <QKeyEvent>

#include <CGAL/Qt/qglviewer.h>
#include <CGAL/Qt/manipulatedFrame.h>
#include <QKeyEvent>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#ifdef __GNUC__
#if __GNUC__ >= 9
#  pragma GCC diagnostic pop
#endif
#endif

#include <vector>
#include <cstdlib>
#include <cfloat>

#include <CGAL/Basic_shaders.h>
#include <CGAL/Qt/CreateOpenGLContext.h>
#include <CGAL/Qt/constraint.h>
#include <CGAL/assertions.h>
#include <CGAL/Qt/init_ogl_context.h>

#define CGAL_BASIC_VIEWER_INIT_SIZE_X 500
#define CGAL_BASIC_VIEWER_INIT_SIZE_Y 450

namespace CGAL {
namespace Qt {
//------------------------------------------------------------------------------
class Basic_viewer : public CGAL::QGLViewer
{
public:
  using BufferType=float;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
  typedef Local_kernel::Point_3  Local_point;
  typedef Local_kernel::Vector_3 Local_vector;
  using GS=Graphics_scene;

  // Constructor/Destructor
  Basic_viewer(QWidget* parent,
               const Graphics_scene& buf,
               const char* title="",
               bool draw_vertices=false,
               bool draw_edges=true,
               bool draw_faces=true,
               bool use_mono_color=false,
               bool inverse_normal=false,
               bool draw_rays=true,
               bool draw_lines=true,
               bool draw_text=true,
               bool no_2D_mode=false) :
    CGAL::QGLViewer(parent),
    gBuffer(buf),
    m_draw_vertices(draw_vertices),
    m_draw_edges(draw_edges),
    m_draw_rays(draw_rays),
    m_draw_lines(draw_lines),
    m_draw_faces(draw_faces),
    m_flatShading(true),
    m_use_mono_color(use_mono_color),
    m_inverse_normal(inverse_normal),
    m_draw_text(draw_text),
    m_no_2D_mode(no_2D_mode),
    m_size_points(7.),
    m_size_edges(3.1),
    m_size_rays(3.1),
    m_size_lines(3.1),
    m_vertices_mono_color(200, 60, 60),
    m_edges_mono_color(0, 0, 0),
    m_rays_mono_color(0, 0, 0),
    m_lines_mono_color(0, 0, 0),
    m_faces_mono_color(60, 60, 200),
    m_ambient_color(0.6f, 0.5f, 0.5f, 0.5f),
    m_are_buffers_initialized(false)
  {
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

    if (title[0]==0)
      setWindowTitle("CGAL Basic Viewer");
    else
      setWindowTitle(title);

    resize(CGAL_BASIC_VIEWER_INIT_SIZE_X, CGAL_BASIC_VIEWER_INIT_SIZE_Y);

    if (inverse_normal)
    { reverse_all_normals(); }
  }

  ~Basic_viewer()
  {
    makeCurrent();
    for (unsigned int i=0; i<NB_GL_BUFFERS; ++i)
      buffers[i].destroy();

    for (unsigned int i=0; i<NB_VAO_BUFFERS; ++i)
      vao[i].destroy();

    Q_FOREACH(QOpenGLShader* shader, rendering_program_p_l.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_face.shaders())
      delete shader;
    Q_FOREACH(QOpenGLShader* shader, rendering_program_clipping_plane.shaders())
      delete shader;
    delete m_frame_plane;
  }

  void draw_vertices(bool b)
  { m_draw_vertices = b; }
  void draw_edges(bool b)
  { m_draw_edges = b; }
  void draw_rays(bool b)
  { m_draw_rays = b; }
  void draw_lines(bool b)
  { m_draw_lines = b; }
  void draw_faces(bool b)
  { m_draw_faces = b; }
  void use_mono_color(bool b)
  { m_use_mono_color = b; }
  void draw_text(bool b)
  { m_draw_text = b; }

  void vertices_mono_color(const CGAL::IO::Color& c)
  { m_vertices_mono_color=c; }
  void edges_mono_color(const CGAL::IO::Color& c)
  { m_edges_mono_color=c; }
  void rays_mono_color(const CGAL::IO::Color& c)
  { m_rays_mono_color=c; }
  void lines_mono_color(const CGAL::IO::Color& c)
  { m_lines_mono_color=c; }
  void faces_mono_color(const CGAL::IO::Color& c)
  { m_faces_mono_color=c; }

  void toggle_draw_vertices()
  { m_draw_vertices = !m_draw_vertices; }
  void toggle_draw_edges()
  { m_draw_edges = !m_draw_edges; }
  void toggle_draw_rays()
  { m_draw_rays = !m_draw_rays; }
  void toggle_draw_lines()
  { m_draw_lines = !m_draw_lines; }
  void toggle_draw_faces()
  { m_draw_faces = !m_draw_faces; }
  void toggle_use_mono_color()
  { m_use_mono_color = !m_use_mono_color; }
  void toggle_draw_text()
  { m_draw_text = !m_draw_text; }

  bool draw_vertices() const
  { return m_draw_vertices; }
  bool draw_edges() const
  { return m_draw_edges; }
  bool draw_rays() const
  { return m_draw_rays; }
  bool draw_lines() const
  { return m_draw_lines; }
  bool draw_faces() const
  { return m_draw_faces; }
  bool use_mono_color() const
  { return m_use_mono_color; }
  bool reverse_normal() const
  { return m_inverse_normal; }
  bool draw_text() const
  { return m_draw_text; }

  const CGAL::IO::Color& vertices_mono_color() const
  { return m_vertices_mono_color; }
  const CGAL::IO::Color& edges_mono_color() const
  { return m_edges_mono_color; }
  const CGAL::IO::Color& rays_mono_color() const
  { return m_rays_mono_color; }
  const CGAL::IO::Color& lines_mono_color() const
  { return m_lines_mono_color; }
  const CGAL::IO::Color& faces_mono_color() const
  { return m_faces_mono_color; }

  Local_kernel::Plane_3 clipping_plane() const
  {
    CGAL::qglviewer::Vec n=m_frame_plane->inverseTransformOf
        (CGAL::qglviewer::Vec(0.f, 0.f, 1.f));
    const CGAL::qglviewer::Vec& pos=m_frame_plane->position();
    return Local_kernel::Plane_3(n[0], n[1], n[2], -n*pos);
  }

  bool clipping_plane_enabled() const
  { return (m_use_clipping_plane!=CLIPPING_PLANE_OFF); }

  const Graphics_scene& graphics_scene() const
  { return gBuffer; }

  virtual void redraw()
  {
    initialize_buffers();
    update();
  }

  void reverse_all_normals()
  {
    m_inverse_normal=!m_inverse_normal;
    gBuffer.reverse_all_normals();
  }

  // Returns true if the data structure lies on a plane
  bool is_two_dimensional() const
  { return !m_no_2D_mode && gBuffer.is_two_dimensional(); }

  virtual void draw()
  {
    glEnable(GL_DEPTH_TEST);

    QMatrix4x4 clipping_mMatrix;
    clipping_mMatrix.setToIdentity();
    for(int i=0; i< 16 ; i++)
    { clipping_mMatrix.data()[i] =  m_frame_plane->matrix()[i]; }
    QVector4D clipPlane = clipping_mMatrix * QVector4D(0.0, 0.0, 1.0, 0.0);
    QVector4D plane_point = clipping_mMatrix * QVector4D(0,0,0,1);
    if(!m_are_buffers_initialized)
    { initialize_buffers(); }

    QColor color;
    attrib_buffers(this);

    if(m_draw_vertices)
    {
      rendering_program_p_l.bind();

      // rendering_mode == -1: draw all
      // rendering_mode == 0: draw inside clipping plane
      // rendering_mode == 1: draw outside clipping plane
      auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
        vao[VAO_MONO_POINTS].bind();
        color.setRgbF((double)m_vertices_mono_color.red()/(double)255,
                      (double)m_vertices_mono_color.green()/(double)255,
                      (double)m_vertices_mono_color.blue()/(double)255);
        rendering_program_p_l.setAttributeValue("color",color);
        rendering_program_p_l.setUniformValue("point_size", GLfloat(m_size_points));
        rendering_program_p_l.setUniformValue("clipPlane", clipPlane);
        rendering_program_p_l.setUniformValue("pointPlane", plane_point);
        rendering_program_p_l.setUniformValue("rendering_mode", rendering_mode);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_MONO_POINTS)));
        vao[VAO_MONO_POINTS].release();

        vao[VAO_COLORED_POINTS].bind();
        if (m_use_mono_color)
        {
          color.setRgbF((double)m_vertices_mono_color.red()/(double)255,
                        (double)m_vertices_mono_color.green()/(double)255,
                      (double)m_vertices_mono_color.blue()/(double)255);
          rendering_program_p_l.disableAttributeArray("color");
          rendering_program_p_l.setAttributeValue("color",color);
        }
        else
        {
          rendering_program_p_l.enableAttributeArray("color");
        }
        rendering_program_p_l.setUniformValue("point_size", GLfloat(m_size_points));
        rendering_program_p_l.setUniformValue("clipPlane", clipPlane);
        rendering_program_p_l.setUniformValue("pointPlane", plane_point);
        rendering_program_p_l.setUniformValue("rendering_mode", rendering_mode);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_COLORED_POINTS)));
        vao[VAO_COLORED_POINTS].release();
      };

      enum {
        DRAW_ALL = -1, // draw all
        DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
        DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
      };

      if (m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY)
      {
        renderer(DRAW_INSIDE_ONLY);
      }
      else
      {
        renderer(DRAW_ALL);
      }

      rendering_program_p_l.release();
    }

    if(m_draw_edges)
    {
      rendering_program_p_l.bind();

      // rendering_mode == -1: draw all
      // rendering_mode == 0: draw inside clipping plane
      // rendering_mode == 1: draw outside clipping plane
      auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {
        vao[VAO_MONO_SEGMENTS].bind();
        color.setRgbF((double)m_edges_mono_color.red()/(double)255,
                      (double)m_edges_mono_color.green()/(double)255,
                      (double)m_edges_mono_color.blue()/(double)255);
        rendering_program_p_l.setAttributeValue("color",color);
        rendering_program_p_l.setUniformValue("clipPlane", clipPlane);
        rendering_program_p_l.setUniformValue("pointPlane", plane_point);
        rendering_program_p_l.setUniformValue("rendering_mode", rendering_mode);
        glLineWidth(m_size_edges);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_MONO_SEGMENTS)));
        vao[VAO_MONO_SEGMENTS].release();

        vao[VAO_COLORED_SEGMENTS].bind();
        if (m_use_mono_color)
        {
          color.setRgbF((double)m_edges_mono_color.red()/(double)255,
                        (double)m_edges_mono_color.green()/(double)255,
                        (double)m_edges_mono_color.blue()/(double)255);
          rendering_program_p_l.disableAttributeArray("color");
          rendering_program_p_l.setAttributeValue("color",color);
        }
        else
        {
          rendering_program_p_l.enableAttributeArray("color");
        }
        rendering_program_p_l.setUniformValue("clipPlane", clipPlane);
        rendering_program_p_l.setUniformValue("pointPlane", plane_point);
        rendering_program_p_l.setUniformValue("rendering_mode", rendering_mode);
        glLineWidth(m_size_edges);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_COLORED_SEGMENTS)));
        vao[VAO_COLORED_SEGMENTS].release();
      };

      enum {
        DRAW_ALL = -1, // draw all
        DRAW_INSIDE_ONLY, // draw only the part inside the clipping plane
        DRAW_OUTSIDE_ONLY // draw only the part outside the clipping plane
      };

      if (m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY)
      {
        renderer(DRAW_INSIDE_ONLY);
      }
      else
      {
        renderer(DRAW_ALL);
      }

      rendering_program_p_l.release();
    }

    if(m_draw_rays)
    {
        rendering_program_p_l.bind();

        vao[VAO_MONO_RAYS].bind();
        color.setRgbF((double)m_rays_mono_color.red()/(double)255,
                      (double)m_rays_mono_color.green()/(double)255,
                      (double)m_rays_mono_color.blue()/(double)255);
        rendering_program_p_l.setAttributeValue("color",color);
        glLineWidth(m_size_rays);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_MONO_RAYS)));
        vao[VAO_MONO_RAYS].release();

        vao[VAO_COLORED_RAYS].bind();
        if (m_use_mono_color)
        {
            color.setRgbF((double)m_rays_mono_color.red()/(double)255,
                          (double)m_rays_mono_color.green()/(double)255,
                          (double)m_rays_mono_color.blue()/(double)255);
            rendering_program_p_l.disableAttributeArray("color");
            rendering_program_p_l.setAttributeValue("color",color);
        }
        else
        {
            rendering_program_p_l.enableAttributeArray("color");
        }
        glLineWidth(m_size_rays);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_COLORED_RAYS)));
        vao[VAO_COLORED_RAYS].release();

        rendering_program_p_l.release();
    }

    if(m_draw_lines)
    {
        rendering_program_p_l.bind();

        vao[VAO_MONO_LINES].bind();
        color.setRgbF((double)m_lines_mono_color.red()/(double)255,
                      (double)m_lines_mono_color.green()/(double)255,
                      (double)m_lines_mono_color.blue()/(double)255);
        rendering_program_p_l.setAttributeValue("color",color);
        glLineWidth(m_size_lines);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_MONO_LINES)));
        vao[VAO_MONO_LINES].release();

        rendering_program_p_l.release();

        vao[VAO_COLORED_LINES].bind();
        if (m_use_mono_color)
        {
            color.setRgbF((double)m_lines_mono_color.red()/(double)255,
                          (double)m_lines_mono_color.green()/(double)255,
                          (double)m_lines_mono_color.blue()/(double)255);
            rendering_program_p_l.disableAttributeArray("color");
            rendering_program_p_l.setAttributeValue("color",color);
        }
        else
        {
            rendering_program_p_l.enableAttributeArray("color");
        }
        glLineWidth(m_size_lines);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_COLORED_LINES)));
        vao[VAO_COLORED_LINES].release();

        rendering_program_p_l.release();
    }

    // Fix Z-fighting by drawing faces at a depth
    GLfloat offset_factor;
    GLfloat offset_units;
    if (is_two_dimensional()) {
      glGetFloatv(GL_POLYGON_OFFSET_FACTOR, &offset_factor);
      glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
      glPolygonOffset(0.1f, 0.9f);
    }

    if (m_draw_faces)
    {
      rendering_program_face.bind();

      // reference: https://stackoverflow.com/questions/37780345/opengl-how-to-create-order-independent-transparency
      // rendering_mode == -1: draw all as solid;
      // rendering_mode == 0: draw solid only;
      // rendering_mode == 1: draw transparent only;
      auto renderer = [this, &color, &clipPlane, &plane_point](float rendering_mode) {

      vao[VAO_MONO_FACES].bind();
      color.setRgbF((double)m_faces_mono_color.red()/(double)255,
                    (double)m_faces_mono_color.green()/(double)255,
                    (double)m_faces_mono_color.blue()/(double)255);
      rendering_program_face.setAttributeValue("color",color);
      rendering_program_face.setUniformValue("rendering_mode", rendering_mode);
      rendering_program_face.setUniformValue("rendering_transparency", clipping_plane_rendering_transparency);
      rendering_program_face.setUniformValue("clipPlane", clipPlane);
      rendering_program_face.setUniformValue("pointPlane", plane_point);
      glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_MONO_FACES)));
      vao[VAO_MONO_FACES].release();

        vao[VAO_COLORED_FACES].bind();
        if (m_use_mono_color)
        {
          color.setRgbF((double)m_faces_mono_color.red()/(double)255,
                        (double)m_faces_mono_color.green()/(double)255,
                        (double)m_faces_mono_color.blue()/(double)255);
          rendering_program_face.disableAttributeArray("color");
          rendering_program_face.setAttributeValue("color",color);
        }
        else
        {
          rendering_program_face.enableAttributeArray("color");
        }
        rendering_program_face.setUniformValue("rendering_mode", rendering_mode);
        rendering_program_face.setUniformValue("rendering_transparency", clipping_plane_rendering_transparency);
        rendering_program_face.setUniformValue("clipPlane", clipPlane);
        rendering_program_face.setUniformValue("pointPlane", plane_point);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(gBuffer.number_of_elements(GS::POS_COLORED_FACES)));
        vao[VAO_COLORED_FACES].release();
      };

      auto renderer_clipping_plane = [this](bool clipping_plane_rendering) {
        if (!isOpenGL_4_3()) return;
        if (!clipping_plane_rendering) return;
        // render clipping plane here
        rendering_program_clipping_plane.bind();
        vao[VAO_CLIPPING_PLANE].bind();
        glLineWidth(0.1f);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>((m_array_for_clipping_plane.size()/3)));
        glLineWidth(1.0f);
        vao[VAO_CLIPPING_PLANE].release();
        rendering_program_clipping_plane.release();
      };

      enum {
        DRAW_SOLID_ALL = -1, // draw all mesh in solid mode
        DRAW_SOLID_HALF, // draw only the mesh inside the clipping plane as solid
        DRAW_TRANSPARENT_HALF // draw only the mesh outside the clipping plane as transparent
      };

      if (m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF)
      {
        // The z-buffer will prevent transparent objects from being displayed behind other transparent objects.
        // Before rendering all transparent objects, disable z-testing first.

        // 1. draw solid first
        renderer(DRAW_SOLID_HALF);

        // 2. draw transparent layer second with back face culling to avoid messy triangles
        glDepthMask(false); //disable z-testing
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CW);
        renderer(DRAW_TRANSPARENT_HALF);

        // 3. draw solid again without culling and blend to make sure the solid mesh is visible
        glDepthMask(true); //enable z-testing
        glDisable(GL_CULL_FACE);
        glDisable(GL_BLEND);
        renderer(DRAW_SOLID_HALF);

        // 4. render clipping plane here
        renderer_clipping_plane(clipping_plane_rendering);
      }
      else if (m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_WIRE_HALF ||
               m_use_clipping_plane == CLIPPING_PLANE_SOLID_HALF_ONLY)
      {
        // 1. draw solid HALF
        renderer(DRAW_SOLID_HALF);

        // 2. render clipping plane here
        renderer_clipping_plane(clipping_plane_rendering);
      }
      else
      {
        // 1. draw solid FOR ALL
        renderer(DRAW_SOLID_ALL);
      }

      if (is_two_dimensional())
        glPolygonOffset(offset_factor, offset_units);

      rendering_program_face.release();
    }

    if (m_draw_text)
    {
      glDisable(GL_LIGHTING);
      for (std::size_t i=0; i<gBuffer.m_texts_size(); ++i)
      {
        auto& m_texts_vec = gBuffer.get_m_texts();
        CGAL::qglviewer::Vec screenPos=camera()->projectedCoordinatesOf
          (CGAL::qglviewer::Vec(std::get<0>(m_texts_vec[i]).x(),
                                std::get<0>(m_texts_vec[i]).y(),
                                std::get<0>(m_texts_vec[i]).z()));

        drawText((int)screenPos[0], (int)screenPos[1],
                 QString(std::get<1>(m_texts_vec[i]).c_str()));
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
  void compile_shaders()
  {
    rendering_program_face.removeAllShaders();
    rendering_program_p_l.removeAllShaders();
    rendering_program_clipping_plane.removeAllShaders();

    // Create the buffers
    for (unsigned int i=0; i<NB_GL_BUFFERS; ++i)
    {
      if(!buffers[i].isCreated() && !buffers[i].create())
      { std::cerr<<"VBO Creation number "<<i<<" FAILED"<<std::endl; }
    }

    for (unsigned int i=0; i<NB_VAO_BUFFERS; ++i)
    {
      if(!vao[i].isCreated() && !vao[i].create())
      { std::cerr<<"VAO Creation number "<<i<<" FAILED"<<std::endl; }
    }

    // Vertices and segments shader

    const char* source_ = isOpenGL_4_3()
        ? vertex_source_p_l
        : vertex_source_p_l_comp;

    QOpenGLShader *vertex_shader_p_l = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_p_l->compileSourceCode(source_))
    { std::cerr<<"Compiling vertex source FAILED"<<std::endl; }

    source_ = isOpenGL_4_3()
        ? fragment_source_p_l
        : fragment_source_p_l_comp;

    QOpenGLShader *fragment_shader_p_l= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader_p_l->compileSourceCode(source_))
    { std::cerr<<"Compiling fragmentsource FAILED"<<std::endl; }

    if(!rendering_program_p_l.addShader(vertex_shader_p_l))
    { std::cerr<<"adding vertex shader FAILED"<<std::endl; }
    if(!rendering_program_p_l.addShader(fragment_shader_p_l))
    { std::cerr<<"adding fragment shader FAILED"<<std::endl; }
    if(!rendering_program_p_l.link())
    { std::cerr<<"linking Program FAILED"<<std::endl; }

    // Faces shader

    source_ = isOpenGL_4_3()
            ? vertex_source_color
            : vertex_source_color_comp;

    QOpenGLShader *vertex_shader_face = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_face->compileSourceCode(source_))
    { std::cerr<<"Compiling vertex source FAILED"<<std::endl; }

    source_ = isOpenGL_4_3()
            ? fragment_source_color
            : fragment_source_color_comp;

    QOpenGLShader *fragment_shader_face= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader_face->compileSourceCode(source_))
    { std::cerr<<"Compiling fragment source FAILED"<<std::endl; }

    if(!rendering_program_face.addShader(vertex_shader_face))
    { std::cerr<<"adding vertex shader FAILED"<<std::endl; }
    if(!rendering_program_face.addShader(fragment_shader_face))
    { std::cerr<<"adding fragment shader FAILED"<<std::endl; }
    if(!rendering_program_face.link())
    { std::cerr<<"linking Program FAILED"<<std::endl; }

    if (isOpenGL_4_3())
    {
      // clipping plane shader
      source_ = vertex_source_clipping_plane;

      QOpenGLShader *vertex_shader_clipping_plane = new QOpenGLShader(QOpenGLShader::Vertex);
      if (!vertex_shader_clipping_plane->compileSourceCode(source_))
      { std::cerr << "Compiling vertex source for clipping plane FAILED" << std::endl; }

      source_ = fragment_source_clipping_plane;

      QOpenGLShader *fragment_shader_clipping_plane = new QOpenGLShader(QOpenGLShader::Fragment);
      if (!fragment_shader_clipping_plane->compileSourceCode(source_))
      { std::cerr << "Compiling fragment source for clipping plane FAILED" << std::endl; }

      if (!rendering_program_clipping_plane.addShader(vertex_shader_clipping_plane))
      { std::cerr << "Adding vertex shader for clipping plane FAILED" << std::endl;}
      if (!rendering_program_clipping_plane.addShader(fragment_shader_clipping_plane))
      { std::cerr << "Adding fragment shader for clipping plane FAILED" << std::endl; }
      if (!rendering_program_clipping_plane.link())
      { std::cerr << "Linking Program for clipping plane FAILED" << std::endl; }

    }

    // source_ = isOpenGL_4_3()
    //         ? vertex_source_clipping_plane
    //         : vertex_source_clipping_plane_comp;

    // QOpenGLShader *vertex_shader_clipping_plane = new QOpenGLShader(QOpenGLShader::Vertex);
    // if (!vertex_shader_clipping_plane->compileSourceCode(source_))
    // { std::cerr << "Compiling vertex source for clipping plane FAILED" << std::endl; }

    // source_ = isOpenGL_4_3()
    //         ? fragment_source_clipping_plane
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
  }

  void initialize_buffers()
  {
    set_camera_mode();
    rendering_program_p_l.bind();

    // 1) POINT SHADER

    // 1.1) Mono points
    vao[VAO_MONO_POINTS].bind();

    unsigned int bufn = 0;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_MONO_POINTS).data(),
                           gBuffer.get_size_of_index(GS::POS_MONO_POINTS));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);
    // TODO/QUESTION can we use double for BufferType?
    // if not remove the template parameter and use float everywhere.

    buffers[bufn].release();

    rendering_program_p_l.disableAttributeArray("color");

    vao[VAO_MONO_POINTS].release();

    // 1.2) Color points
    vao[VAO_COLORED_POINTS].bind();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_COLORED_POINTS).data(),
                           gBuffer.get_size_of_index(GS::POS_COLORED_POINTS));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[bufn].release();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::COLOR_POINTS).data(),
                           gBuffer.get_size_of_index(GS::COLOR_POINTS));
    rendering_program_p_l.enableAttributeArray("color");
    rendering_program_p_l.setAttributeBuffer("color",GL_FLOAT,0,3);
    buffers[bufn].release();

    vao[VAO_COLORED_POINTS].release();

    // 2) SEGMENT SHADER

    // 2.1) Mono segments
    vao[VAO_MONO_SEGMENTS].bind();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_MONO_SEGMENTS).data(),
                           gBuffer.get_size_of_index(GS::POS_MONO_SEGMENTS));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    rendering_program_p_l.disableAttributeArray("color");

    vao[VAO_MONO_SEGMENTS].release();

    // 2.1) Color segments
    vao[VAO_COLORED_SEGMENTS].bind();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_COLORED_SEGMENTS).data(),
                           gBuffer.get_size_of_index(GS::POS_COLORED_SEGMENTS));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::COLOR_SEGMENTS).data(),
                           gBuffer.get_size_of_index(GS::COLOR_SEGMENTS));
    rendering_program_p_l.enableAttributeArray("color");
    rendering_program_p_l.setAttributeBuffer("color",GL_FLOAT,0,3);
    buffers[bufn].release();

    vao[VAO_COLORED_SEGMENTS].release();

    rendering_program_p_l.release();

    // 3) RAYS SHADER

    // 3.1) Mono rays
    vao[VAO_MONO_RAYS].bind();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_MONO_RAYS).data(),
                           gBuffer.get_size_of_index(GS::POS_MONO_RAYS));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeArray("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    rendering_program_p_l.disableAttributeArray("color");

    vao[VAO_MONO_RAYS].release();

    // 3.2) Color rays

    vao[VAO_COLORED_RAYS].bind();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_COLORED_RAYS).data(),
                           gBuffer.get_size_of_index(GS::POS_COLORED_RAYS));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::COLOR_RAYS).data(),
                           gBuffer.get_size_of_index(GS::COLOR_RAYS));
    rendering_program_p_l.enableAttributeArray("color");
    rendering_program_p_l.setAttributeBuffer("color",GL_FLOAT,0,3);
    buffers[bufn].release();

    vao[VAO_COLORED_RAYS].release();

    rendering_program_p_l.release();

    // 4) LINES SHADER
    // 4.1) Mono lines
    vao[VAO_MONO_LINES].bind();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_MONO_LINES).data(),
                           gBuffer.get_size_of_index(GS::POS_MONO_LINES));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeArray("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    rendering_program_p_l.disableAttributeArray("color");

    vao[VAO_MONO_LINES].release();

    // 4.2 Color lines

    vao[VAO_COLORED_LINES].bind();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_COLORED_LINES).data(),
                           gBuffer.get_size_of_index(GS::POS_COLORED_LINES));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::COLOR_LINES).data(),
                           gBuffer.get_size_of_index(GS::COLOR_LINES));
    rendering_program_p_l.enableAttributeArray("color");
    rendering_program_p_l.setAttributeBuffer("color",GL_FLOAT,0,3);
    buffers[bufn].release();

    vao[VAO_COLORED_LINES].release();

    rendering_program_p_l.release();

    // 5) FACE SHADER
    rendering_program_face.bind();

    // 5.1) Mono faces
    vao[VAO_MONO_FACES].bind();

    // 5.1.1) points of the mono faces
    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_MONO_FACES).data(),
                           gBuffer.get_size_of_index(GS::POS_MONO_FACES));
    rendering_program_face.enableAttributeArray("vertex");
    rendering_program_face.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    // 5.1.2) normals of the mono faces
    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    if (m_flatShading)
    {
      buffers[bufn].allocate(gBuffer.get_array_of_index(GS::FLAT_NORMAL_MONO_FACES).data(),
                             gBuffer.get_size_of_index(GS::FLAT_NORMAL_MONO_FACES));
    }
    else
    {
      buffers[bufn].allocate(gBuffer.get_array_of_index(GS::SMOOTH_NORMAL_MONO_FACES).data(),
                             gBuffer.get_size_of_index(GS::SMOOTH_NORMAL_MONO_FACES));
    }
    rendering_program_face.enableAttributeArray("normal");
    rendering_program_face.setAttributeBuffer("normal",GL_FLOAT,0,3);

    buffers[bufn].release();

    // 5.1.3) color of the mono faces
    rendering_program_face.disableAttributeArray("color");
    vao[VAO_MONO_FACES].release();

    // 5.2) Color faces
    vao[VAO_COLORED_FACES].bind();

    // 5.2.1) points of the color faces
    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::POS_COLORED_FACES).data(),
                           gBuffer.get_size_of_index(GS::POS_COLORED_FACES));
    rendering_program_face.enableAttributeArray("vertex");
    rendering_program_face.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    // 5.2.2) normals of the color faces
    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    if (m_flatShading)
    {
      buffers[bufn].allocate(gBuffer.get_array_of_index(GS::FLAT_NORMAL_COLORED_FACES).data(),
                             gBuffer.get_size_of_index(GS::FLAT_NORMAL_COLORED_FACES));
    }
    else
    {
      buffers[bufn].allocate(gBuffer.get_array_of_index(GS::SMOOTH_NORMAL_COLORED_FACES).data(),
                             gBuffer.get_size_of_index(GS::SMOOTH_NORMAL_COLORED_FACES));
    }
    rendering_program_face.enableAttributeArray("normal");
    rendering_program_face.setAttributeBuffer("normal",GL_FLOAT,0,3);

    buffers[bufn].release();

    // 5.2.3) colors of the faces
    ++bufn;
    CGAL_assertion(bufn<NB_GL_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(gBuffer.get_array_of_index(GS::COLOR_FACES).data(),
                           gBuffer.get_size_of_index(GS::COLOR_FACES));
    rendering_program_face.enableAttributeArray("color");
    rendering_program_face.setAttributeBuffer("color",GL_FLOAT,0,3);

    buffers[bufn].release();

    vao[VAO_COLORED_FACES].release();

    rendering_program_face.release();


    // 6) clipping plane shader
    if (isOpenGL_4_3())
    {
      generate_clipping_plane();

      rendering_program_clipping_plane.bind();

      vao[VAO_CLIPPING_PLANE].bind();
      ++bufn;
      CGAL_assertion(bufn < NB_GL_BUFFERS);
      buffers[bufn].bind();
      buffers[bufn].allocate(m_array_for_clipping_plane.data(),
                             static_cast<int>(m_array_for_clipping_plane.size()*sizeof(BufferType)));
      rendering_program_clipping_plane.enableAttributeArray("vertex");
      rendering_program_clipping_plane.setAttributeBuffer("vertex", GL_FLOAT, 0, 3);

      buffers[bufn].release();

      rendering_program_clipping_plane.release();
    }

    m_are_buffers_initialized = true;
  }

  void attrib_buffers(CGAL::QGLViewer* viewer)
  {
    QMatrix4x4 mvpMatrix;
    QMatrix4x4 mvMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(unsigned int i=0; i < 16; i++)
    {
      mvpMatrix.data()[i] = (float)mat[i];
    }
    viewer->camera()->getModelViewMatrix(mat);
    for(unsigned int i=0; i < 16; i++)
    {
      mvMatrix.data()[i] = (float)mat[i];
    }
    // define material
    QVector4D diffuse( 0.9f,
                       0.9f,
                       0.9f,
                       0.9f );

    QVector4D specular( 0.0f,
                        0.0f,
                        0.0f,
                        1.0f );

    CGAL::Bbox_3 bb;
    if (bb==gBuffer.bounding_box()) // Case of "empty" bounding box
    {
      bb=Local_point(CGAL::ORIGIN).bbox();
      bb=bb + Local_point(1,1,1).bbox(); // To avoid a warning from Qglviewer
    }
    else
    { bb=gBuffer.bounding_box(); }

    QVector4D position((bb.xmax()-bb.xmin())/2,
                       (bb.ymax()-bb.ymin())/2,
                       bb.zmax(), 0.0);
    GLfloat shininess =  1.0f;

    rendering_program_face.bind();
    int mvpLocation = rendering_program_face.uniformLocation("mvp_matrix");
    int mvLocation = rendering_program_face.uniformLocation("mv_matrix");
    int lightLocation[5];
    lightLocation[0] = rendering_program_face.uniformLocation("light_pos");
    lightLocation[1] = rendering_program_face.uniformLocation("light_diff");
    lightLocation[2] = rendering_program_face.uniformLocation("light_spec");
    lightLocation[3] = rendering_program_face.uniformLocation("light_amb");
    lightLocation[4] = rendering_program_face.uniformLocation("spec_power");

    rendering_program_face.setUniformValue(lightLocation[0], position);
    rendering_program_face.setUniformValue(lightLocation[1], diffuse);
    rendering_program_face.setUniformValue(lightLocation[2], specular);
    rendering_program_face.setUniformValue(lightLocation[3], m_ambient_color);
    rendering_program_face.setUniformValue(lightLocation[4], shininess);
    rendering_program_face.setUniformValue(mvpLocation, mvpMatrix);
    rendering_program_face.setUniformValue(mvLocation, mvMatrix);
    rendering_program_face.release();

    rendering_program_p_l.bind();
    int mvpLocation2 = rendering_program_p_l.uniformLocation("mvp_matrix");
    rendering_program_p_l.setUniformValue(mvpLocation2, mvpMatrix);
    rendering_program_p_l.release();


    if (isOpenGL_4_3())
    {
      QMatrix4x4 clipping_mMatrix;
      clipping_mMatrix.setToIdentity();
      for(int i=0; i< 16 ; i++)
      { clipping_mMatrix.data()[i] =  m_frame_plane->matrix()[i]; }

      rendering_program_clipping_plane.bind();
      int vpLocation = rendering_program_clipping_plane.uniformLocation("vp_matrix");
      int mLocation = rendering_program_clipping_plane.uniformLocation("m_matrix");
      rendering_program_clipping_plane.setUniformValue(vpLocation, mvpMatrix);
      rendering_program_clipping_plane.setUniformValue(mLocation, clipping_mMatrix);
      rendering_program_clipping_plane.release();
    }
  }

  void set_camera_mode()
  {
    if (is_two_dimensional())
    {
      camera()->setType(CGAL::qglviewer::Camera::ORTHOGRAPHIC);
      //      Camera Constraint:
      constraint.setRotationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::AXIS);
      constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);

      double cx=0., cy=0., cz=0.;
      if (gBuffer.has_zero_x())      { cx=1.; }
      else if (gBuffer.has_zero_y()) { cy=1.; }
      else                           { cz=1.; }

      camera()->setViewDirection(CGAL::qglviewer::Vec(-cx,-cy,-cz));
      constraint.setRotationConstraintDirection(CGAL::qglviewer::Vec(cx, cy, cz));
      camera()->frame()->setConstraint(&constraint);
    }
    else
    {
      camera()->setType(CGAL::qglviewer::Camera::PERSPECTIVE);
      camera()->frame()->setConstraint(nullptr);
    }
  }

  virtual void init()
  {
    set_camera_mode();
    initializeOpenGLFunctions();

    // Light default parameters
    glLineWidth(m_size_edges);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.f,1.f);
    glClearColor(1.0f,1.0f,1.0f,0.0f);
    glDisable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glDisable(GL_POLYGON_SMOOTH_HINT);
    glBlendFunc(GL_ONE, GL_ZERO);
    glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);

    compile_shaders();

    CGAL::Bbox_3 bb;
    if (bb==gBuffer.bounding_box()) // Case of "empty" bounding box
    {
      bb=Local_point(CGAL::ORIGIN).bbox();
      bb=bb + Local_point(1,1,1).bbox(); // To avoid a warning from Qglviewer
    }
    else
    { bb=gBuffer.bounding_box(); }
    this->camera()->setSceneBoundingBox(CGAL::qglviewer::Vec(bb.xmin(),
                                                       bb.ymin(),
                                                       bb.zmin()),
                                        CGAL::qglviewer::Vec(bb.xmax(),
                                                       bb.ymax(),
                                                       bb.zmax()));

    m_frame_plane=new CGAL::qglviewer::ManipulatedFrame;

    this->showEntireScene();
  }

  void generate_clipping_plane()
  {
    qreal size=((gBuffer.bounding_box().xmax()-gBuffer.bounding_box().xmin()) +
                (gBuffer.bounding_box().ymax()-gBuffer.bounding_box().ymin()) +
                (gBuffer.bounding_box().zmax()-gBuffer.bounding_box().zmin()));
    const unsigned int nbSubdivisions=30;

    auto& array = m_array_for_clipping_plane;
    array.clear();
    for (unsigned int i=0; i<=nbSubdivisions; ++i)
    {
      const float pos = float(size*(2.0*i/nbSubdivisions-1.0));
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
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    if(!on_key_pressed || !on_key_pressed(e, this))
    {
      const ::Qt::KeyboardModifiers modifiers = e->modifiers();
      if ((e->key()==::Qt::Key_C) && (modifiers==::Qt::NoButton))
      {
        if (!isOpenGL_4_3()) return;
        if (!is_two_dimensional())
        {
          // toggle clipping plane
          m_use_clipping_plane = (m_use_clipping_plane + 1) % CLIPPING_PLANE_END_INDEX;
          if (m_use_clipping_plane==CLIPPING_PLANE_OFF)
          { setManipulatedFrame(nullptr); }
          else
          { setManipulatedFrame(m_frame_plane); }

          switch(m_use_clipping_plane)
          {
          case CLIPPING_PLANE_OFF: displayMessage(QString("Draw clipping = false")); break;
          case CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF: clipping_plane_rendering=true; displayMessage(QString("Draw clipping = solid half & transparent half")); break;
          case CLIPPING_PLANE_SOLID_HALF_WIRE_HALF: displayMessage(QString("Draw clipping = solid half & wireframe half")); break;
          case CLIPPING_PLANE_SOLID_HALF_ONLY: displayMessage(QString("Draw clipping = solid half only")); break;
          default: break;
          }
          update();
        }
      }

      else if ((e->key()==::Qt::Key_C) && (modifiers==::Qt::AltModifier))
      {
        if (!isOpenGL_4_3()) return;
        if (m_use_clipping_plane!=CLIPPING_PLANE_OFF)
        {
          clipping_plane_rendering = !clipping_plane_rendering;
          displayMessage(QString("Draw clipping plane=%1.").arg(clipping_plane_rendering?"true":"false"));
          update();
        }
      }
      else if ((e->key()==::Qt::Key_E) && (modifiers==::Qt::NoButton))
      {
        m_draw_edges=!m_draw_edges;
        displayMessage(QString("Draw edges=%1.").arg(m_draw_edges?"true":"false"));
        update();
      }
      else if ((e->key()==::Qt::Key_M) && (modifiers==::Qt::NoButton))
      {
        m_use_mono_color=!m_use_mono_color;
        displayMessage(QString("Mono color=%1.").arg(m_use_mono_color?"true":"false"));
        update();
      }
      else if ((e->key()==::Qt::Key_N) && (modifiers==::Qt::NoButton))
      {
        reverse_all_normals();
        displayMessage(QString("Inverse normal=%1.").arg(m_inverse_normal?"true":"false"));
        redraw();
      }
      else if ((e->key()==::Qt::Key_S) && (modifiers==::Qt::NoButton))
      {
        m_flatShading=!m_flatShading;
        if (m_flatShading)
          displayMessage("Flat shading.");
        else
          displayMessage("Gouraud shading.");
        redraw();
      }
      else if ((e->key()==::Qt::Key_T) && (modifiers==::Qt::NoButton))
      {
        m_draw_text=!m_draw_text;
        displayMessage(QString("Draw text=%1.").arg(m_draw_text?"true":"false"));
        update();
      }
      else if ((e->key()==::Qt::Key_U) && (modifiers==::Qt::NoButton))
      {
        if (is_two_dimensional())
        {
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
        //camera()->frame()->setConstraint(NULL);
        // camera()->setViewDirection(CGAL::qglviewer::Vec(-cx,-cy,-cz));
        //constraint.setRotationConstraintDirection(CGAL::qglviewer::Vec(cx, cy, cz));
        //camera()->frame()->setConstraint(&constraint);
        //update();
          redraw();
        }
      }
      else if ((e->key()==::Qt::Key_V) && (modifiers==::Qt::NoButton))
      {
        m_draw_vertices=!m_draw_vertices;
        displayMessage(QString("Draw vertices=%1.").arg(m_draw_vertices?"true":"false"));
        update();
      }
      else if ((e->key()==::Qt::Key_W) && (modifiers==::Qt::NoButton))
      {
        m_draw_faces=!m_draw_faces;
        displayMessage(QString("Draw faces=%1.").arg(m_draw_faces?"true":"false"));
        update();
      }
      else if ((e->key()==::Qt::Key_Plus) && (!modifiers.testFlag(::Qt::ControlModifier))) // No ctrl
      {
        m_size_edges+=.5;
        displayMessage(QString("Size of edges=%1.").arg(m_size_edges));
        update();
      }
      else if ((e->key()==::Qt::Key_Minus) && (!modifiers.testFlag(::Qt::ControlModifier))) // No ctrl
      {
        if (m_size_edges>.5) m_size_edges-=.5;
        displayMessage(QString("Size of edges=%1.").arg(m_size_edges));
        update();
      }
      else if ((e->key()==::Qt::Key_Plus) && (modifiers.testFlag(::Qt::ControlModifier)))
      {
        m_size_points+=.5;
        displayMessage(QString("Size of points=%1.").arg(m_size_points));
        update();
      }
      else if ((e->key()==::Qt::Key_Minus) && (modifiers.testFlag(::Qt::ControlModifier)))
      {
        if (m_size_points>.5) m_size_points-=.5;
        displayMessage(QString("Size of points=%1.").arg(m_size_points));
        update();
      }
      else if ((e->key()==::Qt::Key_PageUp) && (modifiers==::Qt::NoButton))
      {
        m_ambient_color.setX(m_ambient_color.x()+.1);
        if (m_ambient_color.x()>1.) m_ambient_color.setX(1.);
        m_ambient_color.setY(m_ambient_color.x()+.1);
        if (m_ambient_color.y()>1.) m_ambient_color.setY(1.);
        m_ambient_color.setZ(m_ambient_color.x()+.1);
        if (m_ambient_color.z()>1.) m_ambient_color.setZ(1.);
        displayMessage(QString("Light color=(%1 %2 %3).").
                       arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
        update();
      }
      else if ((e->key()==::Qt::Key_PageDown) && (modifiers==::Qt::NoButton))
      {
        m_ambient_color.setX(m_ambient_color.x()-.1);
        if (m_ambient_color.x()<0.) m_ambient_color.setX(0.);
        m_ambient_color.setY(m_ambient_color.y()-.1);
        if (m_ambient_color.y()<0.) m_ambient_color.setY(0.);
        m_ambient_color.setZ(m_ambient_color.z()-.1);
        if (m_ambient_color.z()<0.) m_ambient_color.setZ(0.);
        displayMessage(QString("Light color=(%1 %2 %3).").
                       arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
        update();
      }
      else if ((e->key()==::Qt::Key_PageUp) && (modifiers==::Qt::ShiftModifier))
      {
        m_ambient_color.setX(m_ambient_color.x()+.1);
        if (m_ambient_color.x()>1.) m_ambient_color.setX(1.);
        displayMessage(QString("Light color=(%1 %2 %3).").
                       arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
        update();
      }
      else if ((e->key()==::Qt::Key_PageUp) && (modifiers==::Qt::AltModifier))
      {
        m_ambient_color.setY(m_ambient_color.y()+.1);
        if (m_ambient_color.y()>1.) m_ambient_color.setY(1.);
        displayMessage(QString("Light color=(%1 %2 %3).").
                       arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
        update();
      }
      else if ((e->key()==::Qt::Key_PageUp) && (modifiers==::Qt::ControlModifier))
      {
        m_ambient_color.setZ(m_ambient_color.z()+.1);
        if (m_ambient_color.z()>1.) m_ambient_color.setZ(1.);
        displayMessage(QString("Light color=(%1 %2 %3).").
                       arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
        update();
      }
      else if ((e->key()==::Qt::Key_PageDown) && (modifiers==::Qt::ShiftModifier))
      {
        m_ambient_color.setX(m_ambient_color.x()-.1);
        if (m_ambient_color.x()<0.) m_ambient_color.setX(0.);
        displayMessage(QString("Light color=(%1 %2 %3).").
                       arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
        update();
      }
      else if ((e->key()==::Qt::Key_PageDown) && (modifiers==::Qt::AltModifier))
      {
        m_ambient_color.setY(m_ambient_color.y()-.1);
        if (m_ambient_color.y()<0.) m_ambient_color.setY(0.);
        displayMessage(QString("Light color=(%1 %2 %3).").
                       arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
        update();
      }
      else if ((e->key()==::Qt::Key_PageDown) && (modifiers==::Qt::ControlModifier))
      {
        m_ambient_color.setZ(m_ambient_color.z()-.1);
        if (m_ambient_color.z()<0.) m_ambient_color.setZ(0.);
        displayMessage(QString("Light color=(%1 %2 %3).").
                       arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
        update();
      }
      else if ((e->key()==::Qt::Key_O) && (modifiers==::Qt::NoButton))
      {
        bool old_2D=is_two_dimensional();
        m_no_2D_mode=!m_no_2D_mode;
        if (old_2D!=is_two_dimensional())
        {
          if (is_two_dimensional())
          { displayMessage(QString("Viewer is in 2D mode.")); }
          else { displayMessage(QString("Viewer is in 3D mode.")); }
          set_camera_mode();
          update();
        }
      }
      else
      { CGAL::QGLViewer::keyPressEvent(e); } // By default call QGLViewer key press
    }
  }

  virtual QString helpString() const
  { return helpString("CGAL Basic Viewer"); }

  virtual QString helpString(const char* title) const
  {
    QString text(QString("<h2>")+QString(title)+QString("</h2>"));
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
public:
  std::function<bool(QKeyEvent *, CGAL::Qt::Basic_viewer *)> on_key_pressed;

protected:
  const Graphics_scene& gBuffer;

  bool m_draw_vertices;
  bool m_draw_edges;
  bool m_draw_rays;
  bool m_draw_lines;
  bool m_draw_faces;
  bool m_flatShading;
  bool m_use_mono_color;
  bool m_inverse_normal;
  bool m_draw_text;
  bool m_no_2D_mode;

  enum {
    CLIPPING_PLANE_OFF = 0,
    CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF,
    CLIPPING_PLANE_SOLID_HALF_WIRE_HALF,
    CLIPPING_PLANE_SOLID_HALF_ONLY,
    CLIPPING_PLANE_END_INDEX
  };

  int m_use_clipping_plane=CLIPPING_PLANE_OFF;
  CGAL::qglviewer::ManipulatedFrame* m_frame_plane=nullptr;

  // Buffer for clipping plane is not stored in the scene because it is not
  // filled by users but by the basic viewer.
  std::vector<BufferType> m_array_for_clipping_plane;

  double m_size_points;
  double m_size_edges;
  double m_size_rays;
  double m_size_lines;

  CGAL::IO::Color m_vertices_mono_color;
  CGAL::IO::Color m_edges_mono_color;
  CGAL::IO::Color m_rays_mono_color;
  CGAL::IO::Color m_lines_mono_color;
  CGAL::IO::Color m_faces_mono_color;
  QVector4D   m_ambient_color;

  bool m_are_buffers_initialized;

  // CGAL::qglviewer::LocalConstraint constraint;
  CGAL::qglviewer::WorldConstraint constraint;

  static const unsigned int NB_GL_BUFFERS=(GS::END_POS-GS::BEGIN_POS)+
    (GS::END_COLOR-GS::BEGIN_COLOR)+3; // +2 for normals (mono and color), +1 for clipping plane

  QOpenGLBuffer buffers[NB_GL_BUFFERS]; // +1 for the buffer of clipping plane

  // The following enum gives the indices of the different vao.
  enum
    { VAO_MONO_POINTS=0,
      VAO_COLORED_POINTS,
      VAO_MONO_SEGMENTS,
      VAO_COLORED_SEGMENTS,
      VAO_MONO_RAYS,
      VAO_COLORED_RAYS,
      VAO_MONO_LINES,
      VAO_COLORED_LINES,
      VAO_MONO_FACES,
      VAO_COLORED_FACES,
      VAO_CLIPPING_PLANE,
      NB_VAO_BUFFERS
    };
  QOpenGLVertexArrayObject vao[NB_VAO_BUFFERS];

  QOpenGLShaderProgram rendering_program_face;
  QOpenGLShaderProgram rendering_program_p_l;
  QOpenGLShaderProgram rendering_program_clipping_plane;

  // variables for clipping plane
  bool clipping_plane_rendering = true; // will be toggled when alt+c is pressed, which is used for indicating whether or not to render the clipping plane ;
  float clipping_plane_rendering_transparency = 0.5f; // to what extent the transparent part should be rendered;

};

//------------------------------------------------------------------------------
class QApplication_and_basic_viewer
{
public:
  QApplication_and_basic_viewer(const CGAL::Graphics_scene& buffer,
                                const char* title="CGAL Basic Viewer"):
    m_application(nullptr),
    m_basic_viewer(nullptr),
    m_argc(1)
  {
    m_argv[0]=new char[strlen(title)+1];
    memcpy(m_argv[0], title, strlen(title)+1);
    m_argv[1]=nullptr;

#if defined(CGAL_TEST_SUITE)
    bool cgal_test_suite = true;
#else
    bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

    if (cgal_test_suite)
    { return; }

    Qt::init_ogl_context(4, 3);
    m_application=new QApplication(m_argc, const_cast<char **>(m_argv));
    m_basic_viewer=new Basic_viewer(m_application->activeWindow(),
                                    buffer, title);
  }

  ~QApplication_and_basic_viewer()
  {
    delete[] m_argv[0];
    delete m_basic_viewer;
    delete m_application;
  }

  operator bool() const
  { return m_application!=nullptr; }

  void run()
  {
    if (m_application!=nullptr)
    {
      m_basic_viewer->show();
      m_application->exec();
    }
  }

  Basic_viewer& basic_viewer()
  {
    CGAL_assertion(m_basic_viewer!=nullptr);
    return *m_basic_viewer;
  }

protected:
  QApplication* m_application;
  Basic_viewer* m_basic_viewer;
  char *m_argv[2];
  int m_argc;
};

} // End namespace Qt

// A shortcut to use directly CGAL::Basic_viewer instead of CGAL::Qt::Basic_viewer.
// Can be changed later if we have several viewers.
using Qt::Basic_viewer;

inline
void draw_graphics_scene(const Graphics_scene& graphics_scene,
                         const char *title="CGAL Basic Viewer")
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    Qt::init_ogl_context(4, 3);

    int argc = 1;
    const char *argv[2] = {title, nullptr};
    QApplication app(argc, const_cast<char **>(argv));
    Basic_viewer basic_viewer(app.activeWindow(), graphics_scene, title);

    basic_viewer.show();
    app.exec();
  }
}

} // End namespace CGAL

#else // CGAL_USE_BASIC_VIEWER

namespace CGAL
{

  template<class ... T>
  void draw(T...)
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }

  inline
  void draw_graphics_scene(const Graphics_scene&,
                           const char* ="CGAL Basic Viewer")
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_BASIC_VIEWER_H
