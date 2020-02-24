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

#ifndef CGAL_BASIC_VIEWER_QT_H
#define CGAL_BASIC_VIEWER_QT_H

#include <CGAL/license/GraphicsView.h>
#include <iostream>
#include <tuple>
#include <string>

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
#include <QKeyEvent>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QGLBuffer>
#include <QOpenGLShaderProgram>

#ifdef __GNUC__
#if __GNUC__ >= 9
#  pragma GCC diagnostic pop
#endif
#endif

#include <vector>
#include <cstdlib>
#include <cfloat>

#include <CGAL/Buffer_for_vao.h>
#include <CGAL/Qt/CreateOpenGLContext.h>
#include <CGAL/Qt/constraint.h>
#include <CGAL/Random.h>

namespace CGAL
{

//------------------------------------------------------------------------------
const char vertex_source_color[] =
  {
    "#version 120 \n"
    "attribute highp vec4 vertex;\n"
    "attribute highp vec3 normal;\n"
    "attribute highp vec3 color;\n"

    "uniform highp mat4 mvp_matrix;\n"
    "uniform highp mat4 mv_matrix; \n"

    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "varying highp vec4 fColor; \n"

    "uniform highp float point_size; \n"
    "void main(void)\n"
    "{\n"
    "   fP = mv_matrix * vertex; \n"
    "   fN = mat3(mv_matrix)* normal; \n"
    "   fColor = vec4(color, 1.0); \n"
    "   gl_PointSize = point_size;\n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
  };

const char fragment_source_color[] =
  {
    "#version 120 \n"
    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "varying highp vec4 fColor; \n"
    "uniform highp vec4 light_pos;  \n"
    "uniform highp vec4 light_diff; \n"
    "uniform highp vec4 light_spec; \n"
    "uniform highp vec4 light_amb;  \n"
    "uniform float spec_power ; \n"

    "void main(void) { \n"
    "   highp vec3 L = light_pos.xyz - fP.xyz; \n"
    "   highp vec3 V = -fP.xyz; \n"
    
    "   highp vec3 N = normalize(fN); \n"
    "   L = normalize(L); \n"
    "   V = normalize(V); \n"
    
    "   highp vec3 R = reflect(-L, N); \n"
    "   highp vec4 diffuse = max(dot(N,L), 0.0) * light_diff * fColor; \n"
    "   highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"
    "   gl_FragColor = light_amb*fColor + diffuse  ; \n"
    "} \n"
    "\n"
  };

const char vertex_source_p_l[] =
  {
    "#version 120 \n"
    "attribute highp vec4 vertex;\n"
    "attribute highp vec3 color;\n"
    "uniform highp mat4 mvp_matrix;\n"
    "varying highp vec4 fColor; \n"
    "uniform highp float point_size; \n"
    "void main(void)\n"
    "{\n"
    "   gl_PointSize = point_size;\n"
    "   fColor = vec4(color, 1.0); \n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
  };

const char fragment_source_p_l[] =
  {
    "#version 120 \n"
    "varying highp vec4 fColor; \n"
    "void main(void) { \n"
    "gl_FragColor = fColor; \n"
    "} \n"
    "\n"
  };

//------------------------------------------------------------------------------
//  compatibility shaders

const char vertex_source_color_comp[] =
  {
    "attribute highp vec4 vertex;\n"
    "attribute highp vec3 normal;\n"
    "attribute highp vec3 color;\n"
    
    "uniform highp mat4 mvp_matrix;\n"
    "uniform highp mat4 mv_matrix; \n"
    
    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "varying highp vec4 fColor; \n"
    
    "uniform highp float point_size; \n"
    "void main(void)\n"
    "{\n"
    "   fP = mv_matrix * vertex; \n"
    "   highp mat3 mv_matrix_3; \n"
    "   mv_matrix_3[0] = mv_matrix[0].xyz; \n"
    "   mv_matrix_3[1] = mv_matrix[1].xyz; \n"
    "   mv_matrix_3[2] = mv_matrix[2].xyz; \n"
    "   fN = mv_matrix_3* normal;  \n"
    "   fColor = vec4(color, 1.0); \n"
    "   gl_PointSize = point_size;\n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
  };

const char fragment_source_color_comp[] =
  {
    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "varying highp vec4 fColor; \n"
    "uniform highp vec4 light_pos;  \n"
    "uniform highp vec4 light_diff; \n"
    "uniform highp vec4 light_spec; \n"
    "uniform highp vec4 light_amb;  \n"
    "uniform highp float spec_power ; \n"
    
    "void main(void) { \n"
    
    "   highp vec3 L = light_pos.xyz - fP.xyz; \n"
    "   highp vec3 V = -fP.xyz; \n"
    
    "   highp vec3 N = normalize(fN); \n"
    "   L = normalize(L); \n"
    "   V = normalize(V); \n"
    
    "   highp vec3 R = reflect(-L, N); \n"
    "   highp vec4 diffuse = max(dot(N,L), 0.0) * light_diff * fColor; \n"
    "   highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"
    
    "gl_FragColor = light_amb*fColor + diffuse  ; \n"
    "} \n"
    "\n"
  };

const char vertex_source_p_l_comp[] =
  {
    "attribute highp vec4 vertex;\n"
    "attribute highp vec3 color;\n"
    "uniform highp mat4 mvp_matrix;\n"
    "varying highp vec4 fColor; \n"
    "uniform highp float point_size; \n"
    "void main(void)\n"
    "{\n"
    "   gl_PointSize = point_size;\n"
    "   fColor = vec4(color, 1.0); \n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
  };

const char fragment_source_p_l_comp[] =
  {
    "varying highp vec4 fColor; \n"
    "void main(void) { \n"
    "gl_FragColor = fColor; \n"
    "} \n"
    "\n"
  };



//------------------------------------------------------------------------------
inline CGAL::Color get_random_color(CGAL::Random& random)
{
  CGAL::Color res;
  do
  {
    res=CGAL::Color(random.get_int(0,256),
                    random.get_int(0,256),
                    random.get_int(0,256));
  }
  while(res.red()==255 && res.green()==255 && res.blue()==255);
  return res;
}
//------------------------------------------------------------------------------
class Basic_viewer_qt : public CGAL::QGLViewer
{
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
  typedef Local_kernel::Point_3  Local_point;
  typedef Local_kernel::Vector_3 Local_vector;

  // Constructor/Destructor
  Basic_viewer_qt(QWidget* parent,
                  const char* title="",
                  bool draw_vertices=false,
                  bool draw_edges=true,
                  bool draw_faces=true,
                  bool use_mono_color=false,
                  bool inverse_normal=false,
                  bool draw_rays=true,
                  bool draw_lines=true,
                  bool draw_text=true) :
    CGAL::QGLViewer(parent),
    m_draw_vertices(draw_vertices),
    m_draw_edges(draw_edges),
    m_draw_rays(draw_rays),
    m_draw_lines(draw_lines),
    m_draw_faces(draw_faces),
    m_flatShading(true),
    m_use_mono_color(use_mono_color),
    m_inverse_normal(inverse_normal),
    m_draw_text(draw_text),
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
    m_are_buffers_initialized(false),
    m_buffer_for_mono_points(&arrays[POS_MONO_POINTS],
                             nullptr,
                             &m_bounding_box,
                             nullptr, nullptr, nullptr),
    m_buffer_for_colored_points(&arrays[POS_COLORED_POINTS],
                                nullptr,
                                &m_bounding_box,
                                &arrays[COLOR_POINTS],
                                nullptr, nullptr),
    m_buffer_for_mono_segments(&arrays[POS_MONO_SEGMENTS],
                               nullptr,
                               &m_bounding_box,
                               nullptr, nullptr, nullptr),
    m_buffer_for_colored_segments(&arrays[POS_COLORED_SEGMENTS],
                                  nullptr,
                                  &m_bounding_box,
                                  &arrays[COLOR_SEGMENTS],
                                  nullptr, nullptr),
    m_buffer_for_mono_rays(&arrays[POS_MONO_RAYS],
                            nullptr,
                            &m_bounding_box,
                            nullptr, nullptr),
    m_buffer_for_colored_rays(&arrays[POS_COLORED_RAYS],
                               nullptr,
                               &m_bounding_box,
                               &arrays[COLOR_RAYS],
                               nullptr, nullptr),
    m_buffer_for_mono_lines(&arrays[POS_MONO_RAYS],
                           nullptr,
                           &m_bounding_box,
                           nullptr, nullptr),
    m_buffer_for_colored_lines(&arrays[POS_COLORED_LINES],
                              nullptr,
                              &m_bounding_box,
                              &arrays[COLOR_LINES],
                              nullptr, nullptr),
    m_buffer_for_mono_faces(&arrays[POS_MONO_FACES],
                            nullptr,
                            &m_bounding_box,
                            nullptr,
                            &arrays[FLAT_NORMAL_MONO_FACES],
                            &arrays[SMOOTH_NORMAL_MONO_FACES]),
    m_buffer_for_colored_faces(&arrays[POS_COLORED_FACES],
                               nullptr,
                               &m_bounding_box,
                               &arrays[COLOR_FACES],
                               &arrays[FLAT_NORMAL_COLORED_FACES],
                               &arrays[SMOOTH_NORMAL_COLORED_FACES])
  {
    if (title[0]==0)
      setWindowTitle("CGAL Basic Viewer");
    else
      setWindowTitle(title);

    resize(500, 450);

    if (inverse_normal)
    { negate_all_normals(); }
  }

  ~Basic_viewer_qt()
  {
    for (unsigned int i=0; i<NB_VBO_BUFFERS; ++i)
      buffers[i].destroy();

    for (unsigned int i=0; i<NB_VAO_BUFFERS; ++i)
      vao[i].destroy();
  }

  void clear()
  {
    for (unsigned int i=0; i<LAST_INDEX; ++i)
    { arrays[i].clear(); }

    m_bounding_box=CGAL::Bbox_3();
    m_texts.clear();
  }

  bool is_empty() const
  {
    return (m_buffer_for_mono_points.is_empty() &&
            m_buffer_for_colored_points.is_empty() &&
            m_buffer_for_mono_segments.is_empty() &&
            m_buffer_for_colored_segments.is_empty() &&
            m_buffer_for_mono_rays.is_empty() &&
            m_buffer_for_colored_rays.is_empty() &&
            m_buffer_for_mono_lines.is_empty() &&
            m_buffer_for_colored_lines.is_empty() &&
            m_buffer_for_mono_faces.is_empty() &&
            m_buffer_for_colored_faces.is_empty());
  }

  const CGAL::Bbox_3& bounding_box() const
  { return m_bounding_box; }

  bool has_zero_x() const
  {
    return
      m_buffer_for_mono_points.has_zero_x() && 
      m_buffer_for_colored_points.has_zero_x() &&
      m_buffer_for_mono_segments.has_zero_x() &&
      m_buffer_for_colored_segments.has_zero_x() &&
      m_buffer_for_mono_faces.has_zero_x() &&
      m_buffer_for_colored_faces.has_zero_x() &&
      m_buffer_for_mono_rays.has_zero_x() &&
      m_buffer_for_colored_rays.has_zero_x() &&
      m_buffer_for_mono_lines.has_zero_x() &&
      m_buffer_for_colored_lines.has_zero_x();
  }  

  bool has_zero_y() const
  {
    return
      m_buffer_for_mono_points.has_zero_y() && 
      m_buffer_for_colored_points.has_zero_y() &&
      m_buffer_for_mono_segments.has_zero_y() &&
      m_buffer_for_colored_segments.has_zero_y() &&
      m_buffer_for_mono_faces.has_zero_y() &&
      m_buffer_for_colored_faces.has_zero_y() &&
      m_buffer_for_mono_rays.has_zero_y() &&
      m_buffer_for_colored_rays.has_zero_y() &&
      m_buffer_for_mono_lines.has_zero_y() &&
      m_buffer_for_colored_lines.has_zero_y();
  }
  
  bool has_zero_z() const
  {
    return
      m_buffer_for_mono_points.has_zero_z() && 
      m_buffer_for_colored_points.has_zero_z() &&
      m_buffer_for_mono_segments.has_zero_z() &&
      m_buffer_for_colored_segments.has_zero_z() &&
      m_buffer_for_mono_faces.has_zero_z() &&
      m_buffer_for_colored_faces.has_zero_z() &&
      m_buffer_for_mono_rays.has_zero_z() &&
      m_buffer_for_colored_rays.has_zero_z() &&
      m_buffer_for_mono_lines.has_zero_z() &&
      m_buffer_for_colored_lines.has_zero_z();
  }
  
  template<typename KPoint>
  void add_point(const KPoint& p)
  { m_buffer_for_mono_points.add_point(p); }

  template<typename KPoint>
  void add_point(const KPoint& p, const CGAL::Color& acolor)
  { m_buffer_for_colored_points.add_point(p, acolor); }

  template<typename KPoint>
  void add_segment(const KPoint& p1, const KPoint& p2)
  { m_buffer_for_mono_segments.add_segment(p1, p2); }

  template<typename KPoint>
  void add_segment(const KPoint& p1, const KPoint& p2,
                   const CGAL::Color& acolor)
  { m_buffer_for_colored_segments.add_segment(p1, p2, acolor); }

  template <typename KPoint, typename KVector>
  void update_bounding_box_for_ray(const KPoint &p, const KVector &v)
  {
    Local_point lp = get_local_point(p);
    Local_vector lv = get_local_vector(v);
    CGAL::Bbox_3 b = (lp + lv).bbox();
    m_bounding_box += b;
  }

  template <typename KPoint, typename KVector>
  void update_bounding_box_for_line(const KPoint &p, const KVector &v,
                                    const KVector &pv)
  {
    Local_point lp = get_local_point(p);
    Local_vector lv = get_local_vector(v);
    Local_vector lpv = get_local_vector(pv);

    CGAL::Bbox_3 b = lp.bbox() + (lp + lv).bbox() + (lp + lpv).bbox();
    m_bounding_box += b;
  }

  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v)
  {
    double bigNumber = 1e30;
    m_buffer_for_mono_rays.add_ray_segment(p, (p + (bigNumber)*v));
  }

  template <typename KPoint, typename KVector>
  void add_ray(const KPoint &p, const KVector &v, const CGAL::Color &acolor)
  {
    double bigNumber = 1e30;
    m_buffer_for_colored_rays.add_ray_segment(p, (p + (bigNumber)*v), acolor);
  }

  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v)
  {
    double bigNumber = 1e30;
    m_buffer_for_mono_lines.add_line_segment((p - (bigNumber)*v),
                                             (p + (bigNumber)*v));
  }

  template <typename KPoint, typename KVector>
  void add_line(const KPoint &p, const KVector &v, const CGAL::Color &acolor)
  {
    double bigNumber = 1e30;
    m_buffer_for_colored_lines.add_line_segment((p - (bigNumber)*v),
                                                (p + (bigNumber)*v), acolor);
  }

  template<typename KPoint>
  void add_text(const KPoint& kp, const QString& txt)
  {
    Local_point p=get_local_point(kp);
    m_texts.push_back(std::make_tuple(p, txt));
  }

  template<typename KPoint>
  void add_text(const KPoint& kp, const char* txt)
  { add_text(kp, QString(txt)); }

  template<typename KPoint>
  void add_text(const KPoint& kp, const std::string& txt)
  { add_text(kp, txt.c_str()); }

  bool is_a_face_started() const
  {
    return m_buffer_for_mono_faces.is_a_face_started() ||
      m_buffer_for_colored_faces.is_a_face_started();
  }

  void face_begin()
  {
    if (is_a_face_started())
    {
      std::cerr<<"You cannot start a new face before to finish the previous one."<<std::endl;
    }
    else
    { m_buffer_for_mono_faces.face_begin(); }
  }

  void face_begin(const CGAL::Color& acolor)
  {
    if (is_a_face_started())
    {
      std::cerr<<"You cannot start a new face before to finish the previous one."<<std::endl;
    }
    else
    { m_buffer_for_colored_faces.face_begin(acolor); }
  }

  template<typename KPoint>
  bool add_point_in_face(const KPoint& kp)
  {
    if (m_buffer_for_mono_faces.is_a_face_started())
    { return m_buffer_for_mono_faces.add_point_in_face(kp); }
    else if (m_buffer_for_colored_faces.is_a_face_started())
    { return m_buffer_for_colored_faces.add_point_in_face(kp); }
    return false;
  }

  template<typename KPoint, typename KVector>
  bool add_point_in_face(const KPoint& kp, const KVector& p_normal)
  {
    if (m_buffer_for_mono_faces.is_a_face_started())
    { return m_buffer_for_mono_faces.add_point_in_face(kp, p_normal); }
    else if (m_buffer_for_colored_faces.is_a_face_started())
    { return m_buffer_for_colored_faces.add_point_in_face(kp, p_normal); }
    return false;
  }

  void face_end()
  {
    if (m_buffer_for_mono_faces.is_a_face_started())
    { m_buffer_for_mono_faces.face_end(); }
    else if (m_buffer_for_colored_faces.is_a_face_started())
    { return m_buffer_for_colored_faces.face_end(); }
  }

protected:
  // Shortcuts to simplify function calls.
  template<typename KPoint>
  static Local_point get_local_point(const KPoint& p)
  {
    return internal::Geom_utils<typename CGAL::Kernel_traits<KPoint>::Kernel, Local_kernel>::
      get_local_point(p);
  }
  template<typename KVector>
  static Local_vector get_local_vector(const KVector& v)
  {
    return internal::Geom_utils<typename CGAL::Kernel_traits<KVector>::Kernel, Local_kernel>::
      get_local_vector(v);
  }

  void compile_shaders()
  {
    rendering_program_face.removeAllShaders();
    rendering_program_p_l.removeAllShaders();

    // Create the buffers
    for (unsigned int i=0; i<NB_VBO_BUFFERS; ++i)
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
    { std::cerr<<"Compiling fragmentsource FAILED"<<std::endl; }

    if(!rendering_program_face.addShader(vertex_shader_face))
    { std::cerr<<"adding vertex shader FAILED"<<std::endl; }
    if(!rendering_program_face.addShader(fragment_shader_face))
    { std::cerr<<"adding fragment shader FAILED"<<std::endl; }
    if(!rendering_program_face.link())
    { std::cerr<<"linking Program FAILED"<<std::endl; }
  }

  void initialize_buffers()
  {
    rendering_program_p_l.bind();

    // 1) POINT SHADER

    // 1.1) Mono points
    vao[VAO_MONO_POINTS].bind();

    unsigned int bufn = 0;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_MONO_POINTS].data(),
                           static_cast<int>(arrays[POS_MONO_POINTS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    rendering_program_p_l.disableAttributeArray("color");

    vao[VAO_MONO_POINTS].release();

    // 1.2) Color points
    vao[VAO_COLORED_POINTS].bind();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_COLORED_POINTS].data(),
                           static_cast<int>(arrays[POS_COLORED_POINTS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[bufn].release();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[COLOR_POINTS].data(),
                           static_cast<int>(arrays[COLOR_POINTS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("color");
    rendering_program_p_l.setAttributeBuffer("color",GL_FLOAT,0,3);
    buffers[bufn].release();

    vao[VAO_COLORED_POINTS].release();

    // 2) SEGMENT SHADER

    // 2.1) Mono segments
    vao[VAO_MONO_SEGMENTS].bind();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_MONO_SEGMENTS].data(),
                           static_cast<int>(arrays[POS_MONO_SEGMENTS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    rendering_program_p_l.disableAttributeArray("color");

    vao[VAO_MONO_SEGMENTS].release();

    // 2.1) Color segments
    vao[VAO_COLORED_SEGMENTS].bind();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_COLORED_SEGMENTS].data(),
                           static_cast<int>(arrays[POS_COLORED_SEGMENTS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[COLOR_SEGMENTS].data(),
                           static_cast<int>(arrays[COLOR_SEGMENTS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("color");
    rendering_program_p_l.setAttributeBuffer("color",GL_FLOAT,0,3);
    buffers[bufn].release();

    vao[VAO_COLORED_SEGMENTS].release();

    rendering_program_p_l.release();

    // 3) RAYS SHADER

    // 3.1) Mono rays
    vao[VAO_MONO_RAYS].bind();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_MONO_RAYS].data(),
                           static_cast<int>(arrays[POS_MONO_RAYS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeArray("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    rendering_program_p_l.disableAttributeArray("color");

    vao[VAO_MONO_RAYS].release();

    // 3.2) Color rays

    vao[VAO_COLORED_RAYS].bind();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_COLORED_RAYS].data(),
                           static_cast<int>(arrays[POS_COLORED_RAYS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[COLOR_RAYS].data(),
                           static_cast<int>(arrays[COLOR_RAYS].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("color");
    rendering_program_p_l.setAttributeBuffer("color",GL_FLOAT,0,3);
    buffers[bufn].release();

    vao[VAO_COLORED_RAYS].release();

    rendering_program_p_l.release();

    // 4) LINES SHADER
    // 4.1) Mono lines
    vao[VAO_MONO_LINES].bind();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_MONO_LINES].data(),
                           static_cast<int>(arrays[POS_MONO_LINES].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeArray("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    rendering_program_p_l.disableAttributeArray("color");

    vao[VAO_MONO_LINES].release();

    // 4.2 Color lines

    vao[VAO_COLORED_LINES].bind();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_COLORED_LINES].data(),
                           static_cast<int>(arrays[POS_COLORED_LINES].size()*sizeof(float)));
    rendering_program_p_l.enableAttributeArray("vertex");
    rendering_program_p_l.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[COLOR_LINES].data(),
                           static_cast<int>(arrays[COLOR_LINES].size()*sizeof(float)));
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
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_MONO_FACES].data(),
                           static_cast<int>(arrays[POS_MONO_FACES].size()*sizeof(float)));
    rendering_program_face.enableAttributeArray("vertex");
    rendering_program_face.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    // 5.1.2) normals of the mono faces
    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    if (m_flatShading)
    {
      buffers[bufn].allocate(arrays[FLAT_NORMAL_MONO_FACES].data(),
                                    static_cast<int>(arrays[FLAT_NORMAL_MONO_FACES].size()*
                                                     sizeof(float)));
    }
    else
    {
      buffers[bufn].allocate(arrays[SMOOTH_NORMAL_MONO_FACES].data(),
                                    static_cast<int>(arrays[SMOOTH_NORMAL_MONO_FACES].size()*
                                                       sizeof(float)));
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
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_COLORED_FACES].data(),
                           static_cast<int>(arrays[POS_COLORED_FACES].size()*sizeof(float)));
    rendering_program_face.enableAttributeArray("vertex");
    rendering_program_face.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();

    // 5.2.2) normals of the color faces
    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    if (m_flatShading)
    {
      buffers[bufn].allocate(arrays[FLAT_NORMAL_COLORED_FACES].data(),
                                    static_cast<int>(arrays[FLAT_NORMAL_COLORED_FACES].size()*
                                                     sizeof(float)));
    }
    else
    {
      buffers[bufn].allocate(arrays[SMOOTH_NORMAL_COLORED_FACES].data(),
                                    static_cast<int>(arrays[SMOOTH_NORMAL_COLORED_FACES].size()*
                                                     sizeof(float)));
    }
    rendering_program_face.enableAttributeArray("normal");
    rendering_program_face.setAttributeBuffer("normal",GL_FLOAT,0,3);

    buffers[bufn].release();

    // 5.2.3) colors of the faces
    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[COLOR_FACES].data(),
                           static_cast<int>(arrays[COLOR_FACES].size()*sizeof(float)));
    rendering_program_face.enableAttributeArray("color");
    rendering_program_face.setAttributeBuffer("color",GL_FLOAT,0,3);

    buffers[bufn].release();

    vao[VAO_COLORED_FACES].release();

    rendering_program_face.release();

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
    if (bb==bounding_box()) // Case of "empty" bounding box
    {
      bb=Local_point(CGAL::ORIGIN).bbox();
      bb=bb + Local_point(1,1,1).bbox(); // To avoid a warning from Qglviewer
    }
    else
    { bb=bounding_box(); }

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
  }

  // Returns true if the data structure lies on a plane
  bool is_two_dimensional() {
    return (!is_empty() && (has_zero_x() || has_zero_y() || has_zero_z()));
  }

  virtual void draw()
  {
    glEnable(GL_DEPTH_TEST);
    if(!m_are_buffers_initialized)
    { initialize_buffers(); }

    QColor color;
    attrib_buffers(this);

    if(m_draw_vertices)
    {
      rendering_program_p_l.bind();

      vao[VAO_MONO_POINTS].bind();
      color.setRgbF((double)m_vertices_mono_color.red()/(double)255,
                    (double)m_vertices_mono_color.green()/(double)255,
                    (double)m_vertices_mono_color.blue()/(double)255);
      rendering_program_p_l.setAttributeValue("color",color);
      rendering_program_p_l.setUniformValue("point_size", GLfloat(m_size_points));
      glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(arrays[POS_MONO_POINTS].size()/3));
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
      glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(arrays[POS_COLORED_POINTS].size()/3));
      vao[VAO_COLORED_POINTS].release();

      rendering_program_p_l.release();
    }

    if(m_draw_edges)
    {
      rendering_program_p_l.bind();

      vao[VAO_MONO_SEGMENTS].bind();
      color.setRgbF((double)m_edges_mono_color.red()/(double)255,
                    (double)m_edges_mono_color.green()/(double)255,
                    (double)m_edges_mono_color.blue()/(double)255);
      rendering_program_p_l.setAttributeValue("color",color);
      glLineWidth(m_size_edges);
      glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(arrays[POS_MONO_SEGMENTS].size()/3));
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
      glLineWidth(m_size_edges);
      glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(arrays[POS_COLORED_SEGMENTS].size()/3));
      vao[VAO_COLORED_SEGMENTS].release();

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
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(arrays[POS_MONO_RAYS].size()/3));
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
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(arrays[POS_COLORED_RAYS].size()/3));
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
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(arrays[POS_MONO_LINES].size()/3));
        vao[VAO_MONO_LINES].release();

        rendering_program_p_l.release();

        vao[VAO_COLORED_LINES].bind();
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
        glLineWidth(m_size_lines);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(arrays[POS_COLORED_LINES].size()/3));
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

      vao[VAO_MONO_FACES].bind();
      color.setRgbF((double)m_faces_mono_color.red()/(double)255,
                    (double)m_faces_mono_color.green()/(double)255,
                    (double)m_faces_mono_color.blue()/(double)255);
      rendering_program_face.setAttributeValue("color",color);
      glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(arrays[POS_MONO_FACES].size()/3));
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
      glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(arrays[POS_COLORED_FACES].size()/3));
      vao[VAO_COLORED_FACES].release();

      if (is_two_dimensional())
        glPolygonOffset(offset_factor, offset_units);

      rendering_program_face.release();
    }

    if (is_two_dimensional())
    {
      camera()->setType(CGAL::qglviewer::Camera::ORTHOGRAPHIC);
      //      Camera Constraint:
      constraint.setRotationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::AXIS);
      constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);

      double cx=0., cy=0., cz=0.;
      if (has_zero_x())      { cx=1.; }
      else if (has_zero_y()) { cy=1.; }
      else                   { cz=1.; }
    
      camera()->setViewDirection(CGAL::qglviewer::Vec(-cx,-cy,-cz));
      constraint.setRotationConstraintDirection(CGAL::qglviewer::Vec(cx, cy, cz));
      camera()->frame()->setConstraint(&constraint);
    }

    if (m_draw_text)
    {
      glDisable(GL_LIGHTING);
      for (std::size_t i=0; i<m_texts.size(); ++i)
      {
        CGAL::qglviewer::Vec screenPos=camera()->projectedCoordinatesOf
          (CGAL::qglviewer::Vec(std::get<0>(m_texts[i]).x(),
                                std::get<0>(m_texts[i]).y(),
                                std::get<0>(m_texts[i]).z()));
        
        drawText((int)screenPos[0], (int)screenPos[1], std::get<1>(m_texts[i]));
      }
      glEnable(GL_LIGHTING);
    }
  }

  virtual void redraw()
  {
    initialize_buffers();
    update();
  }

  virtual void init()
  {
    // Restore previous viewer state.
    restoreStateFromFile();
    initializeOpenGLFunctions();

    // Define 'Control+Q' as the new exit shortcut (default was 'Escape')
    setShortcut(qglviewer::EXIT_VIEWER, ::Qt::CTRL+::Qt::Key_Q);

    // Add custom key description (see keyPressEvent).
    setKeyDescription(::Qt::Key_E, "Toggles edges display");
    setKeyDescription(::Qt::Key_F, "Toggles faces display");
    setKeyDescription(::Qt::Key_G, "Switch between flat/Gouraud shading display");
    setKeyDescription(::Qt::Key_M, "Toggles mono color");
    setKeyDescription(::Qt::Key_N, "Inverse direction of normals");
    setKeyDescription(::Qt::Key_T, "Toggles text display");
    setKeyDescription(::Qt::Key_U, "Move camera direction upside down");
    setKeyDescription(::Qt::Key_V, "Toggles vertices display");
    setKeyDescription(::Qt::Key_Plus, "Increase size of edges");
    setKeyDescription(::Qt::Key_Minus, "Decrease size of edges");
    setKeyDescription(::Qt::Key_Plus+::Qt::ControlModifier, "Increase size of vertices");
    setKeyDescription(::Qt::Key_Minus+::Qt::ControlModifier, "Decrease size of vertices");
    setKeyDescription(::Qt::Key_PageDown, "Increase light (all colors, use shift/alt/ctrl for one rgb component)");
    setKeyDescription(::Qt::Key_PageUp, "Decrease light (all colors, use shift/alt/ctrl for one rgb component)");

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
    if (bb==bounding_box()) // Case of "empty" bounding box
    {
      bb=Local_point(CGAL::ORIGIN).bbox();
      bb=bb + Local_point(1,1,1).bbox(); // To avoid a warning from Qglviewer
    }
    else
    { bb=bounding_box(); }
    this->camera()->setSceneBoundingBox(CGAL::qglviewer::Vec(bb.xmin(),
                                                       bb.ymin(),
                                                       bb.zmin()),
                                        CGAL::qglviewer::Vec(bb.xmax(),
                                                       bb.ymax(),
                                                       bb.zmax()));

    this->showEntireScene();
  }

  void negate_all_normals()
  {
    m_buffer_for_mono_faces.negate_normals();
    m_buffer_for_colored_faces.negate_normals();
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    const ::Qt::KeyboardModifiers modifiers = e->modifiers();

    if ((e->key()==::Qt::Key_E) && (modifiers==::Qt::NoButton))
    {
      m_draw_edges=!m_draw_edges;
      displayMessage(QString("Draw edges=%1.").arg(m_draw_edges?"true":"false"));
      update();
    }else if ((e->key()==::Qt::Key_F) && (modifiers==::Qt::NoButton))
    {
      m_draw_faces=!m_draw_faces;
      displayMessage(QString("Draw faces=%1.").arg(m_draw_faces?"true":"false"));
      update();
    }
    else if ((e->key()==::Qt::Key_G) && (modifiers==::Qt::NoButton))
    {
      m_flatShading=!m_flatShading;
      if (m_flatShading)
        displayMessage("Flat shading.");
      else
        displayMessage("Gouraud shading.");
      redraw();
    }
    else if ((e->key()==::Qt::Key_M) && (modifiers==::Qt::NoButton))
    {
      m_use_mono_color=!m_use_mono_color;
      displayMessage(QString("Mono color=%1.").arg(m_use_mono_color?"true":"false"));
      update();
    }
    else if ((e->key()==::Qt::Key_N) && (modifiers==::Qt::NoButton))
    {
      m_inverse_normal=!m_inverse_normal;
      displayMessage(QString("Inverse normal=%1.").arg(m_inverse_normal?"true":"false"));
      negate_all_normals();
      redraw();
    }
    else if ((e->key()==::Qt::Key_U) && (modifiers==::Qt::NoButton))
    {
      if (is_two_dimensional())
      {
        displayMessage(QString("Move camera direction upside down."));
        /* CGAL::qglviewer::Vec cur=camera()->viewDirection();
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
    else if ((e->key()==::Qt::Key_T) && (modifiers==::Qt::NoButton))
    {
      m_draw_text=!m_draw_text;
      displayMessage(QString("Draw text=%1.").arg(m_draw_text?"true":"false"));
      update();
    }
    else if ((e->key()==::Qt::Key_V) && (modifiers==::Qt::NoButton))
    {
      m_draw_vertices=!m_draw_vertices;
      displayMessage(QString("Draw vertices=%1.").arg(m_draw_vertices?"true":"false"));
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
    else
      CGAL::QGLViewer::keyPressEvent(e);
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

protected:
  bool m_draw_vertices;
  bool m_draw_edges;
  bool m_draw_rays;
  bool m_draw_lines;
  bool m_draw_faces;
  bool m_flatShading;
  bool m_use_mono_color;
  bool m_inverse_normal;
  bool m_draw_text;
  
  double m_size_points;
  double m_size_edges;
  double m_size_rays;
  double m_size_lines;

  CGAL::Color m_vertices_mono_color;
  CGAL::Color m_edges_mono_color;
  CGAL::Color m_rays_mono_color;
  CGAL::Color m_lines_mono_color;
  CGAL::Color m_faces_mono_color;
  QVector4D   m_ambient_color;

  bool m_are_buffers_initialized;
  CGAL::Bbox_3 m_bounding_box;

  // CGAL::qglviewer::LocalConstraint constraint;
  CGAL::qglviewer::WorldConstraint constraint;

  // The following enum gives the indices of different elements of arrays vectors.
  enum
  {
    BEGIN_POS=0,
    POS_MONO_POINTS=BEGIN_POS,
    POS_COLORED_POINTS,
    POS_MONO_SEGMENTS,
    POS_COLORED_SEGMENTS,
    POS_MONO_RAYS,
    POS_COLORED_RAYS,
    POS_MONO_LINES,
    POS_COLORED_LINES,
    POS_MONO_FACES,
    POS_COLORED_FACES,
    END_POS,
    BEGIN_COLOR=END_POS,
    COLOR_POINTS=BEGIN_COLOR,
    COLOR_SEGMENTS,
    COLOR_RAYS,
    COLOR_LINES,
    COLOR_FACES,
    END_COLOR,
    BEGIN_NORMAL=END_COLOR,
    SMOOTH_NORMAL_MONO_FACES=BEGIN_NORMAL,
    FLAT_NORMAL_MONO_FACES,
    SMOOTH_NORMAL_COLORED_FACES,
    FLAT_NORMAL_COLORED_FACES,
    END_NORMAL,
    LAST_INDEX=END_NORMAL
  };
  std::vector<float> arrays[LAST_INDEX];

  Buffer_for_vao<float> m_buffer_for_mono_points;
  Buffer_for_vao<float> m_buffer_for_colored_points;
  Buffer_for_vao<float> m_buffer_for_mono_segments;
  Buffer_for_vao<float> m_buffer_for_colored_segments;
  Buffer_for_vao<float> m_buffer_for_mono_rays;
  Buffer_for_vao<float> m_buffer_for_colored_rays;
  Buffer_for_vao<float> m_buffer_for_mono_lines;
  Buffer_for_vao<float> m_buffer_for_colored_lines;
  Buffer_for_vao<float> m_buffer_for_mono_faces;
  Buffer_for_vao<float> m_buffer_for_colored_faces;

  static const unsigned int NB_VBO_BUFFERS=(END_POS-BEGIN_POS)+
    (END_COLOR-BEGIN_COLOR)+2; // +2 for 2 vectors of normals

  QGLBuffer buffers[NB_VBO_BUFFERS];

  // The following enum gives the indices of the differents vao.
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
      NB_VAO_BUFFERS
    };
  QOpenGLVertexArrayObject vao[NB_VAO_BUFFERS];

  QOpenGLShaderProgram rendering_program_face;
  QOpenGLShaderProgram rendering_program_p_l;

  std::vector<std::tuple<Local_point, QString> > m_texts;
};

} // End namespace CGAL

#else // CGAL_USE_BASIC_VIEWER

namespace CGAL
{

  template<class T>
  void draw(const T&, const char* ="", bool=false)
  { 
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_BASIC_VIEWER_QT_H
