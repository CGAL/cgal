// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_BASIC_VIEWER_H
#define CGAL_BASIC_VIEWER_H

#include <QApplication>
#include <QKeyEvent>

#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <QOpenGLFunctions_2_1>
#include <QOpenGLVertexArrayObject>
#include <QGLBuffer>
#include <QOpenGLShaderProgram>

#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Qt/CreateOpenGLContext.h>

#include <vector>
#include <cstdlib>
#include <queue>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3  Local_point;
typedef Local_kernel::Vector_3 Local_vector;

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
    "void main(void)\n"
    "{\n"
    "   fP = mv_matrix * vertex; \n"
    "   fN = mat3(mv_matrix)* normal; \n"
    "   fColor = vec4(color, 1.0); \n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
  };

const char fragment_source_color[] =
  {
    "#version 120 \n"
    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "varying highp vec4 fColor; \n"
    "uniform vec4 light_pos;  \n"
    "uniform vec4 light_diff; \n"
    "uniform vec4 light_spec; \n"
    "uniform vec4 light_amb;  \n"
    "uniform float spec_power ; \n"
    
    "void main(void) { \n"
    
    "   vec3 L = light_pos.xyz - fP.xyz; \n"
    "   vec3 V = -fP.xyz; \n"
    
    "   vec3 N = normalize(fN); \n"
    "   L = normalize(L); \n"
    "   V = normalize(V); \n"
    
    "   vec3 R = reflect(-L, N); \n"
    "   vec4 diffuse = max(dot(N,L), 0.0) * light_diff * fColor; \n"
    "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"
    
    "gl_FragColor = light_amb*fColor + diffuse  ; \n"
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
    "void main(void)\n"
    "{\n"
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

namespace internal {
  template <class Point, class Vector>
  void newell_single_step_3(const Point& p, const Point& q, Vector& n)
  {
    // Compute normal of the face by using Newell's method: for each edge PQ
    // Nx += (Py - Qy) * (Pz + Qz);
    // Ny += (Pz - Qz) * (Px + Qx);
    // Nz += (Px - Qx) * (Py + Qy);
    n = Vector(n.x()+((p.y()-q.y())*(p.z()+q.z())),
               n.y()+((p.z()-q.z())*(p.x()+q.x())),
               n.z()+((p.x()-q.x())*(p.y()+q.y())));
    }
} // End namespace internal

template <class K>
typename K::Vector_3 compute_normal_of_face(const std::vector<typename K::Point_3>& points)
{
  typename K::Vector_3 normal(CGAL::NULL_VECTOR);
  unsigned int nb = 0;
  for (std::size_t i=0; i<points.size(); ++i)
  {
    internal::newell_single_step_3(points[i], points[(i+1)%points.size()], normal);
    ++nb;
  }

  assert(nb>0);
  return (typename K::Construct_scaled_vector_3()(normal, 1.0/nb));
}

class Basic_viewer : public QGLViewer, public QOpenGLFunctions_2_1
{
  struct Vertex_info
  {
    Local_vector v;
  };

  struct Face_info
  {
    bool exist_edge[3];
    bool is_external;
    bool is_process;
  };

  typedef CGAL::Triangulation_2_projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
  typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, P_traits> Vb;

  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, P_traits> Fb1;

  typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>    Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                   TDS;
  // typedef CGAL::No_intersection_tag                                     Itag;
  typedef CGAL::Exact_predicates_tag                                    Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS,
                                                     Itag>              CDT;

public:
  // Constructor/Destructor
  Basic_viewer(const char* title="") :
    QGLViewer(CGAL::Qt::createOpenGLContext()),
    m_draw_vertices(true),
    m_draw_edges(true),
    m_draw_faces(true),
    m_flatShading(true),
    m_use_mono_color(false),
    m_inverse_normal(false),
    m_empty(true),
    m_size_points(7.),
    m_size_edges(3.1),    
    m_vertices_mono_color(51, 51, 178),
    m_edges_mono_color(51, 51, 148),
    m_faces_mono_color(180, 125, 200),
    m_ambient_color(0.6f, 0.5f, 0.5f, 0.5f),
    m_are_buffers_initialized(false),
    m_face_started(false)
  {
    if (title[0]==0)
      setWindowTitle("CGAL Basic Viewer");
    else
      setWindowTitle(title);

    resize(500, 450);

    if ( is_empty() )
    {
      bb=Local_point(CGAL::ORIGIN).bbox();
      bb=bb + Local_point(1,1,1).bbox(); // To avoid a warning from Qglviewer
    }
  }

  ~Basic_viewer()
  {
    for (int i=0; i<NB_VBO_BUFFERS; ++i)
      buffers[i].destroy();

    for (int i=0; i<NB_VAO_BUFFERS; ++i)
      vao[i].destroy();
  }

  void clear()
  {
    for (unsigned int i=0; i<LAST_INDEX; ++i)
    { arrays[i].clear(); }    
  }

  bool is_empty() const
  { return m_empty; }
  
  void add_point(const Local_point& p, std::vector<float>& point_vector)
  {
    point_vector.push_back(p.x());
    point_vector.push_back(p.y());
    point_vector.push_back(p.z());

    if (is_empty())
    { bb=p.bbox(); m_empty=false; }
    else
    { bb=bb+p.bbox(); }
  }

  void add_color(const CGAL::Color& acolor, std::vector<float>& color_vector)
  {
    color_vector.push_back((float)acolor.red()/(float)255);
    color_vector.push_back((float)acolor.green()/(float)255);
    color_vector.push_back((float)acolor.blue()/(float)255);
  }

  void add_normal(const Local_vector& n, std::vector<float>& normal_vector)
  {
    normal_vector.push_back(n.x());
    normal_vector.push_back(n.y());
    normal_vector.push_back(n.z());
  }
  
  void add_mono_point(const Local_point& p)
  { add_point(p, arrays[POS_MONO_POINTS]); }

  void add_colored_point(const Local_point& p, const CGAL::Color& acolor)
  {
    add_point(p, arrays[POS_COLORED_POINTS]);
    add_color(acolor, arrays[COLOR_POINTS]);
  }
  
  void add_mono_segment(const Local_point& p1, const Local_point& p2)
  {
    add_point(p1, arrays[POS_MONO_SEGMENTS]);
    add_point(p2, arrays[POS_MONO_SEGMENTS]);
  }
  
  void add_colored_segment(const Local_point& p1, const Local_point& p2,
                           const CGAL::Color& acolor)
  {
    add_point(p1, arrays[POS_COLORED_SEGMENTS]);
    add_point(p2, arrays[POS_COLORED_SEGMENTS]);
    add_color(acolor, arrays[COLOR_SEGMENTS]);
    add_color(acolor, arrays[COLOR_SEGMENTS]);
  }

  void face_begin()
  {
    if (m_face_started)
    {
      std::cerr<<"You cannot start a new face before to finish the previous one."<<std::endl;
      return;
    }
    
    m_face_started=true;
  }

  void mono_face_begin()
  {
    m_started_face_is_colored=false;
    face_begin();
  }

  /// Start a new face, with a given color.
  void colored_face_begin(const CGAL::Color& acolor)
  {
    color_of_face=acolor;
    m_started_face_is_colored=true;
    face_begin();
  }

  /// Add a point at the end of the current face
  /// With this method, it is not possible to use the Gourod shading.
  /// @param p the point to add
  bool add_point_in_face(const Local_point& p)
  {
    if (!m_face_started) return false;
    
    if (points_of_face.empty() || points_of_face.back()!=p)
    {
      points_of_face.push_back(p);
      return true;
    }
    return false;
  }
  
  /// Add a point at the end of the current face
  /// @param p the point to add
  /// @p_normal the vertex normal in this point (for Gourod shading)
  void add_point_in_face(const Local_point& p, const Local_vector& p_normal)
  {
    if (!m_face_started) return;
    
    if (add_point_in_face(p))
    {
      vertex_normals_for_face.push_back(p_normal);
    }
  }

  /// End the face: compute the triangulation.
  void face_end()
  {
    if (points_of_face.size()<3)
    {
      std::cout<<"PB: you try to triangulate a face with "<<points_of_face.size()<<" vertices."
               <<std::endl;
      
      m_face_started=false;
      points_of_face.clear();
      vertex_normals_for_face.clear();

      return;
    }
    
    Local_vector normal=compute_normal_of_face<Local_kernel>(points_of_face);

    if (points_of_face.size()==3) // Triangle: no need to triangulate
    {
      for (int i=0; i<3; ++i)
      {
        // The point
        add_point(points_of_face[i], arrays[m_started_face_is_colored?
                                            POS_COLORED_FACES:
                                            POS_MONO_FACES]);

        // Its color
        if (m_started_face_is_colored)
        { add_color(color_of_face, arrays[COLOR_FACES]); }

        // Its flat normal
        add_normal(normal, arrays[m_started_face_is_colored?
                                            FLAT_NORMAL_COLORED_FACES:
                                            FLAT_NORMAL_MONO_FACES]);

        // Its smoth normal (if given by the user)
        if (vertex_normals_for_face.size()==3)
        { // Here we have 3 vertex normals; we can use Gourod
          add_normal(vertex_normals_for_face[i], arrays[m_started_face_is_colored?
                                                        SMOOTH_NORMAL_COLORED_FACES:
                                                        SMOOTH_NORMAL_MONO_FACES]);
        }
        else
        { // Here user does not provide all vertex normals: we use face normal istead
          // and thus we will not be able to use Gourod
         add_normal(normal, arrays[m_started_face_is_colored?
                                   SMOOTH_NORMAL_COLORED_FACES:
                                   SMOOTH_NORMAL_MONO_FACES]);
        }
      }
    }
    // TODO CASE OF 4 POINTS ? PB HOW TO FIND (EASILY) THE TWO POINTS TO LINK ?
    // else if (points_of_face.size()==4)
    else
    { // More than 3 points: we triangulate
      try
      {
        P_traits cdt_traits(normal);
        CDT cdt(cdt_traits);

        bool with_vertex_normal=(vertex_normals_for_face.size()==points_of_face.size());

        // (1) We insert all the edges as contraint in the CDT.
        typename CDT::Vertex_handle previous=NULL, first=NULL;
        for (int i=0; i<points_of_face.size(); ++i)
        {
          typename CDT::Vertex_handle vh = cdt.insert(points_of_face[i]);
          if(first==NULL)
          { first=vh; }

          if (with_vertex_normal)
          { vh->info().v=vertex_normals_for_face[i]; }
          else
          { vh->info().v=normal; }

          if(previous!=NULL && previous!=vh)
          { cdt.insert_constraint(previous, vh); }
          previous=vh;
        }

        if (previous!=NULL && previous!=first)
          cdt.insert_constraint(previous, first);

        // (2) We mark all external triangles
        // (2.1) We initialize is_external and is_process values 
        for(typename CDT::All_faces_iterator fit = cdt.all_faces_begin(),
              fitend = cdt.all_faces_end(); fit!=fitend; ++fit)
        {
          fit->info().is_external = true;
          fit->info().is_process = false;
        }
        // (2.2) We check if the facet is external or internal
        std::queue<typename CDT::Face_handle> face_queue;
        typename CDT::Face_handle face_internal = NULL;
        if (cdt.infinite_vertex()->face()!=NULL)
          face_queue.push(cdt.infinite_vertex()->face());
        while(! face_queue.empty() )
        {
          typename CDT::Face_handle fh = face_queue.front();
          face_queue.pop();
          if(!fh->info().is_process)
          {
            fh->info().is_process = true;
            for(int i=0; i<3; ++i)
            {
              if(!cdt.is_constrained(std::make_pair(fh, i)))
              {
                if (fh->neighbor(i)!=NULL)
                  face_queue.push(fh->neighbor(i));
              }
              else if (face_internal==NULL)
              {
                face_internal = fh->neighbor(i);
              }
            }
          }
        }

        if ( face_internal!=NULL )
          face_queue.push(face_internal);
        
        while(! face_queue.empty() )
        {
          typename CDT::Face_handle fh = face_queue.front();
          face_queue.pop();
          if(!fh->info().is_process)
          {
            fh->info().is_process = true;
            fh->info().is_external = false;
            for(int i=0; i<3; ++i)
            {
              if(!cdt.is_constrained(std::make_pair(fh, i)))
              {
                if (fh->neighbor(i)!=NULL)
                  face_queue.push(fh->neighbor(i));
              }
            }
          }
        }

        // (3) Now we iterates on the internal faces to add the vertices to the
        //     positions and the normals to the appropriate vectors
        for(typename CDT::Finite_faces_iterator ffit=cdt.finite_faces_begin(),
              ffitend = cdt.finite_faces_end(); ffit!=ffitend; ++ffit)
        {
          if(!ffit->info().is_external)
          {
            for(int i=0; i<3; ++i)
            {
              // The point
              add_point(ffit->vertex(i)->point(), arrays[m_started_face_is_colored?
                                                         POS_COLORED_FACES:
                                                         POS_MONO_FACES]);

              // Its color
              if (m_started_face_is_colored)
              { add_color(color_of_face, arrays[COLOR_FACES]); }
              
              // Its flat normal
              add_normal(normal, arrays[m_started_face_is_colored?
                                        FLAT_NORMAL_COLORED_FACES:
                                        FLAT_NORMAL_MONO_FACES]);

              // Its smoth normal (if given by the user)
              add_normal(ffit->vertex(i)->info().v, arrays[m_started_face_is_colored?
                                                           SMOOTH_NORMAL_COLORED_FACES:
                                                           SMOOTH_NORMAL_MONO_FACES]);
            }
          }
        }
      }
      catch(...)
      { // Triangulation crash: the face is not filled
        std::cout<<"Catch: face not filled."<<std::endl;
      }
    }

    m_face_started=false;
    points_of_face.clear();
    vertex_normals_for_face.clear();
  }
  
protected:
  
  void compile_shaders()
  {
    rendering_program_face.removeAllShaders();
    rendering_program_p_l.removeAllShaders();
    
    // Create the buffers
    for (int i=0; i<NB_VBO_BUFFERS; ++i)
    {
      if(!buffers[i].isCreated() && !buffers[i].create())
      { std::cerr<<"VBO Creation number "<<i<<" FAILED"<<std::endl; }
    }
    
    for (int i=0; i<NB_VAO_BUFFERS; ++i)
    {
      if(!vao[i].isCreated() && !vao[i].create())
      { std::cerr<<"VAO Creation number "<<i<<" FAILED"<<std::endl; }
    }
    
    // Vertices and segments shader
    QOpenGLShader *vertex_shader_p_l = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_p_l->compileSourceCode(vertex_source_p_l))
    { std::cerr<<"Compiling vertex source FAILED"<<std::endl; }

    QOpenGLShader *fragment_shader_p_l= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader_p_l->compileSourceCode(fragment_source_p_l))
    { std::cerr<<"Compiling fragmentsource FAILED"<<std::endl; }

    if(!rendering_program_p_l.addShader(vertex_shader_p_l))
    { std::cerr<<"adding vertex shader FAILED"<<std::endl; }
    if(!rendering_program_p_l.addShader(fragment_shader_p_l))
    { std::cerr<<"adding fragment shader FAILED"<<std::endl; }
    if(!rendering_program_p_l.link())
    { std::cerr<<"linking Program FAILED"<<std::endl; }

    // Faces shader
    QOpenGLShader *vertex_shader_face = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_face->compileSourceCode(vertex_source_color))
    { std::cerr<<"Compiling vertex source FAILED"<<std::endl; }

    QOpenGLShader *fragment_shader_face= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader_face->compileSourceCode(fragment_source_color))
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
    
    int bufn = 0;
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

    // 1.2) Color segments
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
    
    // 3) FACE SHADER
    rendering_program_face.bind();

    // 3.1) Mono faces
    vao[VAO_MONO_FACES].bind();

    // 3.1.1) points of the mono faces
    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_MONO_FACES].data(),
                           static_cast<int>(arrays[POS_MONO_FACES].size()*sizeof(float)));
    rendering_program_face.enableAttributeArray("vertex");
    rendering_program_face.setAttributeBuffer("vertex",GL_FLOAT,0,3);

    buffers[bufn].release();
    
    // 3.1.2) normals of the mono faces
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

    // 3.1.3) color of the mono faces
    rendering_program_face.disableAttributeArray("color");
    vao[VAO_MONO_FACES].release();

    // 3.2) Color faces
    vao[VAO_COLORED_FACES].bind();
  
    // 3.2.1) points of the color faces    
    ++bufn;
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(arrays[POS_COLORED_FACES].data(),
                           static_cast<int>(arrays[POS_COLORED_FACES].size()*sizeof(float)));
    rendering_program_face.enableAttributeArray("vertex");
    rendering_program_face.setAttributeBuffer("vertex",GL_FLOAT,0,3);
      
    buffers[bufn].release();
    
    // 3.2.2) normals of the color faces
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
    
    // 3.2.3) colors of the faces
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

  void attrib_buffers(QGLViewer* viewer)
  {
    QMatrix4x4 mvpMatrix;
    QMatrix4x4 mvMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
      mvpMatrix.data()[i] = (float)mat[i];
    }
    viewer->camera()->getModelViewMatrix(mat);
    for(int i=0; i < 16; i++)
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

    QVector4D position((bb.xmax()-bb.xmin())/2, (bb.ymax()-bb.ymin())/2,bb.zmax(), 0.0 );
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

  virtual void draw()
  {
    glEnable(GL_DEPTH_TEST);
    if(!m_are_buffers_initialized)
      initialize_buffers();

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
      ::glPointSize(m_size_points);
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
      ::glPointSize(m_size_points);
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
      ::glLineWidth(m_size_edges);
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
      ::glLineWidth(m_size_edges);
      glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(arrays[POS_COLORED_SEGMENTS].size()/3));
      vao[VAO_COLORED_SEGMENTS].release();

      rendering_program_p_l.release();
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

      rendering_program_face.release();
      }
  }

  virtual void init()
  {
    // Restore previous viewer state.
    restoreStateFromFile();
    initializeOpenGLFunctions();

    // Define 'Control+Q' as the new exit shortcut (default was 'Escape')
    setShortcut(EXIT_VIEWER, Qt::CTRL+Qt::Key_Q);

    // Add custom key description (see keyPressEvent).
    setKeyDescription(Qt::Key_E, "Toggles edges display");
    setKeyDescription(Qt::Key_F, "Toggles faces display");
    setKeyDescription(Qt::Key_G, "Switch between flat/Gouraud shading display");
    setKeyDescription(Qt::Key_M, "Toggles mono color for all faces");
    setKeyDescription(Qt::Key_N, "Inverse direction of normals");
    setKeyDescription(Qt::Key_V, "Toggles vertices display");
    setKeyDescription(Qt::Key_Plus, "Increase size of edges");
    setKeyDescription(Qt::Key_Minus, "Decrease size of edges");
    setKeyDescription(Qt::Key_Plus+Qt::ControlModifier, "Increase size of vertices");
    setKeyDescription(Qt::Key_Minus+Qt::ControlModifier, "Decrease size of vertices");
    setKeyDescription(Qt::Key_PageDown, "Increase light (all colors, use shift/alt/ctrl for one rgb component)");
    setKeyDescription(Qt::Key_PageUp, "Decrease light (all colors, use shift/alt/ctrl for one rgb component)");

    // Light default parameters
    ::glLineWidth(m_size_edges);
    ::glPointSize(m_size_points);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(1.f,1.f);
    ::glClearColor(1.0f,1.0f,1.0f,0.0f);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    ::glEnable(GL_LIGHTING);
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    ::glShadeModel(GL_FLAT);
    ::glDisable(GL_BLEND);
    ::glEnable(GL_LINE_SMOOTH);
    ::glDisable(GL_POLYGON_SMOOTH_HINT);
    ::glBlendFunc(GL_ONE, GL_ZERO);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);

    compile_shaders();

    this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(),
                                                       bb.ymin(),
                                                       bb.zmin()),
                                        qglviewer::Vec(bb.xmax(),
                                                       bb.ymax(),
                                                       bb.zmax()));

    this->showEntireScene();
  }

  void negate_all_normals()
  {
    for (unsigned int k=BEGIN_NORMAL; k<END_NORMAL; ++k)
    {
      for (std::size_t i=0; i<arrays[k].size(); ++i)
      { arrays[k][i]=-arrays[k][i]; }
    }
  }
  
  virtual void keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    if ((e->key()==Qt::Key_E) && (modifiers==Qt::NoButton))
    {
      m_draw_edges=!m_draw_edges;
      displayMessage(QString("Draw edges=%1.").arg(m_draw_edges?"true":"false"));
      update();
    }
    else if ((e->key()==Qt::Key_F) && (modifiers==Qt::NoButton))
    {
      m_draw_faces=!m_draw_faces;
      displayMessage(QString("Draw faces=%1.").arg(m_draw_faces?"true":"false"));
      update();
    }
    else if ((e->key()==Qt::Key_G) && (modifiers==Qt::NoButton))
    {
      m_flatShading=!m_flatShading;
      if (m_flatShading)
        displayMessage("Flat shading.");
      else
        displayMessage("Gouraud shading.");
      initialize_buffers();
      update();
    }
    else if ((e->key()==Qt::Key_M) && (modifiers==Qt::NoButton))
    {
      m_use_mono_color=!m_use_mono_color;
      displayMessage(QString("Mono color=%1.").arg(m_use_mono_color?"true":"false"));
      update();
    }
    else if ((e->key()==Qt::Key_N) && (modifiers==Qt::NoButton))
    {
      m_inverse_normal=!m_inverse_normal;
      displayMessage(QString("Inverse normal=%1.").arg(m_inverse_normal?"true":"false"));
      negate_all_normals();
      initialize_buffers();
      update();
    }
    else if ((e->key()==Qt::Key_V) && (modifiers==Qt::NoButton))
    {
      m_draw_vertices=!m_draw_vertices;
      displayMessage(QString("Draw vertices=%1.").arg(m_draw_vertices?"true":"false"));
      update();
    }
    else if ((e->key()==Qt::Key_Plus) && (!modifiers.testFlag(Qt::ControlModifier))) // No ctrl
    {
      m_size_edges+=.5;
      displayMessage(QString("Size of edges=%1.").arg(m_size_edges));
      update();
    }
    else if ((e->key()==Qt::Key_Minus) && (!modifiers.testFlag(Qt::ControlModifier))) // No ctrl
    {
      if (m_size_edges>.5) m_size_edges-=.5;
      displayMessage(QString("Size of edges=%1.").arg(m_size_edges));
      update();
    }
    else if ((e->key()==Qt::Key_Plus) && (modifiers.testFlag(Qt::ControlModifier)))
    {
      m_size_points+=.5;
      displayMessage(QString("Size of points=%1.").arg(m_size_points));
      update();
    }
    else if ((e->key()==Qt::Key_Minus) && (modifiers.testFlag(Qt::ControlModifier)))
    {
      if (m_size_points>.5) m_size_points-=.5;
      displayMessage(QString("Size of points=%1.").arg(m_size_points));
      update();
    }
    else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton))
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
    else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::NoButton))
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
    else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::ShiftModifier))
    {
      m_ambient_color.setX(m_ambient_color.x()+.1);
      if (m_ambient_color.x()>1.) m_ambient_color.setX(1.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
      update();
    }
    else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::AltModifier))
    {
      m_ambient_color.setY(m_ambient_color.y()+.1);
      if (m_ambient_color.y()>1.) m_ambient_color.setY(1.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
      update();
    }
    else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::ControlModifier))
    {
      m_ambient_color.setZ(m_ambient_color.z()+.1);
      if (m_ambient_color.z()>1.) m_ambient_color.setZ(1.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
      update();
    }
    else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::ShiftModifier))
    {
      m_ambient_color.setX(m_ambient_color.x()-.1);
      if (m_ambient_color.x()<0.) m_ambient_color.setX(0.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
      update();
    }
    else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::AltModifier))
    {
      m_ambient_color.setY(m_ambient_color.y()-.1);
      if (m_ambient_color.y()<0.) m_ambient_color.setY(0.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
      update();
    }
    else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::ControlModifier))
    {
      m_ambient_color.setZ(m_ambient_color.z()-.1);
      if (m_ambient_color.z()<0.) m_ambient_color.setZ(0.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(m_ambient_color.x()).arg(m_ambient_color.y()).arg(m_ambient_color.z()));
      update();
    }
    else
      QGLViewer::keyPressEvent(e);
  }

  virtual QString helpString() const
  {
    QString text("<h2>C G A L   B a s i c   V i e w e r</h2>");
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

private:
  bool m_draw_vertices;
  bool m_draw_edges;
  bool m_draw_faces;
  bool m_flatShading;
  bool m_use_mono_color;
  bool m_inverse_normal;
  bool m_empty;
  
  double m_size_points;
  double m_size_edges;

  CGAL::Color m_vertices_mono_color;
  CGAL::Color m_edges_mono_color;
  CGAL::Color m_faces_mono_color;
  QVector4D   m_ambient_color;

  bool m_are_buffers_initialized;
  CGAL::Bbox_3 bb;
  
  // The following enum gives the indices of different elements of arrays vectors.
  enum
  {
    BEGIN_POS=0,
    POS_MONO_POINTS=BEGIN_POS,
    POS_COLORED_POINTS,
    POS_MONO_SEGMENTS,
    POS_COLORED_SEGMENTS,
    POS_MONO_FACES,
    POS_COLORED_FACES,
    END_POS,
    BEGIN_COLOR=END_POS,
    COLOR_POINTS=BEGIN_COLOR,
    COLOR_SEGMENTS,
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

  static const unsigned int NB_VBO_BUFFERS=(END_POS-BEGIN_POS)+
    (END_COLOR-BEGIN_COLOR)+2; // +2 for 2 vectors of normals

  QGLBuffer buffers[NB_VBO_BUFFERS];

  // The following enum gives the indices of differents vao.
  enum 
    { VAO_MONO_POINTS=0,
      VAO_COLORED_POINTS,
      VAO_MONO_SEGMENTS,
      VAO_COLORED_SEGMENTS,
      VAO_MONO_FACES,
      VAO_COLORED_FACES,
      NB_VAO_BUFFERS
    };
  QOpenGLVertexArrayObject vao[NB_VAO_BUFFERS];

  QOpenGLShaderProgram rendering_program_face;
  QOpenGLShaderProgram rendering_program_p_l;

  // Local variables, used when we started a new face.
  bool m_face_started;
  bool m_started_face_is_colored;
  std::vector<Local_point> points_of_face;
  std::vector<Local_vector> vertex_normals_for_face;
  CGAL::Color color_of_face;
};

#endif // CGAL_BASIC_VIEWER_H
