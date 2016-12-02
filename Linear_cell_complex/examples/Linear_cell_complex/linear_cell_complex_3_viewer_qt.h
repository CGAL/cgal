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

#ifndef CGAL_LCC_3_VIEWER_QT_H
#define CGAL_LCC_3_VIEWER_QT_H

#include <QApplication>
#include <QKeyEvent>

#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <QOpenGLFunctions_2_1>
#include <QOpenGLVertexArrayObject>
#include <QGLBuffer>
#include <QOpenGLShaderProgram>

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Qt/CreateOpenGLContext.h>
#include <CGAL/Random.h>

#include <vector>
#include <cstdlib>

#define NB_VBO_BUFFERS 8
#define NB_VAO_BUFFERS 4

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3  Local_point;
typedef Local_kernel::Vector_3 Local_vector;

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctor
{
  template<typename LCC>
  static CGAL::Color run(const LCC& alcc,
                         typename LCC::Dart_const_handle dh)
  {
    if (dh==alcc.null_handle) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    // Here dh is the smaller dart of its face
    CGAL::Random random(alcc.darts().index(dh));
    return CGAL::Color(random.get_int(0,256),
                       random.get_int(0,256),
                       random.get_int(0,256));
  }
};

//Vertex source code
const char vertex_source_mono[] =
  {
    "#version 120 \n"
    "attribute highp vec4 vertex;\n"
    "attribute highp vec3 normal;\n"
    
    "uniform highp mat4 mvp_matrix;\n"
    "uniform highp mat4 mv_matrix; \n"
    
    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "void main(void)\n"
    "{\n"
    "   fP = mv_matrix * vertex; \n"
    "   fN = mat3(mv_matrix)* normal; \n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
  };

const char vertex_source_multi[] =
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

//Vertex source code
const char fragment_source_mono[] =
  {
    "#version 120 \n"
    "varying highp vec4 fP; \n"
    "varying highp vec3 fN; \n"
    "uniform highp vec4 color; \n"
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
    "   vec4 diffuse = max(dot(N,L), 0.0) * light_diff * color; \n"
    "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"
    
    "gl_FragColor = light_amb*color + diffuse  ; \n"
    "} \n"
    "\n"
  };

const char fragment_source_multi[] =
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

//Vertex source code
const char vertex_source_p_l[] =
  {
    "#version 120 \n"
    "attribute highp vec4 vertex;\n"
    "uniform highp mat4 mvp_matrix;\n"
    "void main(void)\n"
    "{\n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
  };
//Vertex source code
const char fragment_source_p_l[] =
  {
    "#version 120 \n"
    "uniform highp vec4 color; \n"
    "void main(void) { \n"
    "gl_FragColor = color; \n"
    "} \n"
    "\n"
  };

template<class LCC, int dim=LCC::ambient_dimension>
struct Geom_utils;

template<class LCC>
struct Geom_utils<LCC,3>
{
  Local_point get_point(const LCC& lcc,
                        typename LCC::Vertex_attribute_const_handle vh)
  { return converter(lcc.point_of_vertex_attribute(vh)); }

  Local_point get_point(const LCC& lcc, typename LCC::Dart_const_handle dh)
  { return converter(lcc.point(dh)); }

  Local_vector get_facet_normal(const LCC& lcc,
                                typename LCC::Dart_const_handle dh)
  {
    Local_vector n = converter(CGAL::compute_normal_of_cell_2<LCC>(lcc,dh));
    n = n/(CGAL::sqrt(n*n));
    return n;
  }

  Local_vector get_vertex_normal(const LCC& lcc,
                                 typename LCC::Dart_const_handle dh)
  {
    Local_vector n = converter(CGAL::compute_normal_of_cell_0<LCC>(lcc,dh));
    n = n/(CGAL::sqrt(n*n));
    return n;
  }
protected:
  CGAL::Cartesian_converter<typename LCC::Traits, Local_kernel> converter;
};

template<class LCC>
struct Geom_utils<LCC,2>
{
  Local_point get_point(const LCC& lcc,
                        typename LCC::Vertex_attribute_const_handle vh)
  {
    Local_point p(converter(lcc.point_of_vertex_attribute(vh).x()),0,
                  converter(lcc.point_of_vertex_attribute(vh).y()));
    return p;
  }

  Local_point get_point(const LCC& lcc, typename LCC::Dart_const_handle dh)
  { return get_point(lcc, lcc.vertex_attribute(dh)); }

  Local_vector get_facet_normal(const LCC&, typename LCC::Dart_const_handle)
  {
    Local_vector n(0,1,0);
    return n;
  }

  Local_vector get_vertex_normal(const LCC&, typename LCC::Dart_const_handle)
  {
    Local_vector n(0,1,0);
    return n;
  }
protected:
  CGAL::Cartesian_converter<typename LCC::Traits, Local_kernel> converter;
};

template<class LCC>
CGAL::Bbox_3 bbox(const LCC& lcc)
{
  CGAL::Bbox_3 bb;
  Geom_utils<LCC> geomutils;

  typename LCC::Vertex_attribute_range::const_iterator
    it=lcc.vertex_attributes().begin(), itend=lcc.vertex_attributes().end();
  if ( it!=itend )
  {
    bb = geomutils.get_point(lcc, it).bbox();
    for( ++it; it!=itend; ++it)
    {
      bb = bb + geomutils.get_point(lcc, it).bbox();
    }
  }

  return bb;
}

template<class LCC, class ColorFunctor>
class SimpleLCCViewerQt : public QGLViewer, public QOpenGLFunctions_2_1
{
  typedef typename LCC::Dart_const_handle Dart_const_handle;

  struct Vertex_info
  {
    Dart_const_handle dh;
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
  SimpleLCCViewerQt(const LCC& alcc) :
    QGLViewer(CGAL::Qt::createOpenGLContext()),
    lcc(alcc),
    wireframe(false),
    flatShading(true),
    edges(true),
    vertices(true),
    mono(false),
    inverse_normal(false),
    size_points(7.),
    size_edges(3.1),
    ambient(0.6f, 0.5f, 0.5f, 0.5f),
    are_buffers_initialized(false)
  {
    setWindowTitle("3D lcc viewer");
    resize(500, 450);

    QGLFormat newFormat = this->format();
    newFormat.setSampleBuffers(true);
    newFormat.setSamples(16);
    this->setFormat(newFormat);
}

  ~SimpleLCCViewerQt()
  {
    for (int i=0; i<NB_VBO_BUFFERS; ++i)
      buffers[i].destroy();

    for (int i=0; i<NB_VAO_BUFFERS; ++i)
      vao[i].destroy();
  }

protected:
  void compile_shaders()
  {
    rendering_program.removeAllShaders();
    rendering_program_p_l.removeAllShaders();
    
    // Create the buffers
    for (int i=0; i<NB_VBO_BUFFERS; ++i)
      if(!buffers[i].isCreated() && !buffers[i].create())
      {
        std::cerr<<"VBO Creation number "<<i<<" FAILED"<<std::endl;
      }
    
    for (int i=0; i<NB_VAO_BUFFERS; ++i)
      if(!vao[i].isCreated() && !vao[i].create())
      {
        std::cerr<<"VAO Creation number "<<i<<" FAILED"<<std::endl;
      }

    //The Facets
    QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if (mono)
    {
      if(!vertex_shader->compileSourceCode(vertex_source_mono))
      {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
      }
    }
    else
    {
      if(!vertex_shader->compileSourceCode(vertex_source_multi))
      {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
      }
    }
    
    QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if (mono)
    {
      if(!fragment_shader->compileSourceCode(fragment_source_mono))
      {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
      }
    }
    else
    {
      if(!fragment_shader->compileSourceCode(fragment_source_multi))
      {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
      }
    }
    
    if(!rendering_program.addShader(vertex_shader))
    {
      std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program.addShader(fragment_shader))
    {
      std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    
    if(!rendering_program.link())
    {
      std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program.bind();

    vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader->compileSourceCode(vertex_source_p_l))
    {
      std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader->compileSourceCode(fragment_source_p_l))
    {
      std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program_p_l.addShader(vertex_shader))
    {
      std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program_p_l.addShader(fragment_shader))
    {
      std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program_p_l.link())
    {
      std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program_p_l.bind();
  }

  void initialize_buffers()
  {
    int bufn = 0;
    int vaon = 0;
    
    //points of the facets
    assert(vaon<NB_VAO_BUFFERS);
    vao[vaon].bind();
    
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(pos_facets.data(),
                        static_cast<int>(pos_facets.size()*sizeof(float)));
    vertexLocation[vaon] = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.enableAttributeArray(vertexLocation[vaon]);
    rendering_program.setAttributeBuffer(vertexLocation[vaon],GL_FLOAT,0,3);
    rendering_program.release();

    buffers[bufn].release();
    ++bufn;
    
    //normals of the facets
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(flat_normals.data(),
                        static_cast<int>(flat_normals.size()*sizeof(float)));
    normalsLocation = rendering_program.attributeLocation("normal");
    rendering_program.bind();
    rendering_program.enableAttributeArray(normalsLocation);
    rendering_program.setAttributeBuffer(normalsLocation,GL_FLOAT,0,3);
    rendering_program.release();

    buffers[bufn].release();
    ++bufn;

    if (!mono)
    {
      //colors of the facets
      assert(bufn<NB_VBO_BUFFERS);
      buffers[bufn].bind();
      buffers[bufn].allocate(colors.data(),
                             static_cast<int>(colors.size()*sizeof(float)));
      colorsLocation = rendering_program.attributeLocation("color");
      rendering_program.bind();
      rendering_program.enableAttributeArray(colorsLocation);
      rendering_program.setAttributeBuffer(colorsLocation,GL_FLOAT,0,3);
      rendering_program.release();

      buffers[bufn].release();      
      ++bufn;
    }
    
    vao[vaon].release();
    ++vaon;
    
    assert(vaon<NB_VAO_BUFFERS);
    vao[vaon].bind();

    //points of the facets
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(pos_facets.data(),
                        static_cast<int>(pos_facets.size()*sizeof(float)));
    vertexLocation[vaon] = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.enableAttributeArray(vertexLocation[vaon]);
    rendering_program.setAttributeBuffer(vertexLocation[vaon],GL_FLOAT,0,3);
    rendering_program.release();

    buffers[bufn].release();
    ++bufn;

    //normals of the facets
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(smooth_normals.data(),
                        static_cast<int>(smooth_normals.size()*sizeof(float)));
    normalsLocation = rendering_program.attributeLocation("normal");
    rendering_program.bind();
    rendering_program.enableAttributeArray(normalsLocation);
    rendering_program.setAttributeBuffer(normalsLocation,GL_FLOAT,0,3);
    rendering_program.release();

    buffers[bufn].release();
    ++bufn;

    if (!mono)
    {
      //colors of the facets
      assert(bufn<NB_VBO_BUFFERS);
      buffers[bufn].bind();
      buffers[bufn].allocate(colors.data(), static_cast<int>(colors.size()*sizeof(float)));
      colorsLocation = rendering_program.attributeLocation("color");
      rendering_program.bind();
      rendering_program.enableAttributeArray(colorsLocation);
      rendering_program.setAttributeBuffer(colorsLocation,GL_FLOAT,0,3);
      buffers[bufn].release();

      rendering_program.release();
      ++bufn;
    }
    
    vao[vaon].release();
    ++vaon;
    
    //The lines
    assert(vaon<NB_VAO_BUFFERS);
    vao[vaon].bind();
    
    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(pos_lines.data(), static_cast<int>(pos_lines.size()*sizeof(float)));
    vertexLocation[vaon] = rendering_program_p_l.attributeLocation("vertex");
    rendering_program_p_l.bind();
    rendering_program_p_l.enableAttributeArray(vertexLocation[vaon]);
    rendering_program_p_l.setAttributeBuffer(vertexLocation[vaon],GL_FLOAT,0,3);
    rendering_program_p_l.release();

    buffers[bufn].release();
    ++bufn;
    
    vao[vaon].release();
    ++vaon;
    
    //The points
    assert(vaon<NB_VAO_BUFFERS);
    vao[vaon].bind();

    assert(bufn<NB_VBO_BUFFERS);
    buffers[bufn].bind();
    buffers[bufn].allocate(pos_points.data(), static_cast<int>(pos_points.size()*sizeof(float)));
    vertexLocation[vaon] = rendering_program_p_l.attributeLocation("vertex");
    rendering_program_p_l.bind();
    rendering_program_p_l.enableAttributeArray(vertexLocation[vaon]);
    rendering_program_p_l.setAttributeBuffer(vertexLocation[vaon],GL_FLOAT,0,3);
    rendering_program_p_l.release();

    buffers[bufn].release();
    ++bufn;

    vao[vaon].release();
    ++vaon;
    
    are_buffers_initialized = true;
  }

  void compute_face(Dart_const_handle dh)
  {
    // We fill only closed faces.
    Dart_const_handle cur=dh;
    Dart_const_handle min=dh;
    do
    {
      if (!lcc.is_next_exist(cur)) return; // open face=>not filled
      if (cur<min) min=cur;
      cur=lcc.next(cur);
    }
    while(cur!=dh);
    
    double r=0.;
    double g=0.;
    double b=0.;

    CGAL::Color c=ColorFunctor::run(lcc, dh);
    r=(double)c.red()/(double)255;
    g=(double)c.green()/(double)255;
    b=(double)c.blue()/(double)255;
    
    //compute flat normals
    Local_vector normal = geomutils.get_facet_normal(lcc, dh);
    normal = normal/(CGAL::sqrt(normal*normal));
    if (inverse_normal) normal=normal*-1;
    
    if (lcc.next(lcc.next(lcc.next(dh)))!=dh)
    {
      try
      {
        P_traits cdt_traits(normal);
        CDT cdt(cdt_traits);
        
        // Iterates on the vector of facet handles
        typename CDT::Vertex_handle previous = NULL, first = NULL;
        cur=dh;
        do
        {
          typename CDT::Vertex_handle vh = cdt.insert(geomutils.get_point(lcc, cur));
          if(first==NULL)
          { first=vh; }
          vh->info().v = geomutils.get_vertex_normal(lcc, cur);
          if (inverse_normal) vh->info().v=vh->info().v*-1;
          
          if(previous!=NULL && previous!=vh)
          { cdt.insert_constraint(previous, vh); }
          previous=vh;
          cur=lcc.next(cur);
        }
        while(cur!=dh);
        
        if (previous!=NULL)
          cdt.insert_constraint(previous, first);
        
        // sets mark is_external
        for(typename CDT::All_faces_iterator fit = cdt.all_faces_begin(),
              fitend = cdt.all_faces_end(); fit!=fitend; ++fit)
        {
          fit->info().is_external = true;
          fit->info().is_process = false;
        }
        //check if the facet is external or internal
        std::queue<typename CDT::Face_handle> face_queue;
        typename CDT::Face_handle face_internal = NULL;
        face_queue.push(cdt.infinite_vertex()->face());
        while(! face_queue.empty() )
        {
          typename CDT::Face_handle fh = face_queue.front();
          face_queue.pop();
          if(!fh->info().is_process)
          {
            fh->info().is_process = true;
            for(int i = 0; i <3; ++i)
            {
              if(!cdt.is_constrained(std::make_pair(fh, i)))
              {
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
            for(int i = 0; i <3; ++i)
            {
              if(!cdt.is_constrained(std::make_pair(fh, i)))
              {
                face_queue.push(fh->neighbor(i));
              }
            }
          }
        }

        //iterates on the internal faces to add the vertices to the positions
        //and the normals to the appropriate vectors
        for(typename CDT::Finite_faces_iterator ffit = cdt.finite_faces_begin(),
              ffitend = cdt.finite_faces_end(); ffit != ffitend; ++ffit)
        {
          if(!ffit->info().is_external)
          {
            flat_normals.push_back(normal.x());
            flat_normals.push_back(normal.y());
            flat_normals.push_back(normal.z());
            
            flat_normals.push_back(normal.x());
            flat_normals.push_back(normal.y());
            flat_normals.push_back(normal.z());
            
            flat_normals.push_back(normal.x());
            flat_normals.push_back(normal.y());
            flat_normals.push_back(normal.z());
            
            smooth_normals.push_back(ffit->vertex(0)->info().v.x());
            smooth_normals.push_back(ffit->vertex(0)->info().v.y());
            smooth_normals.push_back(ffit->vertex(0)->info().v.z());
            
            smooth_normals.push_back(ffit->vertex(1)->info().v.x());
            smooth_normals.push_back(ffit->vertex(1)->info().v.y());
            smooth_normals.push_back(ffit->vertex(1)->info().v.z());
            
            smooth_normals.push_back(ffit->vertex(2)->info().v.x());
            smooth_normals.push_back(ffit->vertex(2)->info().v.y());
            smooth_normals.push_back(ffit->vertex(2)->info().v.z());
            
            pos_facets.push_back(ffit->vertex(0)->point().x());
            pos_facets.push_back(ffit->vertex(0)->point().y());
            pos_facets.push_back(ffit->vertex(0)->point().z());
            
            pos_facets.push_back(ffit->vertex(1)->point().x());
            pos_facets.push_back(ffit->vertex(1)->point().y());
            pos_facets.push_back(ffit->vertex(1)->point().z());
            
            pos_facets.push_back(ffit->vertex(2)->point().x());
            pos_facets.push_back(ffit->vertex(2)->point().y());
            pos_facets.push_back(ffit->vertex(2)->point().z());
            
            colors.push_back(r);colors.push_back(g);colors.push_back(b);
            colors.push_back(r);colors.push_back(g);colors.push_back(b);
            colors.push_back(r);colors.push_back(g);colors.push_back(b);
          }
        }
      }
      catch(...)
      { // Triangulation crash: the face is not filled
      }
    }
    else
    {
      colors.push_back(r);colors.push_back(g);colors.push_back(b);
      colors.push_back(r);colors.push_back(g);colors.push_back(b);
      colors.push_back(r);colors.push_back(g);colors.push_back(b);
      
      flat_normals.push_back(normal.x());
      flat_normals.push_back(normal.y());
      flat_normals.push_back(normal.z());

      flat_normals.push_back(normal.x());
      flat_normals.push_back(normal.y());
      flat_normals.push_back(normal.z());

      flat_normals.push_back(normal.x());
      flat_normals.push_back(normal.y());
      flat_normals.push_back(normal.z());

      //compute Smooth normals
      cur=dh;
      do
      {
        Local_vector normal = geomutils.get_vertex_normal(lcc, cur);
        normal = normal/(CGAL::sqrt(normal*normal));
        if (inverse_normal) normal=normal*-1;

        smooth_normals.push_back(normal.x());
        smooth_normals.push_back(normal.y());
        smooth_normals.push_back(normal.z());
        
        Local_point p = geomutils.get_point(lcc, cur);
        pos_facets.push_back(p.x());
        pos_facets.push_back(p.y());
        pos_facets.push_back(p.z());
        cur=lcc.next(cur);
      }
      while(cur!=dh);
    }
  }

  void compute_edge(Dart_const_handle dh)
  {
    Local_point p =  geomutils.get_point(lcc, dh);
    Dart_const_handle d2 = lcc.other_extremity(dh);
    if ( d2!=NULL )
    {
      Local_point p2 = geomutils.get_point(lcc, d2);
      pos_lines.push_back(p.x());
      pos_lines.push_back(p.y());
      pos_lines.push_back(p.z());

      pos_lines.push_back(p2.x());
      pos_lines.push_back(p2.y());
      pos_lines.push_back(p2.z());
    }
  }

  void compute_vertex(Dart_const_handle dh, bool& empty)
  {
    Local_point p =  geomutils.get_point(lcc, dh);
    pos_points.push_back(p.x());
    pos_points.push_back(p.y());
    pos_points.push_back(p.z());

    if ( empty )
    {
      bb = p.bbox();
      empty = false;
    }
    else
      bb = bb + p.bbox();
  }

  void compute_elements()
  {
    pos_facets.clear();
    flat_normals.clear();
    smooth_normals.clear();
    pos_lines.clear();
    pos_points.clear();
    colors.clear();
    
    if ( lcc.is_empty() )
    {
      bb = Local_point(CGAL::ORIGIN).bbox();
      bb = bb + Local_point(1,1,1).bbox(); // To avoid a warning from Qglviewer
      return;
    }

    unsigned int markfaces    = lcc.get_new_mark();
    unsigned int markedges    = lcc.get_new_mark();
    unsigned int markvertices = lcc.get_new_mark();

    bool empty = true;
    for (typename LCC::Dart_range::const_iterator it=lcc.darts().begin(),
           itend=lcc.darts().end(); it!=itend; ++it )
    {
      if ( !lcc.is_marked(it, markfaces) )
      {
        compute_face(it);
        CGAL::mark_cell<LCC, 2>(lcc, it, markfaces);
      }

      if ( !lcc.is_marked(it, markedges) )
      {
        compute_edge(it);
        CGAL::mark_cell<LCC, 1>(lcc, it, markedges);
      }

      if ( !lcc.is_marked(it, markvertices) )
      {
        compute_vertex(it, empty);
        CGAL::mark_cell<LCC, 0>(lcc, it, markvertices);
      }
    }

    lcc.free_mark(markfaces);
    lcc.free_mark(markedges);
    lcc.free_mark(markvertices);
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

    rendering_program.bind();
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    mvLocation = rendering_program.uniformLocation("mv_matrix");
    lightLocation[0] = rendering_program.uniformLocation("light_pos");
    lightLocation[1] = rendering_program.uniformLocation("light_diff");
    lightLocation[2] = rendering_program.uniformLocation("light_spec");
    lightLocation[3] = rendering_program.uniformLocation("light_amb");
    lightLocation[4] = rendering_program.uniformLocation("spec_power");

    rendering_program.setUniformValue(lightLocation[0], position);
    rendering_program.setUniformValue(lightLocation[1], diffuse);
    rendering_program.setUniformValue(lightLocation[2], specular);
    rendering_program.setUniformValue(lightLocation[3], ambient);
    rendering_program.setUniformValue(lightLocation[4], shininess);
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);
    rendering_program.setUniformValue(mvLocation, mvMatrix);

    if (mono)
      colorLocation2 = rendering_program.uniformLocation("color");
    
    rendering_program.release();

    rendering_program_p_l.bind();
    mvpLocation[1] = rendering_program_p_l.uniformLocation("mvp_matrix");
    colorLocation = rendering_program_p_l.uniformLocation("color");
    rendering_program.setUniformValue(mvpLocation[1], mvpMatrix);
    rendering_program_p_l.release();
  }

  virtual void draw()
  {
    glEnable(GL_DEPTH_TEST);
    if(!are_buffers_initialized)
      initialize_buffers();

    QColor color;
    if ( !wireframe )
    {
      if(flatShading)
      {
        vao[0].bind();
        attrib_buffers(this);
        rendering_program.bind();
        if (mono)
        {
          CGAL::Color c=ColorFunctor::run(lcc, NULL);
          color.setRgbF((double)c.red()/(double)255,
                        (double)c.green()/(double)255,
                        (double)c.blue()/(double)255);
          rendering_program.setUniformValue(colorLocation2,color);
        }      
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_facets.size()/3));
        rendering_program.release();
        vao[0].release();
      }
      else
      {
        vao[1].bind();
        attrib_buffers(this);
        rendering_program.bind();
        if (mono)
        {
          CGAL::Color c=ColorFunctor::run(lcc, NULL);
          color.setRgbF((double)c.red()/(double)255,
                        (double)c.green()/(double)255,
                        (double)c.blue()/(double)255);
          rendering_program.setUniformValue(colorLocation2,color);
        }
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_facets.size()/3));
        rendering_program.release();
        vao[1].release();
      }
    }
    if(edges)
    {
      vao[2].bind();
      attrib_buffers(this);
      color.setRgbF(0.2f, 0.2f, 0.7f);
      rendering_program_p_l.bind();
      rendering_program_p_l.setAttributeValue(colorLocation,color);
      ::glLineWidth(size_edges);
      glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_lines.size()/3));
      rendering_program_p_l.release();
      vao[2].release();
    }
    if(vertices)
    {
      ::glPointSize(size_points);
      vao[3].bind();
      attrib_buffers(this);
      color.setRgbF(.2f,.2f,.6f);
      rendering_program_p_l.bind();
      rendering_program_p_l.setAttributeValue(colorLocation,color);
      glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pos_points.size()/3));
      rendering_program_p_l.release();
      vao[3].release();
    }
  }

  virtual void init()
  {
    // Restore previous viewer state.
    restoreStateFromFile();

    // Define 'Control+Q' as the new exit shortcut (default was 'Escape')
    setShortcut(EXIT_VIEWER, Qt::CTRL+Qt::Key_Q);

    // Add custom key description (see keyPressEvent).
    setKeyDescription(Qt::Key_W, "Toggles wire frame display");
    setKeyDescription(Qt::Key_F, "Toggles flat shading display");
    setKeyDescription(Qt::Key_E, "Toggles edges display");
    setKeyDescription(Qt::Key_V, "Toggles vertices display");
    setKeyDescription(Qt::Key_M, "Toggles mono color for all faces");
    setKeyDescription(Qt::Key_N, "Inverse direction of normals");
    setKeyDescription(Qt::Key_Plus, "Increase size of edges");
    setKeyDescription(Qt::Key_Minus, "Decrease size of edges");
    setKeyDescription(Qt::Key_Plus+Qt::ShiftModifier, "Increase size of vertices");
    setKeyDescription(Qt::Key_Minus+Qt::ShiftModifier, "Decrease size of vertices");
    setKeyDescription(Qt::Key_PageDown, "Increase light (all colors, use shift/alt/ctrl for one rgb component)");
    setKeyDescription(Qt::Key_PageUp, "Decrease light (all colors, use shift/alt/ctrl for one rgb component)");

    // Light default parameters
    ::glLineWidth(size_edges);
    ::glPointSize(size_points);
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

    initializeOpenGLFunctions();
    compile_shaders();
    compute_elements();

    this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(),
                                                       bb.ymin(),
                                                       bb.zmin()),
                                        qglviewer::Vec(bb.xmax(),
                                                       bb.ymax(),
                                                       bb.zmax()));

    this->showEntireScene();
  }

  void keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    if ((e->key()==Qt::Key_W) && (modifiers==Qt::NoButton))
    {
      wireframe = !wireframe;
    if (wireframe)
    {
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      displayMessage("Wireframe.");
    }
    else
    {
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      displayMessage("Filled faces.");
    }
      updateGL();
    }
    else if ((e->key()==Qt::Key_F) && (modifiers==Qt::NoButton))
    {
      flatShading = !flatShading;
      if (flatShading)
        displayMessage("Flat shading.");
      else
        displayMessage("Gouraud shading.");
      updateGL();
    }
    else if ((e->key()==Qt::Key_E) && (modifiers==Qt::NoButton))
    {
      edges = !edges;
      displayMessage(QString("Draw edges=%1.").arg(edges?"true":"false"));
      updateGL();
    }
    else if ((e->key()==Qt::Key_V) && (modifiers==Qt::NoButton))
    {
      vertices = !vertices;
      displayMessage(QString("Draw vertices=%1.").arg(vertices?"true":"false"));
      updateGL();
    }
    else if ((e->key()==Qt::Key_M) && (modifiers==Qt::NoButton))
    {
      mono = !mono;
      initialize_buffers();
      compile_shaders();
      displayMessage(QString("Mono color=%1.").arg(mono?"true":"false"));
      updateGL();
    }
    else if ((e->key()==Qt::Key_N) && (modifiers==Qt::NoButton))
    {
      inverse_normal = !inverse_normal;
      compute_elements();
      initialize_buffers();
      compile_shaders();
      displayMessage(QString("Inverse normal=%1.").arg(inverse_normal?"true":"false"));
      updateGL();
    }
    else if ((e->key()==Qt::Key_Plus) && (modifiers==Qt::KeypadModifier))
    {
      size_edges+=.5;
      displayMessage(QString("Size of edges=%1.").arg(size_edges));
      updateGL();
    }
    else if ((e->key()==Qt::Key_Minus) && (modifiers==Qt::KeypadModifier))
    {
      if (size_edges>.5) size_edges-=.5;
      displayMessage(QString("Size of edges=%1.").arg(size_edges));
      updateGL();
    }
    else if ((e->key()==Qt::Key_Plus) && (modifiers==(Qt::ShiftModifier|Qt::KeypadModifier)))
    {
      size_points+=.5;
      displayMessage(QString("Size of points=%1.").arg(size_points));
      updateGL();
    }
    else if ((e->key()==Qt::Key_Minus) && (modifiers==(Qt::ShiftModifier|Qt::KeypadModifier)))
    {
      if (size_points>.5) size_points-=.5;
      displayMessage(QString("Size of points=%1.").arg(size_points));
      updateGL();
    }
    else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton))
    {
      ambient.setX(ambient.x()+.1);
      if (ambient.x()>1.) ambient.setX(1.);
      ambient.setY(ambient.x()+.1);
      if (ambient.y()>1.) ambient.setY(1.);
      ambient.setZ(ambient.x()+.1);
      if (ambient.z()>1.) ambient.setZ(1.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
      updateGL();
    }
    else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::NoButton))
    {
      ambient.setX(ambient.x()-.1);
      if (ambient.x()<0.) ambient.setX(0.);
      ambient.setY(ambient.y()-.1);
      if (ambient.y()<0.) ambient.setY(0.);
      ambient.setZ(ambient.z()-.1);
      if (ambient.z()<0.) ambient.setZ(0.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
      updateGL();
    }
    else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::ShiftModifier))
    {
      ambient.setX(ambient.x()+.1);
      if (ambient.x()>1.) ambient.setX(1.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
      updateGL();
    }
    else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::AltModifier))
    {
      ambient.setY(ambient.y()+.1);
      if (ambient.y()>1.) ambient.setY(1.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
      updateGL();
    }
    else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::ControlModifier))
    {
      ambient.setZ(ambient.z()+.1);
      if (ambient.z()>1.) ambient.setZ(1.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
      updateGL();
    }
    else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::ShiftModifier))
    {
      ambient.setX(ambient.x()-.1);
      if (ambient.x()<0.) ambient.setX(0.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
      updateGL();
    }
    else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::AltModifier))
    {
      ambient.setY(ambient.y()-.1);
      if (ambient.y()<0.) ambient.setY(0.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
      updateGL();
    }
    else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::ControlModifier))
    {
      ambient.setZ(ambient.z()-.1);
      if (ambient.z()<0.) ambient.setZ(0.);
      displayMessage(QString("Light color=(%1 %2 %3).").
                     arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
      updateGL();
    }
    else
      QGLViewer::keyPressEvent(e);
  }


  virtual QString helpString() const
  {
    QString text("<h2>L C C   V i e w e r</h2>");
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
  const LCC& lcc;
  
  bool wireframe;
  bool flatShading;
  bool edges;
  bool vertices;
  bool mono;
  bool inverse_normal;
  
  double size_points;
  double size_edges;

  QVector4D ambient;
  
  bool are_buffers_initialized;

  //Shaders elements
  int vertexLocation[NB_VAO_BUFFERS];
  int normalsLocation;
  int mvpLocation[2];
  int mvLocation;
  int colorLocation;
  int colorLocation2;
  int lightLocation[5];

  std::vector<float> pos_points;
  std::vector<float> pos_lines;
  std::vector<float> pos_facets;
  std::vector<float> smooth_normals;
  std::vector<float> flat_normals;
  std::vector<float> colors;

  QGLBuffer buffers[NB_VBO_BUFFERS];
  QOpenGLVertexArrayObject vao[NB_VAO_BUFFERS];
  int colorsLocation;

  QOpenGLShaderProgram rendering_program;
  QOpenGLShaderProgram rendering_program_p_l;

  CGAL::Bbox_3 bb;
  Geom_utils<LCC> geomutils;
};

template<class LCC, class ColorFunctor=DefaultColorFunctor>
void display_lcc(const LCC& alcc)
{
  int argc=1;

  const char* argv[2]={"lccviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimpleLCCViewerQt<LCC, ColorFunctor> mainwindow(alcc);
  mainwindow.show();

  app.exec();
}

#endif // CGAL_LCC_3_VIEWER_QT_H
