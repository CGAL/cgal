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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Kumar Snehasish <kumar.snehasish@gmail.com>
//
#include "Viewer.h"

#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Qt/CreateOpenGLContext.h>

#include <QGLViewer/vec.h>
#include <QDebug>

//Vertex source code
const char vertex_source[] =
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
const char fragment_source[] =
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

Viewer::Viewer(QWidget* parent)
  : QGLViewer(CGAL::Qt::createOpenGLContext(),parent),
    wireframe(false),
    flatShading(true),
    edges(true),
    vertices(true),
    inverse_normal(false),
    size_points(7.),
    size_edges(3.1),
    ambient(0.6f, 0.5f, 0.5f, 0.5f),
    m_previous_scene_empty(true),
    are_buffers_initialized(false)
{
}

Viewer::~Viewer()
{
  for (int i=0; i<NB_VBO_BUFFERS; ++i)
    buffers[i].destroy();
  
  for (int i=0; i<NB_VAO_BUFFERS; ++i)
    vao[i].destroy();
}

void Viewer::compile_shaders()
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
  if(!vertex_shader->compileSourceCode(vertex_source))
  {
    std::cerr<<"Compiling vertex source FAILED"<<std::endl;
  }

  QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
  if(!fragment_shader->compileSourceCode(fragment_source))
  {
    std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
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

void Viewer::initialize_buffers()
{
  //points of the facets
  vao[0].bind();
  buffers[0].bind();
  buffers[0].allocate(pos_facets.data(),
                      static_cast<int>(pos_facets.size()*sizeof(float)));
  vertexLocation[0] = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.enableAttributeArray(vertexLocation[0]);
  rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
  rendering_program.release();
  buffers[0].release();

  //normals of the facets
  buffers[1].bind();
  buffers[1].allocate(flat_normals.data(),
                      static_cast<int>(flat_normals.size()*sizeof(float)));
  normalsLocation = rendering_program.attributeLocation("normal");
  rendering_program.bind();
  rendering_program.enableAttributeArray(normalsLocation);
  rendering_program.setAttributeBuffer(normalsLocation,GL_FLOAT,0,3);
  buffers[1].release();

  //colors of the facets
  buffers[2].bind();
  buffers[2].allocate(colors.data(),
                      static_cast<int>(colors.size()*sizeof(float)));
  colorsLocation = rendering_program.attributeLocation("color");
  rendering_program.bind();
  rendering_program.enableAttributeArray(colorsLocation);
  rendering_program.setAttributeBuffer(colorsLocation,GL_FLOAT,0,3);
  buffers[2].release();
  rendering_program.release();

  vao[0].release();
  vao[1].bind();

  //points of the facets
  buffers[3].bind();
  buffers[3].allocate(pos_facets.data(), static_cast<int>(pos_facets.size()*sizeof(float)));
  vertexLocation[0] = rendering_program.attributeLocation("vertex");
  rendering_program.bind();
  rendering_program.enableAttributeArray(vertexLocation[0]);
  rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
  rendering_program.release();
  buffers[3].release();

  //normals of the facets
  buffers[4].bind();
  buffers[4].allocate(smooth_normals.data(),
                      static_cast<int>(smooth_normals.size()*sizeof(float)));
  normalsLocation = rendering_program.attributeLocation("normal");
  rendering_program.bind();
  rendering_program.enableAttributeArray(normalsLocation);
  rendering_program.setAttributeBuffer(normalsLocation,GL_FLOAT,0,3);
  buffers[4].release();

  //colors of the facets
  buffers[5].bind();
  buffers[5].allocate(colors.data(), static_cast<int>(colors.size()*sizeof(float)));
  colorsLocation = rendering_program.attributeLocation("color");
  rendering_program.bind();
  rendering_program.enableAttributeArray(colorsLocation);
  rendering_program.setAttributeBuffer(colorsLocation,GL_FLOAT,0,3);
  buffers[5].release();
  rendering_program.release();

  vao[1].release();

  //The lines
  vao[2].bind();
  buffers[6].bind();
  buffers[6].allocate(pos_lines.data(), static_cast<int>(pos_lines.size()*sizeof(float)));
  vertexLocation[2] = rendering_program_p_l.attributeLocation("vertex");
  rendering_program_p_l.bind();
  rendering_program_p_l.enableAttributeArray(vertexLocation[2]);
  rendering_program_p_l.setAttributeBuffer(vertexLocation[2],GL_FLOAT,0,3);
  buffers[6].release();
  rendering_program_p_l.release();
  vao[2].release();

  //The points
  vao[3].bind();
  buffers[7].bind();
  buffers[7].allocate(pos_points.data(), static_cast<int>(pos_points.size()*sizeof(float)));
  vertexLocation[2] = rendering_program_p_l.attributeLocation("vertex");
  rendering_program_p_l.bind();
  rendering_program_p_l.enableAttributeArray(vertexLocation[2]);
  rendering_program_p_l.setAttributeBuffer(vertexLocation[2],GL_FLOAT,0,3);
  buffers[7].release();
  rendering_program_p_l.release();
  vao[3].release();

  are_buffers_initialized = true;
}

void Viewer::compute_face(Dart_handle dh, LCC::size_type markface)
{
  LCC &lcc = *scene->lcc;

  CGAL::mark_cell<LCC, 2>(lcc, dh, markface);

  double r = (double)lcc.info<3>(dh).color().r()/255.0;
  double g = (double)lcc.info<3>(dh).color().g()/255.0;
  double b = (double)lcc.info<3>(dh).color().b()/255.0;
  if ( !lcc.is_free(dh, 3) )
  {
    r += (double)lcc.info<3>(lcc.beta(dh,3)).color().r()/255.0;
    g += (double)lcc.info<3>(lcc.beta(dh,3)).color().g()/255.0;
    b += (double)lcc.info<3>(lcc.beta(dh,3)).color().b()/255.0;
    r /= 2; g /= 2; b /= 2;
  }

  //compute flat normals
  LCC::Vector normal = CGAL::compute_normal_of_cell_2(lcc,dh);
  normal = normal/(CGAL::sqrt(normal*normal));
  if (inverse_normal)
    normal=normal*-1;
  
  if (lcc.beta<1,1,1>(dh)!=dh)
  {
    try // Try catch to avoir crash of triangulation
    {
      P_traits cdt_traits(normal);
      CDT cdt(cdt_traits);
      
      // Iterates on the vector of facet handles
      CDT::Vertex_handle previous = NULL, first = NULL;
      for (LCC::Dart_of_orbit_range<1>::const_iterator
             he_circ = lcc.darts_of_orbit<1>(dh).begin(),
             he_circ_end = lcc.darts_of_orbit<1>(dh).end();
           he_circ!=he_circ_end; ++he_circ)
      {
        CDT::Vertex_handle vh = cdt.insert(lcc.point(he_circ));
        if(first == NULL)
        { first = vh; }
        vh->info().v = CGAL::compute_normal_of_cell_0<LCC>(lcc, he_circ);
        if (inverse_normal) vh->info().v=vh->info().v*-1;
        if(previous!=NULL && previous != vh)
        { cdt.insert_constraint(previous, vh); }
        previous = vh;
      }
      if (previous!=NULL)
        cdt.insert_constraint(previous, first);

      // sets mark is_external
      for(CDT::All_faces_iterator fit = cdt.all_faces_begin(),
            fitend = cdt.all_faces_end(); fit!=fitend; ++fit)
      {
        fit->info().is_external = true;
        fit->info().is_process = false;
      }
      //check if the facet is external or internal
      std::queue<CDT::Face_handle> face_queue;
      CDT::Face_handle face_internal = NULL;
      face_queue.push(cdt.infinite_vertex()->face());
      while(! face_queue.empty() )
      {
        CDT::Face_handle fh = face_queue.front();
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
        CDT::Face_handle fh = face_queue.front();
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
      for(CDT::Finite_faces_iterator ffit = cdt.finite_faces_begin(),
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
  { // The face is a triangle
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

    for (LCC::Dart_of_orbit_range<1>::const_iterator
           orbitIter = lcc.darts_of_orbit<1>(dh).begin();
         orbitIter.cont(); ++orbitIter)
    {
      //compute Smooth normals
      LCC::Vector normal = CGAL::compute_normal_of_cell_0(lcc,orbitIter);
      normal = normal/(CGAL::sqrt(normal*normal));
      if (inverse_normal) normal=normal*-1;

      smooth_normals.push_back(normal.x());
      smooth_normals.push_back(normal.y());
      smooth_normals.push_back(normal.z());

      const LCC::Point& p = lcc.point(orbitIter);
      pos_facets.push_back(p.x());
      pos_facets.push_back(p.y());
      pos_facets.push_back(p.z());
    }
  }
}

void Viewer::compute_edge(Dart_handle dh, LCC::size_type markedge)
{
  LCC &lcc = *scene->lcc;

  CGAL::mark_cell<LCC, 1>(lcc, dh, markedge);

  const LCC::Point& p =  lcc.point(dh);
  Dart_handle d2 = lcc.other_extremity(dh);
  if ( d2!=NULL )
  {
    const LCC::Point& p2 = lcc.point(d2);
    pos_lines.push_back(p.x());
    pos_lines.push_back(p.y());
    pos_lines.push_back(p.z());

    pos_lines.push_back(p2.x());
    pos_lines.push_back(p2.y());
    pos_lines.push_back(p2.z());
  }
}

void Viewer::compute_vertex(Dart_handle dh, LCC::size_type markvertex, bool& empty)
{
  LCC &lcc = *scene->lcc;

  CGAL::mark_cell<LCC, 0>(lcc, dh, markvertex);

  const LCC::Point& p =  lcc.point(dh);
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

void Viewer::compute_elements()
{
  LCC &lcc = *scene->lcc;

  pos_facets.clear();
  flat_normals.clear();
  smooth_normals.clear();
  colors.clear();
  pos_lines.clear();
  pos_points.clear();

  if ( lcc.is_empty() )
  {
    bb = LCC::Point(CGAL::ORIGIN).bbox();
    bb = bb + LCC::Point(1,1,1).bbox(); // To avoid a warning from Qglviewer
    return;
  }

  LCC::size_type markvertex = lcc.get_new_mark();
  LCC::size_type markedge   = lcc.get_new_mark();
  LCC::size_type markface   = lcc.get_new_mark();

  bool empty = true;
  for (LCC::Attribute_range<3>::type::iterator it=lcc.attributes<3>().begin(),
         itend=lcc.attributes<3>().end(); it!=itend; ++it )
  {
    if ( it->info().is_visible() )
    {
      for(LCC::Dart_of_cell_range<3>::iterator
            dartIter=lcc.darts_of_cell<3>(lcc.dart_of_attribute<3>(it)).begin();
          dartIter.cont(); ++dartIter)
      {
        if ( it->info().is_filled() && !lcc.is_marked(dartIter, markface) )
          compute_face(dartIter, markface);

        if ( !lcc.is_marked(dartIter, markedge) )
          compute_edge(dartIter, markedge);

        if ( !lcc.is_marked(dartIter, markvertex) )
        compute_vertex(dartIter, markvertex, empty);
      }
    }
  }

  if ( empty )
  {
    bb = LCC::Point(CGAL::ORIGIN).bbox();
    bb = bb + LCC::Point(1,1,1).bbox(); // To avoid a warning from Qglviewer
  }

  for (LCC::Dart_range::iterator it=lcc.darts().begin(),
         itend=lcc.darts().end(); it!=itend; ++it )
  {
    lcc.unmark(it, markvertex);
    lcc.unmark(it, markedge);
    lcc.unmark(it, markface);
  }

  lcc.free_mark(markvertex);
  lcc.free_mark(markedge);
  lcc.free_mark(markface);
}

void Viewer::attrib_buffers(QGLViewer* viewer)
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

  rendering_program.release();
  rendering_program_p_l.bind();
  mvpLocation[1] = rendering_program_p_l.uniformLocation("mvp_matrix");
  colorLocation = rendering_program_p_l.uniformLocation("color");
  rendering_program.setUniformValue(mvpLocation[1], mvpMatrix);
  rendering_program_p_l.release();
}

void Viewer::sceneChanged()
{
  compute_elements();
  this->camera()->setSceneBoundingBox(qglviewer::Vec(bb.xmin(),
                                                     bb.ymin(),
                                                     bb.zmin()),
                                      qglviewer::Vec(bb.xmax(),
                                                     bb.ymax(),
                                                     bb.zmax()));
  are_buffers_initialized = false;

  if (m_previous_scene_empty)
    this->showEntireScene();
  else
#if QGLVIEWER_VERSION >= 0x020700
    this->update();
#else
    this->updateGL();
#endif

  m_previous_scene_empty = scene->lcc->is_empty(); // for the next call to sceneChanged
}

void Viewer::draw()
{
  if(scene)
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
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_facets.size()/3));
        rendering_program.release();
        vao[0].release();
      }
      else
      {
        vao[1].bind();
        attrib_buffers(this);
        rendering_program.bind();
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
      glLineWidth(size_edges);
      glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_lines.size()/3));
      rendering_program_p_l.release();
      vao[2].release();
    }
    if(vertices)
    {
      glPointSize(size_points);
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
}
void Viewer::init()
{
  // Restore previous viewer state.
  restoreStateFromFile();
  initializeOpenGLFunctions();
  // Define 'Control+Q' as the new exit shortcut (default was 'Escape')
  setShortcut(EXIT_VIEWER, Qt::CTRL+Qt::Key_Q);

  // Add custom key description (see keyPressEvent).
  setKeyDescription(Qt::Key_W, "Toggles wire frame display");
  setKeyDescription(Qt::Key_F, "Toggles flat shading display");
  setKeyDescription(Qt::Key_E, "Toggles edges display");
  setKeyDescription(Qt::Key_V, "Toggles vertices display");
  setKeyDescription(Qt::Key_N, "Inverse direction of normals");
  setKeyDescription(Qt::Key_Plus, "Increase size of edges");
  setKeyDescription(Qt::Key_Minus, "Decrease size of edges");
  setKeyDescription(Qt::Key_Plus+Qt::ShiftModifier, "Increase size of vertices");
  setKeyDescription(Qt::Key_Minus+Qt::ShiftModifier, "Decrease size of vertices");
  setKeyDescription(Qt::Key_PageDown, "Increase light (all colors, use shift/alt/ctrl for one rgb component)");
  setKeyDescription(Qt::Key_PageUp, "Decrease light (all colors, use shift/alt/ctrl for one rgb component)");

  // Light default parameters
  glLineWidth(size_edges);
  glPointSize(size_points);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0f,1.0f);
  glClearColor(1.0f,1.0f,1.0f,0.0f);
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  glEnable(GL_LIGHTING);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glShadeModel(GL_FLAT);
  glDisable(GL_BLEND);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_POLYGON_SMOOTH_HINT);
  glBlendFunc(GL_ONE, GL_ZERO);
  glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
  compile_shaders();
}

void Viewer::keyPressEvent(QKeyEvent *e)
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
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_F) && (modifiers==Qt::NoButton))
  {
    flatShading = !flatShading;
    if (flatShading)
      displayMessage("Flat shading.");
    else
      displayMessage("Gouraud shading.");
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_E) && (modifiers==Qt::NoButton))
  {
    edges = !edges;
    displayMessage(QString("Draw edges=%1.").arg(edges?"true":"false"));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_V) && (modifiers==Qt::NoButton))
  {
    vertices = !vertices;
    displayMessage(QString("Draw vertices=%1.").arg(vertices?"true":"false"));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_N) && (modifiers==Qt::NoButton))
  {
    inverse_normal = !inverse_normal;
    displayMessage(QString("Inverse normal=%1.").arg(inverse_normal?"true":"false"));
    sceneChanged();
  }
  else if ((e->key()==Qt::Key_Plus) && (modifiers==Qt::KeypadModifier))
  {
    size_edges+=.5;
    displayMessage(QString("Size of edges=%1.").arg(size_edges));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_Minus) && (modifiers==Qt::KeypadModifier))
  {
    if (size_edges>.5) size_edges-=.5;
    displayMessage(QString("Size of edges=%1.").arg(size_edges));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_Plus) && (modifiers==(Qt::ShiftModifier|Qt::KeypadModifier)))
  {
    size_points+=.5;
    displayMessage(QString("Size of points=%1.").arg(size_points));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_Minus) && (modifiers==(Qt::ShiftModifier|Qt::KeypadModifier)))
  {
    if (size_points>.5) size_points-=.5;
    displayMessage(QString("Size of points=%1.").arg(size_points));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
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
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
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
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::ShiftModifier))
  {
    ambient.setX(ambient.x()+.1);
    if (ambient.x()>1.) ambient.setX(1.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::AltModifier))
  {
    ambient.setY(ambient.y()+.1);
    if (ambient.y()>1.) ambient.setY(1.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::ControlModifier))
  {
    ambient.setZ(ambient.z()+.1);
    if (ambient.z()>1.) ambient.setZ(1.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::ShiftModifier))
  {
    ambient.setX(ambient.x()-.1);
    if (ambient.x()<0.) ambient.setX(0.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::AltModifier))
  {
    ambient.setY(ambient.y()-.1);
    if (ambient.y()<0.) ambient.setY(0.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::ControlModifier))
  {
    ambient.setZ(ambient.z()-.1);
    if (ambient.z()<0.) ambient.setZ(0.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
#if QGLVIEWER_VERSION >= 0x020700
    update();
#else
    updateGL();
#endif
  }
  else
    QGLViewer::keyPressEvent(e);
}

QString Viewer::helpString() const
{
  QString text("<h2>L C C   V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with "
    "the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the "
    "camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys "
    "(<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several "
    "keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and "
    "restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the "
    "world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to "
    "save a snapshot. ";
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

