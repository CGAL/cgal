/***************************************************************************

    begin                : jan 02
    copyright            : (C) 2002 by Pierre Alliez
    email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "glview.h"
#include <stdio.h>
#include <assert.h>
#include "toolview.h"

QGLView::QGLView(const QGLFormat & format,
                 QWidget *parent,
                 const char *name,
                 int wflags)
    : QGLWidget(format, parent, name, 0, wflags)
{
  m_Culling = false;
  m_Antialiasing = false;
  m_Lighting = false;
  m_Texture = false;
  m_Smooth = false;
  m_Superimpose_edges = false;
  m_Superimpose_vertices = false;
  m_clear_color[0] = 0.0f;
  m_clear_color[1] = 0.0f;
  m_clear_color[2] = 0.0f;
  m_fore_color[0] = 1.0f;
  m_fore_color[1] = 1.0f;
  m_fore_color[2] = 1.0f;
  m_mode = FILL;
}

QGLView::~QGLView()
{
}

void QGLView::set_clear_color(float r,
                              float g,
                              float b)
{
  makeCurrent();
  glClearColor(r,g,b,0.0f);
  m_clear_color[0] = r;
  m_clear_color[1] = g;
  m_clear_color[2] = b;
}

void QGLView::set_fore_color(float r,
                             float g,
                             float b)
{
  makeCurrent();
  glColor3f(r,g,b);
  m_fore_color[0] = r;
  m_fore_color[1] = g;
  m_fore_color[2] = b;
}

void QGLView::initializeGL()
{
  glClearColor(0.0f,0.0f,0.0f,0.0f);
  glShadeModel(GL_FLAT);
  glPolygonMode(GL_FRONT,GL_LINE);
  glPolygonMode(GL_BACK,GL_LINE);

  // Texturing
  glEnable(GL_TEXTURE_2D);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
}

void QGLView::clearColor(float r, float g, float b)
{
  makeCurrent();
  glClearColor(r,g,b,0.0f);
}

void QGLView::changeMode(int mode)
{
  makeCurrent();
  glPolygonMode(GL_FRONT,mode);
  glPolygonMode(GL_BACK,mode);
}

void QGLView::toggleCulling()
{
  makeCurrent();
  m_Culling = !m_Culling;
  if(m_Culling)
    glEnable(GL_CULL_FACE);
  else
    glDisable(GL_CULL_FACE);
}

void QGLView::toggleSmooth()
{
  makeCurrent();
  m_Smooth = !m_Smooth;
  if(m_Smooth)
    glShadeModel(GL_SMOOTH);
  else
    glShadeModel(GL_FLAT);
}

void QGLView::enableTexture(bool enable)
{
  makeCurrent();
  m_Texture = enable;
  if(m_Texture)
  {
    glEnable(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NICEST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NICEST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  }
  else
    glDisable(GL_TEXTURE_2D);
}

void QGLView::toggleLighting()
{
  makeCurrent();
  m_Lighting = !m_Lighting;
  if(m_Lighting)
  {
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    float lpos[4] = { -.2f, .2f, .9797958971f, 0 };
    glLightfv(GL_LIGHT0,GL_POSITION,lpos);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

    glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
  }
  else
  {
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
  }
}

void QGLView::toggleSuperimposing()
{
  makeCurrent();
  m_Superimpose_edges = !m_Superimpose_edges;
}

void QGLView::toggleSuperimposeV()
{
  makeCurrent();
  m_Superimpose_vertices = !m_Superimpose_vertices;
}

void QGLView::toggleMode(Rendering_mode m)
{
  makeCurrent();
  if(m == FILL){
    changeMode(GL_FILL);
    m_mode = FILL;
  } else if(m == LINE){
    changeMode(GL_LINE);
    m_mode = LINE;
  } else {
    changeMode(GL_POINT);
    m_mode = VERTEX;
  }
}
void QGLView::toggleAntialiasing()
{
  makeCurrent();
  m_Antialiasing = !m_Antialiasing;
  if(m_Antialiasing)
  {
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
    glLineWidth(1.5);
  }
  else
  {
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);
    glLineWidth(1.0);
  }
}

void QGLView::paintGL()
{
  glClearColor(m_clear_color[0], m_clear_color[1], m_clear_color[2], 0.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();

  // set up the current transform
  theArcball.glDraw();

  // paint scene
  ToolView *pView = (ToolView *)this;

  if(m_Superimpose_edges || m_Superimpose_vertices)
  {
    if(m_Superimpose_edges && m_Superimpose_vertices)
      pView->drawScene(true, true, true);
    else
      if(m_Superimpose_edges)
        pView->drawScene(true, false, true);
      else
        pView->drawScene(false, true, true);

    bool lighting = ::glIsEnabled(GL_LIGHTING);
    if(lighting)
    {
      glDisable(GL_LIGHTING);
      glDisable(GL_COLOR_MATERIAL);
    }
    
    if(m_Superimpose_edges && m_Superimpose_vertices)
      pView->drawScene(true, true, false);
    else
      if(m_Superimpose_edges)
        pView->drawScene(true, false, false);
      else
        pView->drawScene(false, true, false);
    
    if(lighting)
    {
      glEnable(GL_LIGHTING);
      glEnable(GL_COLOR_MATERIAL);
    }
  }
  else
  {
    pView->drawScene(false, false);
  }

  glPopMatrix();
}

void QGLView::resizeGL( int w, int h )
{
  // set aspect
   double width = (double)w;
   double height = (double)h;
   double aspect = width;
   if(height != 0)
     aspect = width/height;

  // Set OpenGL viewport and perspective
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  //fprintf(stderr,"update viewport...\n");
  theViewport.SetOrigin(0,0);
  theViewport.SetSize(w,h);
  theCamera.SetHeightAngle(45.);
  theCamera.SetPosition(0.,0.,5.);
  theCamera.SetOrientation(0.,1.,0.,0.);
  theCamera.SetNearDistance(0.1f);
  theCamera.SetFarDistance(1000.0);

  theViewport.glDraw();
  theCamera.glDraw(&theViewport);
}

GLuint QGLView::makeObject()
{
  GLuint list = glGenLists( 1 );
  glNewList(list, GL_COMPILE );
  qglColor(white);

  float x = 2.0f;

  glBegin(GL_POLYGON);
    glVertex3d( x,  x, x);
    glVertex3d( x, -x, x);
    glVertex3d(-x, -x, x);
    glVertex3d(-x,  x, x);
  glEnd();

  glBegin(GL_POLYGON);
    glVertex3d(-x,  x, -x);
    glVertex3d(-x, -x, -x);
    glVertex3d( x, -x, -x);
    glVertex3d( x,  x, -x);
  glEnd();

  glBegin(GL_POLYGON);
    glVertex3d( x,  x,  x);
    glVertex3d( x,  x, -x);
    glVertex3d( x, -x, -x);
    glVertex3d( x, -x,  x);
  glEnd();

  glBegin(GL_POLYGON);
    glVertex3d(-x,  x,  x);
    glVertex3d(-x,  x, -x);
    glVertex3d(-x, -x, -x);
    glVertex3d(-x, -x,  x);
  glEnd();

  glBegin(GL_POLYGON);
    glVertex3d(-x, -x,  x);
    glVertex3d( x, -x,  x);
    glVertex3d( x, -x, -x);
    glVertex3d(-x, -x, -x);
  glEnd();


  glBegin(GL_POLYGON);
    glVertex3d(-x,  x,  x);
    glVertex3d( x,  x,  x);
    glVertex3d( x,  x, -x);
    glVertex3d(-x,  x, -x);
  glEnd();

  glEndList();

  return list;
}

#include "glview.moc"