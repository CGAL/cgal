// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
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
//
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#include <CGAL/basic.h>

#include <CGAL/IO/Qt_widget_OpenGL.h>
#include <cmath>

namespace CGAL {

Qt_widget_OpenGL::Qt_widget_OpenGL(int width, int height, double scale) :
  window_width(width),
  window_height(height),
  motion_mode(ROTATE),
  dx(0),
  dy(0),
  dz(0),
  s(scale),
  init_s(scale),
  rotation(CGAL::IDENTITY){}
  
CGAL::OGL::OGL_base_object::Affine_3
Qt_widget_OpenGL::virtual_sphere_transformation( double old_x, double old_y, 
						 double new_x, double new_y) {

  if ( old_x == new_x && old_y == new_y)// zero rotation.
	return Affine_3( CGAL::IDENTITY);
    // Determine the projected vectors on the `sphere'.
    double dd = old_x * old_x + old_y * old_y;
    Double_vector v_old( old_x, old_y, 
                         ((dd < 0.5) ? std::sqrt(1-dd) : 0.5 / std::sqrt(dd)));
    dd = new_x * new_x + new_y * new_y;
    Double_vector v_new( new_x, new_y, 
                         ((dd < 0.5) ? std::sqrt(1-dd) : 0.5 / std::sqrt(dd)));
    Double_vector axis  = cross_product( v_old, v_new);
    double angle = 0.0;
    double norm = std::sqrt( (v_old*v_old)*(v_new*v_new));
    if ( norm != 0) {
        double x = v_old*v_new/ norm;
        if ( x <= -1)
        angle = CGAL_PI;
        if ( x < 1)
            angle = std::acos(x);
    }
    double len = std::sqrt( double(axis * axis));
    double s   = std::sin( angle / 2.0) / len;
    double q1 = axis.x() * s; // quaternion
    double q2 = axis.y() * s;
    double q3 = axis.z() * s;
    double q0 = std::cos( angle / 2.0);
    double a   = q1 * q2;
    double b   = q0 * q3;
    double c   = q1 * q3;
    double d   = q0 * q2;
    double e   = q2 * q3;
    double f   = q0 * q1;
    double qq0 = q0 * q0;
    double qq1 = q1 * q1;
    double qq2 = q2 * q2;
    double qq3 = q3 * q3;
    return Affine_3( qq0 + qq1 - qq2 - qq3, 2 * (a-b), 2 * (c+d),
		     2 * (a+b), qq0 - qq1 + qq2 - qq3, 2 * (e-f),
		     2 * (c-d), 2 * (e+f), qq0 - qq1 - qq2 + qq3);
}

void Qt_widget_OpenGL::mouseMoveEvent(QMouseEvent* event) {
  int x = event->x();
  int y = event->y();
  switch ( interaction) {
  case SCALE:
    s *= std::exp( (x - mouse_x + mouse_y -y) * factor_s );
    break;
  case ROTATE: {
    double old_x =   1.2 * (mouse_x -  window_width/2) / window_radius;
    double old_y = - 1.2 * (mouse_y - window_height/2) / window_radius;
    double new_x =   1.2 * (x -  window_width/2) / window_radius;
    double new_y = - 1.2 * (y - window_height/2) / window_radius;
    rotation =  virtual_sphere_transformation( old_x, old_y, new_x, new_y)
      * rotation;
  }
    break;
  case TRANSLATE:
    dx += (x - mouse_x) * 2.0 / window_radius;
    dy -= (y - mouse_y) * 2.0 / window_radius;
    break;
  case TRANS_Z:
    dz += (x - mouse_x + mouse_y -y) * 2.0 / window_radius;
    break;
  default:
    break;
  }
  mouse_x = x;
  mouse_y = y;

  updateGL();
}

void Qt_widget_OpenGL::mousePressEvent(QMouseEvent* event) {
  mouse_x = event->x();
  mouse_y = event->y();
  interaction = 0;
  if (event->stateAfter() & QMouseEvent::LeftButton) {
    if (event->stateAfter() & QMouseEvent::ShiftButton)
      interaction = SCALE;
    else
      interaction = motion_mode;
  }
  if(event->stateAfter() & QMouseEvent::MidButton) {
    if (event->stateAfter() & QMouseEvent::ShiftButton)
      interaction = TRANS_Z;
    else
      interaction = TRANSLATE;
  }
  if(event->stateAfter() & QMouseEvent::RightButton)
    main->exec(QPoint(event->globalX(),event->globalY()));
}

void Qt_widget_OpenGL::mouseReleaseEvent(QMouseEvent* event) {
  mousePressEvent(event);
}

void Qt_widget_OpenGL::paintGL() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glPushMatrix();
  glTranslated(dx,dy,dz);
  //  glTranslated(0,0,1);
  GLdouble M[16] = { rotation.m(0,0), rotation.m(1,0), rotation.m(2,0), 0.0,
                     rotation.m(0,1), rotation.m(1,1), rotation.m(2,1), 0.0,
                     rotation.m(0,2), rotation.m(1,2), rotation.m(2,2), 0.0,
                     rotation.m(0,3), rotation.m(1,3), rotation.m(2,3), 1.0};
  glMultMatrixd( M);
  glScaled(s,s,s);
  object_->draw();
  glPopMatrix();
}

void Qt_widget_OpenGL::initializeGL() {
  GLfloat mat_diffuse[4] = { 0.7f, 0.7f, 0.7f, 1.0f };
  GLfloat mat_specular[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
  GLfloat mat_shininess[] = { 100.0f };
  GLfloat ambient_light[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse );
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
#ifdef SCREENSHOTS
  GLfloat mat_emission[] = { 0.1, 0.1, 0.2, 0.0 };
  glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
  //for screenshots enable this section
#endif

  GLfloat light0[4] = { 4.0, 4.0, 10.0, 1.0 };
  glLightfv (GL_LIGHT0, GL_POSITION, light0);
  glEnable (GL_LIGHT0);
 
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);
}

void Qt_widget_OpenGL::resizeGL(int width, int height) {
  window_width = width;
  window_height = height;
  window_radius = (std::min)( width, height) / 2;

  glViewport(0, 0, (GLint)width, (GLint)height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (width>height)
    {
      long double w = (long double) width / (long double) height;
      glOrtho( -2*w, 2*w, -w, w, -4.0, 4.0 );
      factor_d =  2.0   / (height/2.0);           
      // halbe Fensterhoehe soll 2 LE entsprechen
      factor_w = 90.0   / (height/2.0);           
      // halbe Fensterhoehe soll 90 Grad entsprechen
      factor_s = std::log(4.0) / (height/2.0);           
      // halbe Fensterhoehe soll Faktor 4 entsprechen
    }  
  else
    {
      long double h = (long double) height / (long double) width;
      glOrtho( -2.0, 2.0, -2*h, 2*h, -4.0, 4.0 );
      factor_d =  2.0   / (width/2.0);            
      // halbe Fensterbreite soll 2 LE entsprechen
      factor_w = 90.0   / (width/2.0);            
      // halbe Fensterbreite soll 90 Grad entsprechen
      factor_s = std::log(4.0) / (height/2.0);           
      // halbe Fensterhoehe soll Faktor 4 entsprechen
    }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void Qt_widget_OpenGL::slotControlMenu(int index) {
  if(index == RESET_CONTROL) {
    dx = dy = dz = 0.0;
    s = init_s;
    rotation = Affine_3(CGAL::IDENTITY);
    motion_mode = ROTATE;
    object_->init();
    updateGL();
  } else 
    motion_mode = index;
}

void Qt_widget_OpenGL::slotRenderMenu(int index) {
  object_->set_style(index);
  updateGL();
}

void Qt_widget_OpenGL::slotOptionsMenu(int index) {
  object_->toggle(index);
  updateGL();
}

void Qt_widget_OpenGL::slotFullscreen() {
  if(fullscreen) {
    showNormal();
  } else {
    showMaximized();
  }
  fullscreen = !fullscreen;
}

void Qt_widget_OpenGL::slotPerspective() {
  perspective = !perspective;
}

} // namespace CGAL
#include "Qt_widget_OpenGL.moc"
