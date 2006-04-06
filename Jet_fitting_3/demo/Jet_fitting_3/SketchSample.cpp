// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#include "SketchSample.h"

SketchSample::SketchSample() {
  highlight = false;
}

SketchSample::~SketchSample() {
}

void SketchSample::buildDisplayList(GLuint surf) {
  static GLfloat ared[4] = {0.8, 0.1, 0.0, 1.0 };
  glNewList(surf, GL_COMPILE);
  glBegin(GL_TRIANGLES);
  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, ared );
  glNormal3d(0., 0., 1.);
  glVertex3d(0., 0., 0.);
  glVertex3d(1., 0., 0.);
  glVertex3d(0., 1., 1.);
  glEnd();
  glEndList();
}

void SketchSample::buildPicking(GLuint pickSurf) {
  int red, green, blue;
  int i = 1;
  glNewList(pickSurf, GL_COMPILE);

  blue = i % 256;
  green = (i/256) % 256;
  red = i/65536;
  glBegin(GL_TRIANGLES);
  glColor3ub((GLubyte) red, (GLubyte) green, (GLubyte) blue);
  glVertex3d(0., 0., 0.);
  glVertex3d(1., 0., 0.);
  glVertex3d(0., 1., 1.);
  glEnd();
  glEndList();
}

bool SketchSample::selectedColor(GLubyte red, GLubyte green, GLubyte blue) {
  highlight = (blue == 1 && red == 0 && green == 0);
  std::cout << highlight << "\n";
  return false;
}

void SketchSample::drawHighlighted() {
  if(highlight) {
    glBegin(GL_TRIANGLES);
    glColor3d(1., 1., 0.);
    glVertex3d(0., 0., 0.);
    glVertex3d(1., 0., 0.);
    glVertex3d(0., 1., 1.);
    glEnd();
  }

}

double SketchSample::rcx() {
  return 0.;
}

double SketchSample::rcy() {
  return 0.;
}

double SketchSample::rcz() {
  return 0.;
}

double SketchSample::rmm() {
  return 2.;
}

const double* SketchSample::rcoord() {
  return 0;
}
