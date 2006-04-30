// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#include "SketchSample.h"

SketchSample::SketchSample(Mesh* mesh, DS* ridge_data) {
  highlight = false;
  p_mesh = mesh;
  p_ridge_data = ridge_data;
}

SketchSample::~SketchSample() {
}

void SketchSample::buildDisplayList(GLuint surf) {

  glNewList(surf, GL_COMPILE);
  glDisable(GL_LIGHTING);

  for (DS_iterator it=p_ridge_data->begin();it!=p_ridge_data->end();it++)
    draw_one_ridge(*it);

  glEnable(GL_LIGHTING);
  glShadeModel(GL_SMOOTH);

  //mesh
  //glMaterialfv is def in the mesh 
  p_mesh->gl_draw_facets(true);

  //additional objects that may be displayed
//   p_mesh->gl_draw_vertices_normal(); 
//   p_mesh->gl_draw_edges();
//   p_mesh->gl_draw_vertices();
//   p_mesh->gl_draw_facets_normal();
//   p_mesh->gl_draw_vertices_normal();
 
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


void SketchSample::draw_one_ridge(data_line* line)
{
  if (line->ridge_type == BE) glColor3f(0.,0.,1.);
  if (line->ridge_type == BH) glColor3f(0.,1.,0.);
  if (line->ridge_type == BC) glColor3f(0.,0.,1.);
  if (line->ridge_type == RE) glColor3f(1.,0.,0.);
  if (line->ridge_type == RH) glColor3f(1.,1.,0.);
  if (line->ridge_type == RC) glColor3f(1.,0.,0.);
  
  std::list<Point>::iterator iter = line->ridge_points.begin(), 
    ite = line->ridge_points.end();

  glLineWidth(3 );
  glBegin(GL_LINES);
  for (;iter!=ite;iter++) glVertex3d(iter->x(), iter->y(), iter->z());
  glEnd();	
}
