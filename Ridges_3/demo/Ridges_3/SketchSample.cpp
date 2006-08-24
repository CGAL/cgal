// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#include "SketchSample.h"

extern double strength_threshold;
extern double sharpness_threshold;

SketchSample::SketchSample(Mesh* mesh, DS* ridge_data) {
  highlight = false;
  p_mesh = mesh;
  p_ridge_data = ridge_data;
}

SketchSample::~SketchSample() {
}

void SketchSample::buildDisplayList(GLuint surf) {
  glNewList(surf, GL_COMPILE);
  
  glEnable(GL_LIGHTING);
  glShadeModel(GL_SMOOTH);
  //glShadeModel(GL_FLAT);
  //mesh glMaterialfv is def in the mesh with offset to see the lines z_buffered
  p_mesh->gl_draw_facets(true);

  //ridges, drawn without light
  glDisable(GL_LIGHTING);
  glColor3d(1., 1., 0.);

  for (DS_iterator it=p_ridge_data->begin();it!=p_ridge_data->end();it++)
    {
      if ( (*it)->strength >= strength_threshold && 
	   (*it)->sharpness >= sharpness_threshold )
	draw_one_ridge(*it);
     }
  glEnable(GL_LIGHTING);

  //additional objects that may be displayed
//   p_mesh->gl_draw_vertices_normal(); 
//   p_mesh->gl_draw_edges();
//   p_mesh->gl_draw_vertices();
//   p_mesh->gl_draw_facets_normal();
//   p_mesh->gl_draw_vertices_normal();
 
  glEndList();
}

void SketchSample::draw_one_ridge(data_line* line)
{
  if (line->ridge_type == CGAL::BLUE_ELLIPTIC_RIDGE)   glColor3f(0.,0.,1.);
  if (line->ridge_type == CGAL::BLUE_HYPERBOLIC_RIDGE) glColor3f(0.,1.,0.);
  if (line->ridge_type == CGAL::BLUE_CREST)            glColor3f(0.,0.,1.);
  if (line->ridge_type == CGAL::RED_ELLIPTIC_RIDGE)    glColor3f(1.,0.,0.);
  if (line->ridge_type == CGAL::RED_HYPERBOLIC_RIDGE)  glColor3f(1.,1.,0.);
  if (line->ridge_type == CGAL::RED_CREST)             glColor3f(1.,0.,0.);
  
  std::vector<Point>::iterator iter = line->ridge_points.begin(), 
    ite = line->ridge_points.end();

  glLineWidth(3 );
  glBegin(GL_LINE_STRIP);
  for (;iter!=ite;iter++)  glVertex3d(iter->x(), iter->y(), iter->z()); 
  glEnd();	
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


