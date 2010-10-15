// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#include "SketchSample.h"

SketchSample::SketchSample(Mesh* mesh, DS_* ppal_data) {
  highlight = false;
  p_mesh = mesh;
  p_ppal_data = ppal_data;
}

SketchSample::~SketchSample() {
}

void SketchSample::buildDisplayList(GLuint surf) {

  glNewList(surf, GL_COMPILE);

  glDisable(GL_LIGHTING);

  //points first
  for (DS_iterator it=p_ppal_data->begin();it!=p_ppal_data->end();it++)
    draw_point((*it)->P1);
  //next, ppal dirs and normal
  for (DS_iterator it=p_ppal_data->begin();it!=p_ppal_data->end();it++)
    {
      glColor3f(0.,0.,1.);//dmax
      draw_vector((*it)->P1,(*it)->D1);
      glColor3f(1.0,0.0,0.);//dmin
      draw_vector((*it)->P1,(*it)->D2);
      glColor3f(0.0,1.,.0);//normal
      Vector normal = CGAL::cross_product( (*it)->D1 ,(*it)->D2 )
	/ CGAL::sqrt( ((*it)->D1) * ((*it)->D1) );
      draw_vector((*it)->P1, normal) ;
    }

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

// //move to sketch
// void draw_point(Point& P);
// void draw_vector(Point& P, Vector& V);
// void MakeCallList(DS_& L);


void SketchSample::draw_point(Point& P)
{
  glPointSize(2.0);
  glBegin(GL_POINTS);
  glColor3f(1.0,1.,1.);
  glVertex3d(P.x(),P.y(),P.z());
  glEnd();
  glPointSize(1.0);
}

void SketchSample::draw_vector(Point& P, Vector& V)
{
  glBegin(GL_LINES);
  glVertex3d(P.x()-V.x()/2.,P.y()-V.y()/2.,P.z()-V.z()/2.);
  glVertex3d(P.x()+V.x()/2.,P.y()+V.y()/2.,P.z()+V.z()/2.);
  glEnd();

  glPointSize(3.0);
  glBegin(GL_POINTS);
  glVertex3d(P.x()+V.x()/2.,P.y()+V.y()/2.,P.z()+V.z()/2.);
  glEnd();
  glPointSize(1.0);
}
