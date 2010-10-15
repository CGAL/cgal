// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#include "SketchSample.h"
#include <CGAL/bounding_box.h>

extern double strength_threshold;
extern double sharpness_threshold;

SketchSample::SketchSample(Mesh* mesh, DS_* ridge_data) {
  highlight = false;
  p_mesh = mesh;
  p_ridge_data = ridge_data;
  this->compute_mesh_position();
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
  if (line->ridge_type == CGAL::MAX_ELLIPTIC_RIDGE)   glColor3f(0.,0.,1.);
  if (line->ridge_type == CGAL::MAX_HYPERBOLIC_RIDGE) glColor3f(0.,1.,0.);
  if (line->ridge_type == CGAL::MAX_CREST_RIDGE)            glColor3f(0.,0.,1.);
  if (line->ridge_type == CGAL::MIN_ELLIPTIC_RIDGE)    glColor3f(1.,0.,0.);
  if (line->ridge_type == CGAL::MIN_HYPERBOLIC_RIDGE)  glColor3f(1.,1.,0.);
  if (line->ridge_type == CGAL::MIN_CREST_RIDGE)             glColor3f(1.,0.,0.);

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
  return mesh_center_x;
}

double SketchSample::rcy() {
  return mesh_center_y;
}

double SketchSample::rcz() {
  return mesh_center_z;
}

double SketchSample::rmm() {
  return mesh_radius;
}

const double* SketchSample::rcoord() {
  return 0;
}

void SketchSample::compute_mesh_position() {
  typedef CGAL::Cartesian<double> Kernel;
  Kernel::Iso_cuboid_3  box = CGAL::bounding_box(p_mesh->points_begin(),
					   p_mesh->points_end());
  Kernel::Point_3 center = CGAL::midpoint(box.min(), box.max());
  Kernel::Vector_3 diameter = box.min() - box.max();
  mesh_center_x = center.x();
  mesh_center_y = center.y();
  mesh_center_z = center.z();
  mesh_radius = CGAL::sqrt( diameter * diameter )/2;
}
