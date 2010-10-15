// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#include "Sketcher.h"
#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Ridges.h>
//marc
#include "visu_poly.h"
#include "enriched_polyhedron.h"

/* extern CGAL::Ridge_type  NONE, BLUE_RIDGE, RED_RIDGE, CREST,  */
/*   BLUE_ELLIPTIC_RIDGE, BLUE_HYPERBOLIC_RIDGE, BLUE_CREST,  */
/*   RED_ELLIPTIC_RIDGE, RED_HYPERBOLIC_RIDGE, RED_CREST; */

//Data structure for one line of the input file .4ogl.txt
struct data_line{
  int ridge_type;
  double strength, sharpness;
  std::vector<Point> ridge_points;
  data_line(int ridge_type, double strength, double sharpness,
	    std::vector<Point> ridge_points):
  ridge_type(ridge_type), strength(strength), sharpness(sharpness),
    ridge_points(ridge_points)
    {};
};

typedef std::list<data_line* > DS_;
typedef std::list<data_line* >::iterator DS_iterator;


class SketchSample : public Sketcher {
 private:
  bool highlight;
  Mesh* p_mesh;
  DS_* p_ridge_data;
  double mesh_center_x, mesh_center_y, mesh_center_z, mesh_radius;

 public:
  SketchSample(Mesh* mesh, DS_* ridge_data);

  ~SketchSample();

  void buildDisplayList(GLuint surf);
  void buildPicking(GLuint pickSurf);
  bool selectedColor(GLubyte red, GLubyte green, GLubyte blue);
  void drawHighlighted();
  double rcx();
  double rcy();
  double rcz();
  double rmm();
  const double* rcoord();
  void draw_one_ridge(data_line* line);
  void compute_mesh_position();
};
