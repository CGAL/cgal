// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#include "Sketcher.h"
#include <iostream>

//marc
#include "visu_poly.h"
#include "enriched_polyhedron.h" 

enum Ridge_type {NONE=0, BLUE_RIDGE, RED_RIDGE, CREST, BE, BH, BC, RE, RH, RC};

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

typedef std::list<data_line* > DS;
typedef std::list<data_line* >::iterator DS_iterator;


class SketchSample : public Sketcher {
 private: 
  bool highlight;
  Mesh* p_mesh;
  DS* p_ridge_data;
  
 public: 
  SketchSample(Mesh* mesh, DS* ridge_data);  
  
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
};
