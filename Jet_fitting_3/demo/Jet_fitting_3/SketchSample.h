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


//Data structure for one line of the input file .4ogl.txt
struct data_line{
  Point P1,P2;
  Vector D1,D2;
  double k1,k2;
  data_line(Point P1,Point P2,
	    Vector  D1,Vector D2, double k1,double k2 ):
    P1(P1),P2(P2),D1(D1),D2(D2),k1(k1),k2(k2)
    {};
};

typedef std::list<data_line* > DS_;
typedef std::list<data_line* >::iterator DS_iterator;



class SketchSample : public Sketcher {
 private:
  bool highlight;
  Mesh* p_mesh;
  DS_* p_ppal_data;

 public:
  SketchSample(Mesh* mesh, DS_* ppal_data);

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
  void SketchSample::draw_point(Point& P);
  void SketchSample::draw_vector(Point& P, Vector& V);

};
