// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#include "Sketcher.h"
#include <iostream>

class SketchSample : public Sketcher {
 private: 
  bool highlight;

 public: 
  SketchSample();

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
};
