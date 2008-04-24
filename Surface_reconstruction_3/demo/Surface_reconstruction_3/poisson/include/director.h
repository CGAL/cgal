// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

#ifndef DIRECTOR_H
#define DIRECTOR_H

#include <math.h>
#include <stdio.h>
#include <iostream>

class Director {
 private:
  float offset[3];
  float vCenter[3];
  float vRadius;
  float vMat[16];
  float tRotation[4];
  float vRotation[4];
  float vCut[4];

  static void qRotation(float x, float y, float xx, float yy, float* q);
  static void qMulti(float* q, float* r);
  static void qMultiVect(float* q, float* v);
  static void qBar(float* r);

  void matTransfo(float* mat, bool avecTransZoom);
 
 public:
  Director();
  ~Director();

  Director(const Director& b);

  // Dans toute la suite, les coordonnées sont données entre
  // -1 et 1. L'orientation des axes est l'orientation
  // mathématique usuelle : x vers le haut, y vers la droite.

  // position du centre avant rotation, intégré à la matrice de
  // transformation
  const float* center();

  const float* getOffset();

  void setOffset(float x, float y, float z);

  void image(float* point);

  void preImageSansTransZoom(float* point);

  // rayon, facteur de zoom, est intégré à la matrice de transformation
  float&       radius();

  // matrice de transformation, c'est le pointeur à donner
  // à OpenGL par glLoadMatrixf avant d'afficher l'objet
  const float* mat();

  const float* matSansTransZoom();

  const float* tempMat();
  
  void tempReset();

  // équation du plan de coupe, solidaire du repère de l'objet,
  // à donner à OpenGL par glClipPlane avant glLoadMatrixf
  const float* cut();

  // hauteur du plan de coupe
  float&       hCut();

  // effectue la rotation autour de O = (0, 0, 0) qui correspond
  // au déplacement de la souris de (x,y) à (xx, yy).
  void rotation(float x, float y, float xx, float yy);

  void rotationObj(float x, float y, float xx, float yy);

  // effectue la rotation du plan de coupe autour de "centre" 
  // qui correspond au déplacement de la souris de (x,y) à (xx, yy).
  void rotationCut(float x, float y, float xx, float yy);

  void rotZtoCut(float* v, float& angle);

  // effectue la translation qui correspond au déplacement
  // à l'écran de (dx, dy, dz). En général dz = 0.
  void translation(float dx, float dy, float dz);

  // remet le centre à la position voulue.
  void recenter(float x, float y, float z);

  // place dans (X,Y,Z) les coordonnées du point sélectionné à la souris
  // en cliquant en (x,y)
  void select(float x, float y, 
	          float& X, float& Y, float& Z,
	          float& XX, float& YY, float& ZZ);
  
  // renvoie la hauteur de la souris sur la boule dans la 
  // direction du plan de coupe
  float height(float x, float y);

};


#endif
