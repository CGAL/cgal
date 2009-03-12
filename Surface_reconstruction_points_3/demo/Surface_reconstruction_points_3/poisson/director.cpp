// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// Author: Camille Wormser
// Version: 1.0

//#include <stdafx.h>

#include <algorithm>
#include "Director.h"

Director::Director() {
  offset[0] = offset[1] = offset[2] = 0.;
  vCenter[0] = vCenter[1] = vCenter[2] = 0.;
  vRadius = 1.;
  vRotation[0] = 1.;
  vRotation[1] = vRotation[2] = vRotation[3] = 0.;
  tRotation[0] = 1.;
  tRotation[1] = tRotation[2] = tRotation[3] = 0.;
  vCut[0] = 1;
  vCut[1] = vCut[2] = 0.;
  vCut[3] = 0.1f;
}

Director::~Director() {
}

Director::Director(const Director& b)
{
  std::copy(b.vCenter, b.vCenter+3, vCenter);

  vRadius = b.vRadius;

  std::copy(b.vRotation, b.vRotation+4, vRotation);

  std::copy(b.vCut, b.vCut+4, vCut);

  std::copy(b.vMat, b.vMat+16, vMat);
}

// calcule le quaternion de la rotation définie par le glissement de souris
// (x,y) -> (xx, yy) 
void Director::qRotation(float x, float y, float xx, float yy, float* q) {
  x = x;
  y = y;
  xx = xx;
  yy = yy;
  float norm = sqrt(x*x + y*y);
  if(norm > .9999f) {
    x = .9999f*x/norm;
    y = .9999f*y/norm;
  }
  norm = sqrt(xx*xx + yy*yy);
  if(norm > .9999f) {
    xx = .9999f*xx/norm;
    yy = .9999f*yy/norm;
  }
  // std::cout << 1. - x*x - y*y << " xx\n";
  float z = sqrt(1.f - x*x - y*y);
  float zz = sqrt(1.f - xx*xx - yy*yy);
  // calcule le vecteur de rotation
  float rx = y*zz - z*yy;
  float ry = z*xx - x*zz;
  float rz = x*yy - y*xx;
  float n = sqrt(rx*rx + ry*ry + rz*rz);
  rx = rx/2.f;
  ry = ry/2.f;
  rz = rz/2.f;
  // en pratique, angle petit :
  // et n = sin t = t
  float cos = (float)sqrt(1. - n*n/4.);
  //  std::cout << "cos t/2 = " << cos << "  et  sin t/2 v = " << rx << ";"<< ry << ";"<< rz << "\n";
  // calcule le quaternion
  q[0] = cos;
  q[1] = rx;
  q[2] = ry;
  q[3] = rz;
}

// effectue r = q*r
void Director::qMulti(float* q, float* r) {
  float rx = q[0]*r[0] - q[1]*r[1] - q[2]*r[2] - q[3]*r[3];
  float ry = q[0]*r[1] + q[1]*r[0] + q[2]*r[3] - q[3]*r[2];
  float rz = q[0]*r[2] + q[2]*r[0] + q[3]*r[1] - q[1]*r[3];
  float rt = q[0]*r[3] + q[3]*r[0] + q[1]*r[2] - q[2]*r[1];
  r[0] = rx;
  r[1] = ry;
  r[2] = rz;
  r[3] = rt;
}

// effectue r = bar(r)
void Director::qBar(float* r) {
  r[1] = -r[1];
  r[2] = -r[2];
  r[3] = -r[3];
}

// effectue v = Rot(q)*v où v est un vecteur
// en calculant q*v*bar(q)
// attention, si q n'est pas unitaire, il faudra
// renormaliser d'un facteur |q|^2
void Director::qMultiVect(float* q, float* v) {
  float w[] = {0., v[0], v[1], v[2]};
  qBar(w);
  qMulti(q,w);
  qBar(w);
  qMulti(q,w);
  v[0] = w[1];
  v[1] = w[2];
  v[2] = w[3];
}

void Director::matTransfo(float* Rot, bool avecTransZoom) {
  float xx = Rot[1]*Rot[1];
  float xy = Rot[1]*Rot[2];
  float xz = Rot[1]*Rot[3];
  float xw = Rot[1]*Rot[0];
  float yy = Rot[2]*Rot[2];
  float yz = Rot[2]*Rot[3];
  float yw = Rot[2]*Rot[0];
  float zz = Rot[3]*Rot[3];
  float zw = Rot[3]*Rot[0];
  if(avecTransZoom) {
  //  std::cout << "fin précalc\n";
  vMat[0]  = vRadius*(1 - 2*(yy + zz));
  vMat[1]  = vRadius*     2*(xy + zw);
  vMat[2]  = vRadius*     2*(xz - yw);
  vMat[4]  = vRadius*     2*(xy - zw);
  vMat[5]  = vRadius*(1 - 2*(xx + zz));
  vMat[6]  = vRadius*     2*(yz + xw);
  vMat[8]  = vRadius*     2*(xz + yw);
  vMat[9]  = vRadius*     2*(yz - xw);
  vMat[10] = vRadius*(1 - 2*(xx + yy));
  vMat[3]  = vMat[7] = vMat[11] = 0;

  vMat[12] = vRadius*vCenter[0];
  vMat[13] = vRadius*vCenter[1];
  vMat[14] = vRadius*vCenter[2];

  // si on veut zoomer par rapport au centre de l'objet
  //  vMat[12] = vCenter[0];
  //  vMat[13] = vCenter[1];
  //  vMat[14] = vCenter[2];
  vMat[15] = 1;
  }
  else {
  vMat[0]  = (1 - 2*(yy + zz));
  vMat[1]  =      2*(xy + zw);
  vMat[2]  =      2*(xz - yw);
  vMat[4]  =      2*(xy - zw);
  vMat[5]  = (1 - 2*(xx + zz));
  vMat[6]  =      2*(yz + xw);
  vMat[8]  =      2*(xz + yw);
  vMat[9]  =      2*(yz - xw);
  vMat[10] = (1 - 2*(xx + yy));
  vMat[3]  = vMat[7] = vMat[11] = 0;
  vMat[12] = 0.;
  vMat[13] = 0.;
  vMat[14] = 0.;
  vMat[15] = 1;
  }
}

void Director::rotZtoCut(float* v, float& angle) {
  float n = 1 - vCut[2]*vCut[2];
  if(n <= 0.0001) {
    v[0] = 1.;
    v[1] = v[2] = 0.;
    if(vCut[2] > 0.) 
      angle = 0.;
    else 
      angle = 3.1416f;
  }
  else {
    n = sqrt(n);
    v[0] = -vCut[1]/n;
    v[1] = vCut[0]/n;
    v[2] = 0.;
    if(vCut[2] > 0.)
      angle = asin(n);
    else 
      angle = 3.1416f - asin(n);
  }
}

const float* Director::center() {
  return vCenter;
}

const float* Director::getOffset() {
  return offset;
}

void Director::preImageSansTransZoom(float* point) {
  qBar(vRotation);
  qMultiVect(vRotation, point);
  qBar(vRotation);
}

void Director::image(float* point) {
  //  std::cout << "point: " << point[0] << " : " << point[1] << " : " << point[2]<< "\n";
  //  std::cout << "offst: " << offset[0] << " : " << offset[1] << " : " << offset[2]<< "\n";
  point[0] += offset[0];
  point[1] += offset[1];
  point[2] += offset[2];
  qMultiVect(vRotation, point);
  //  std::cout << "Rpont: " << point[0] << " : " << point[1] << " : " << point[2]<< "\n";
  point[0] += vCenter[0];
  point[1] += vCenter[1];
  point[2] += vCenter[2];
  //  std::cout << "Tpont: " << point[0] << " : " << point[1] << " : " << point[2]<< "\n";
  point[0] *= vRadius;
  point[1] *= vRadius;
  point[2] *= vRadius;
  //  std::cout << "Zpont: " << point[0] << " : " << point[1] << " : " << point[2]<< "\n";
}

float& Director::radius() {
  return vRadius;
}

const float* Director::mat() {
  matTransfo(vRotation, true);
  return vMat;
}

const float* Director::matSansTransZoom() {
  matTransfo(vRotation, false);
  return vMat;
}

const float* Director::tempMat() {
  matTransfo(tRotation, false);
  return vMat;
}

void Director::tempReset() {
  tRotation[0] = 1.;
  tRotation[1] = tRotation[2] = tRotation[3] = 0.;
}

const float* Director::cut() {
  return vCut;
}

float& Director::hCut() {
  return vCut[3];
}

void Director::rotation(float x, float y, float xx, float yy) {
  float q[] = {0., 0., 0., 0.};
  qRotation(x, y, xx, yy, q);
  qMulti(q, vRotation);
  qMultiVect(q, vCenter);
  qMulti(q, tRotation);
}

void Director::rotationObj(float x, float y, float xx, float yy) {
  float q[] = {0., 0., 0., 0.};
  qRotation(x, y, xx, yy, q);
  qMulti(q, vRotation);
}

void Director::rotationCut(float x, float y, float xx, float yy) {
  float q[] = {0., 0., 0., 0.};
  qRotation(x, y, xx, yy, q);
  qBar(q);
  qBar(vRotation);
  qMulti(vRotation, q);
  qBar(q);
  qMulti(vRotation, q);
  qBar(vRotation);
  qMultiVect(q, vCut);
}

void Director::setOffset(float x, float y, float z) {
  offset[0] = -x;
  offset[1] = -y;
  offset[2] = -z;
}

void Director::recenter(float x, float y, float z) {
  vCenter[0] = -x;
  vCenter[1] = -y;
  vCenter[2] = -z;

  qMultiVect(vRotation, vCenter);
  //  vCenter[0] *= vRadius;
  //  vCenter[1] *= vRadius;
  //  vCenter[2] *= vRadius;
  
}

void Director::translation(float dx, float dy, float dz) {
  vCenter[0] += dx/vRadius;
  vCenter[1] += dy/vRadius;
  vCenter[2] += dz/vRadius;

  // si on veut zoomer par rapport au centre de l'objet
  //  vCenter[0] += dx;
  //  vCenter[1] += dy;
  //  vCenter[2] += dz;
}

void Director::select(float x, float y, 
                      float& X, float& Y, float& Z,
                      float& XX, float& YY, float& ZZ) {
  float norm = sqrt(x*x + y*y);
  if(norm > .9999f) {
    x = .9999f*x/norm;
    y = .9999f*y/norm;
  }
  float w[] = {x, y, (float)sqrt(1. - x*x - y*y)};
  w[0] -= vCenter[0];
  w[1] -= vCenter[1];
  w[2] -= vCenter[2];
  w[0] /= vRadius;
  w[1] /= vRadius;
  w[2] /= vRadius;
  qBar(vRotation);
  qMultiVect(vRotation, w);
  qBar(vRotation);
  X = w[0];
  Y = w[1];
  Z = w[2];
  float h = (w[0]*vCut[0] + w[1]*vCut[1] + w[2]*vCut[2] + vCut[3])/
             (vCut[0]*vCut[0] + vCut[1]*vCut[1] + vCut[2]*vCut[2]);
  XX = w[0] - h*vCut[0];
  YY = w[1] - h*vCut[1];
  ZZ = w[2] - h*vCut[2];
}

float Director::height(float x, float y) {
  float norm = sqrt(x*x + y*y);
  if(norm > .9999f) {
    x = .9999f*x/norm;
    y = .9999f*y/norm;
  }
  //    float w[] = {x, y, -sqrt(1. - x*x - y*y)};

  float w[] = {x, y, 0.};

  w[0] -= vCenter[0];
  w[1] -= vCenter[1];
  w[2] -= vCenter[2];
  w[0] /= vRadius;
  w[1] /= vRadius;
  w[2] /= vRadius;
  qBar(vRotation);
  qMultiVect(vRotation, w);
  qBar(vRotation);

  return w[0]*vCut[0] + w[1]*vCut[1] + w[2]*vCut[2];
}
