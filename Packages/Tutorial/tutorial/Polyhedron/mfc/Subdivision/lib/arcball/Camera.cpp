//********************************************
// Camera.cpp
//********************************************
// class CCamera
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/07/00
// Modified : 09/07/00
//********************************************

#include "stdafx.h"
#include "Camera.h"
#include <math.h>


//////////////////////////////////////////////
// DATA SETTING
//////////////////////////////////////////////

//********************************************
// Set
//********************************************
void
CCamera::Set(CCamera& c)
{
  position     = c.position;
  radAngle     = c.radAngle;
  rotAxis      = c.rotAxis;
  heightAngle  = c. heightAngle;
  nearDistance = c.nearDistance;
  farDistance  = c.farDistance;
}

//********************************************
// Set
//********************************************
void
CCamera::Set(CCamera *pC)
{
  position     = pC->position;
  radAngle     = pC->radAngle;
  rotAxis      = pC->rotAxis;
  heightAngle  = pC->heightAngle;
  nearDistance = pC->nearDistance;
  farDistance  = pC->farDistance;
}

//////////////////////////////////////////////
// COORDINATE VECTORS
//////////////////////////////////////////////

//********************************************
// GetRight
//********************************************
CVector3d
CCamera::GetRight()
{
  CMatrix44 rotMat;
  rotMat.SetRotate(rotAxis[0],rotAxis[1],rotAxis[2],radAngle);
  return CVector3d(rotMat[0][0],rotMat[1][0],rotMat[2][0]);
}

//********************************************
// GetUp
//********************************************
CVector3d
CCamera::GetUp()
{
  CMatrix44 rotMat;
  rotMat.SetRotate(rotAxis[0],rotAxis[1],rotAxis[2],radAngle);
  return CVector3d(rotMat[0][1],rotMat[1][1],rotMat[2][1]);
}

//********************************************
// GetToward
//********************************************
CVector3d
CCamera::GetToward()
{
  CMatrix44 rotMat;
  rotMat.SetRotate(rotAxis[0],rotAxis[1],rotAxis[2],radAngle);
  return CVector3d(-rotMat[0][2],-rotMat[1][2],-rotMat[2][2]);
}

//////////////////////////////////////////////
// RAY CASTING
//////////////////////////////////////////////

//********************************************
// GetRayDirection
//  Returns the ray direction given a point 
//  on the image plane.
//********************************************
CVector3d
CCamera::GetRayDirection(int x, int y, CViewport& vp)
{
  int xRes = vp.xRes();
  int yRes = vp.yRes();
  double aspectRatio = vp.GetAspectRatio();

  double h = 2.f * nearDistance * tan((heightAngle/2.f)*(PI/180.));
  double w = aspectRatio * h;

  double rx = (((double)x)/((double)(xRes-1))-0.5f) * w;
  double ry = (((double)y)/((double)(yRes-1))-0.5f) * h;
  double rz = -nearDistance;

  CVector3d rayDir(rx,ry,rz);
  CMatrix44 rotMat;
  rotMat.SetRotate(rotAxis[0],rotAxis[1],rotAxis[2],radAngle);
  return (rotMat*rayDir);
}

//////////////////////////////////////////////
// VIEWING
//////////////////////////////////////////////

//********************************************
// ViewAll
//  Sets the camera to view the bounding box.
//********************************************
void
CCamera::ViewAll(double xMin, 
								 double xMax,
								 double yMin, 
								 double yMax,
								 double zMin, 
								 double zMax,
								 CViewport& vp)
{
  // Get bounding sphere of bounding box
  CVector3d boxMin(xMin, yMin, zMin);
  CVector3d boxMax(xMax, yMax, zMax);
  CVector3d center = 0.5f * (boxMin + boxMax);
  double radius = (boxMax - center).Length();

  // Get the distance required to fit the object in view
  double hOverD = tan((heightAngle/2.f)*(PI/180.));
  double aspectRatio = vp.GetAspectRatio();
  if(aspectRatio < 1.f)
    hOverD *= aspectRatio;
  double distToCenter = radius / hOverD;

  // Set the internal variables appropriately
  CVector3d toward = GetToward();
  position = center - distToCenter * toward;
  //nearDistance = distToCenter - radius;
  //farDistance  = distToCenter + radius;
  nearDistance = 0.1;
  farDistance  = 100.0;
}

//////////////////////////////////////////////
// OPENGL
//////////////////////////////////////////////

//********************************************
// glDraw
//********************************************
void
CCamera::glDraw(CViewport& vp) const
{
  double aspectRatio = vp.GetAspectRatio();

  // Position and Orientation
  glMatrixMode(GL_MODELVIEW);
  glRotated(radAngle*(360.f/PI),rotAxis[0],rotAxis[1],rotAxis[2]);
  glTranslated(-position.x(),-position.y(),-position.z());

  // Projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(heightAngle,aspectRatio,nearDistance,farDistance);
  glMatrixMode(GL_MODELVIEW);
}

// ** EOF **
