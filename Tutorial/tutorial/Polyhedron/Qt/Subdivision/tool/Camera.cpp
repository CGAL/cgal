//********************************************
// Camera.cpp
//********************************************
// class CCamera
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/07/00
// Modified : 09/07/00
//********************************************

#include "Camera.h"
#include "Viewport.h"
#include <math.h>
#include <qgl.h>

#define PI 3.14159265358979323846

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
CCamera::GetRayDirection(int x,
                         int y,
                         CViewport *vp)
{
  int xRes = vp->xRes();
  int yRes = vp->yRes();
  float aspectRatio = vp-> GetAspectRatio();

  float h = 2.f * nearDistance * tan((heightAngle/2.f)*(PI/180.));
  float w = aspectRatio * h;

  float rx = (((float)x)/((float)(xRes-1))-0.5f) * w;
  float ry = (((float)y)/((float)(yRes-1))-0.5f) * h;
  float rz = -nearDistance;

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
CCamera::ViewAll(float xMin, float xMax,
		 float yMin, float yMax,
		 float zMin, float zMax,
		 CViewport *vp)
{
  // Get bounding sphere of bounding box
  CVector3d boxMin(xMin, yMin, zMin);
  CVector3d boxMax(xMax, yMax, zMax);
  CVector3d center = 0.5f * (boxMin + boxMax);
  float radius = (boxMax - center).Length();

  // Get the distance required to fit the object in view
  float hOverD = tan((heightAngle/2.f)*(PI/180.));
  float aspectRatio = vp->GetAspectRatio();
  if(aspectRatio < 1.f)
    hOverD *= aspectRatio;
  float distToCenter = radius / hOverD;

  // Set the internal variables appropriately
  CVector3d toward = GetToward();
  position = center - distToCenter * toward;
  nearDistance = distToCenter - radius;
  farDistance  = distToCenter + radius;
}

//////////////////////////////////////////////
// OPENGL
//////////////////////////////////////////////

//********************************************
// glDraw
//********************************************
void
CCamera::glDraw(CViewport *vp) const
{
  float aspectRatio = vp->GetAspectRatio();

  // Position and Orientation
  glMatrixMode(GL_MODELVIEW);
  glRotatef(radAngle*(360.f/PI),rotAxis[0],rotAxis[1],rotAxis[2]);
  glTranslatef(-position.x(),-position.y(),-position.z());

  // Projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(heightAngle,aspectRatio,nearDistance,farDistance);
  glMatrixMode(GL_MODELVIEW);
}

// ** EOF **
