//********************************************
// Arcball.cpp
//********************************************
// class CArcball
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/22/00
// Modified : 09/22/00
//********************************************

#include <math.h>
#include <qgl.h>
#include "Arcball.h"

#define PI 3.14159265358979323846
//////////////////////////////////////////////
// CONSTRUCTORS
//////////////////////////////////////////////

//////////////////////////////////////////////
// DATA
//////////////////////////////////////////////

//********************************************
// Clear
//********************************************
void
CArcball::Clear()
{
  m_Center.Set(0.f,0.f,0.f);
  m_Radius = 1.f;
  m_Dragging = 0;
  m_Mode = ARCBALL_ROTATE;
  m_Show = 0;
  m_QuatTotal.Set(1.f,0.f,0.f,0.f);
  m_QuatCurrent.Set(1.f,0.f,0.f,0.f);
}

//********************************************
// Set
//********************************************
void
CArcball::Set(const float x, const float y, const float z, 
	      const float r /* = 1 */)
{
  m_Center.Set(x,y,z);
  m_Radius = r;
}

//********************************************
// Set
//********************************************
void
CArcball::Set(CVector3d center, const float r /* = 1 */)
{
  m_Center = center;
  m_Radius = r;
}

//********************************************
// Set
//********************************************
void
CArcball::Set(CArcball &a)
{
  m_Radius      = a.m_Radius;
  m_Dragging    = a.m_Dragging;
  m_Mode        = a.m_Mode;
  m_Show        = a.m_Show;

  m_VecDown.Copy(a.m_VecDown);
  m_Center.Copy(a.m_Center);
  m_QuatCurrent.Copy(a.m_QuatCurrent);
  m_QuatTotal.Copy(a.m_QuatTotal);
}

//********************************************
// Set
//********************************************
void
CArcball::Set(CArcball *pA)
{
  m_Radius      = pA->m_Radius;
  m_Dragging    = pA->m_Dragging;
  m_Mode        = pA->m_Mode;
  m_Show        = pA->m_Show;

  m_VecDown.Copy(pA->m_VecDown);
  m_Center.Copy(pA->m_Center);
  m_QuatCurrent.Copy(pA->m_QuatCurrent);
  m_QuatTotal.Copy(pA->m_QuatTotal);
}

//********************************************
// GetMatrix
//********************************************
void
CArcball::GetMatrix(float *mat)
{
  CQuaternion q = m_QuatCurrent * m_QuatTotal;
  q.GetMatrix(mat);
  mat[3]  = m_Center[0];
  mat[7]  = m_Center[1];
  mat[11] = m_Center[2];
}

//********************************************
// GetMatrix
//********************************************
CMatrix44
CArcball::GetMatrix()
{
  CQuaternion q = m_QuatCurrent * m_QuatTotal;
  CMatrix44 mat = q.GetMatrix();
  mat[0][3] = m_Center[0];
  mat[1][3] = m_Center[1];
  mat[2][3] = m_Center[2];
  return mat;
}

//********************************************
// BeginDrag
//********************************************
void
CArcball::BeginDrag(CVector3d& currentVec)
{
  if(m_Mode == ARCBALL_NO_MOTION) return;
  m_VecDown  = currentVec;
  m_Dragging = 1;
}

//********************************************
// EndDrag
//********************************************
void
CArcball::EndDrag(CVector3d& currentVec)
{
  if(!m_Dragging) return;
  Motion(currentVec);
  if(m_Mode == ARCBALL_ROTATE)
  {
    m_QuatTotal = m_QuatCurrent * m_QuatTotal;
    m_QuatCurrent.Clear();
  }
  m_VecDown.Clear();
  m_Dragging = 0;
}

//********************************************
// Motion
//********************************************
void
CArcball::Motion(CVector3d& currentVec)
{
  if(!m_Dragging) return;
  if(m_Mode == ARCBALL_ROTATE)
  {
	  m_QuatCurrent = CQuaternion(m_VecDown,currentVec);
  }
  else if(m_Mode == ARCBALL_TRANSLATE_XY)
  {
    m_Center += currentVec - m_VecDown;
    m_VecDown = currentVec;
  }
  else if(m_Mode == ARCBALL_TRANSLATE_Z)
  {
    float zDisp    = (m_VecDown[1] - currentVec[1])/30.0f; // tone down the speed
    float zCurrent = m_Center[2];
    m_Center.z(zCurrent + zDisp);
    m_VecDown = currentVec;
  }
}

//********************************************
// glDraw
//********************************************
int
CArcball::glDraw()
{
  float mat[16];
  CQuaternion q = m_QuatCurrent * m_QuatTotal;
  q.Conjugate().GetMatrix(mat);
  glTranslatef(m_Center[0],m_Center[1],m_Center[2]);
  glMultMatrixf(mat);
  if(m_Show)
  {
    glBegin(GL_POLYGON);
      glVertex3f(-m_Radius, 0.f     ,0.f);
      glVertex3f( 0.f     ,-m_Radius,0.f);
      glVertex3f( m_Radius, 0.f     ,0.f);
      glVertex3f( 0.f     , m_Radius,0.f);
    glEnd();
  }
  return 1;
}

//********************************************
// IntersectSphere
//********************************************
int
CArcball::IntersectSphere(CVector3d& rayStart, CVector3d& rayDir, 
			  CVector3d& result)
{
  CVector3d q = rayStart - m_Center;
  float a = rayDir.Dot(rayDir);
  float b = 2.f*rayDir.Dot(q);
  float c = q.Dot(q) - m_Radius*m_Radius;

  float d = b*b - 4.f*a*c;
  if(d < 0.f)
  {
    IntersectPlane(rayStart,rayDir,result);
    result -= m_Center;
    result.Normalize();
    //TRACE("mVec   : (%g,%g,%g)\n",result.x(),result.y(),result.z());
    return 0;
  }
  float t1  = (-b+sqrt(d))/(2.f*a);
  float t2  = (-b-sqrt(d))/(2.f*a);
  if (t2 < 0.f) t2 = t1;
  result = rayStart + t2 * rayDir - m_Center;
  result.Normalize();
  //TRACE("mVec   : (%g,%g,%g)\n",result.x(),result.y(),result.z());
  return 1;
}

//********************************************
// IntersectPlane
//********************************************
int 
CArcball::IntersectPlane(CVector3d& rayStart, CVector3d& rayDir, 
			 CVector3d& result, int whichPlane /* = 2 */)
{
  float denom = rayDir[whichPlane];
  if(denom == 0.f)
  {
    result = m_VecDown;
    return 0;
  }

  CVector3d numerVec = m_Center - rayStart;
  float t = (numerVec[whichPlane])/denom;
  result = rayStart + t * rayDir;
  if(t < 0.f)
  {
    result = m_VecDown;
    return 0;
  }
  return 1;
}

//********************************************
// Intersect
//********************************************
int
CArcball::Intersect(CVector3d& rayStart, CVector3d& rayDir, CVector3d& result)
{
  if(m_Mode == ARCBALL_ROTATE)
    return IntersectSphere(rayStart,rayDir,result);
  else if(m_Mode == ARCBALL_TRANSLATE_XY)
    return IntersectPlane(rayStart,rayDir,result,2);
  else if(m_Mode == ARCBALL_TRANSLATE_Z)
    return IntersectPlane(rayStart,rayDir,result,2);
  return 0;
}

//********************************************
// Intersect
//********************************************
CVector3d
CArcball::Intersect(int x, int y, CCamera& theCam, CViewport& vp)
{
  int xRes = vp.xRes();
  int yRes = vp.yRes();
  int xc   = (xRes + vp.oX())/2;
  int yc   = (yRes + vp.oY())/2;
  float r  = xRes - xc;
  r = (yRes - yc) < r ? (yRes - yc) : r;

  CVector3d theVec, p;
  float screenWidth;
  if(m_Mode == ARCBALL_ROTATE)
  {
    theVec.Set((x - xc)/r, (y - yc)/r, 0.f);
    float len = theVec.Length();
    if(len >= 1.f)
      theVec /= len;
    else
    {
      theVec.z(sqrt(1.f - len * len));
      theVec.Normalize();
    }
  }
  else if(m_Mode == ARCBALL_TRANSLATE_XY)
  {
    theVec.Set(((float)x)/xRes,((float)y)/yRes,0.f);
    p = m_Center - theCam.GetPosition();
	screenWidth = vp.GetAspectRatio() * 2.f
		*tan((theCam.GetHeightAngle()/2.f)*(PI/180.));
    theVec *= p.Dot(theCam.GetToward())*screenWidth;
  }
  else if(m_Mode == ARCBALL_TRANSLATE_Z)
  {
    theVec.Set(0.f,((float)y)/yRes,0.f);
	theVec *= .3f*(theCam.GetFarDistance() - theCam.GetNearDistance());
  }
  return theVec;
}
// ** EOF **

