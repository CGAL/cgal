//********************************************
// Arcball.cpp
//********************************************
// class CArcball
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/22/00
// 02/20/03 pierre.alliez@sophia.inria.fr
//********************************************

#include "stdafx.h"
#include <math.h>
#include "Arcball.h"

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
  m_Show = 0; // hide by default
  m_QuatTotal.Set(1.f,0.f,0.f,0.f);
  m_QuatCurrent.Set(1.f,0.f,0.f,0.f);
}

//********************************************
// Set
//********************************************
void
CArcball::Set(const double x, const double y, const double z, 
	      const double r /* = 1 */)
{
  m_Center.Set(x,y,z);
  m_Radius = r;
}

//********************************************
// Set
//********************************************
void
CArcball::Set(CVector3d center, const double r /* = 1 */)
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
  m_Center   = a.m_Center;
  m_Radius   = a.m_Radius;
  m_Dragging = a.m_Dragging;
  m_Mode     = a.m_Mode;
  m_Show     = a.m_Show;

  m_VecDown.Set(a.m_VecDown);
  m_Center.Set(a.m_Center);
  m_QuatCurrent.Set(a.m_QuatCurrent);
  m_QuatTotal.Set(a.m_QuatTotal);
}

//********************************************
// Set
//********************************************
void
CArcball::Set(CArcball *pA)
{
  m_Center   = pA->m_Center;
  m_Radius   = pA->m_Radius;
  m_Dragging = pA->m_Dragging;
  m_Mode     = pA->m_Mode;
  m_Show     = pA->m_Show;

  m_VecDown.Set(pA->m_VecDown);
  m_Center.Set(pA->m_Center);
  m_QuatCurrent.Set(pA->m_QuatCurrent);
  m_QuatTotal.Set(pA->m_QuatTotal);
}

//********************************************
// GetMatrix
//********************************************
void
CArcball::GetMatrix(double *mat)
{
  CQuaternion q = m_QuatCurrent * m_QuatTotal;
  q.GetMatrix(mat);
  mat[3]  = m_Center[0];
  mat[7]  = m_Center[1];
  mat[11] = m_Center[2];
}

void
CArcball::GetMatrix(double mat[4][4])
{
  CQuaternion q = m_QuatCurrent * m_QuatTotal;
	double m[16];
  q.GetMatrix(m);
  m[3]  = m_Center[0];
  m[7]  = m_Center[1];
  m[11] = m_Center[2];

	mat[0][0] = m[0];
	mat[0][1] = m[1];
	mat[0][2] = m[2];
	mat[0][3] = m[3];

	mat[1][0] = m[4];
	mat[1][1] = m[5];
	mat[1][2] = m[6];
	mat[1][3] = m[7];

	mat[2][0] = m[8];
	mat[2][1] = m[9];
	mat[2][2] = m[10];
	mat[2][3] = m[11];

	mat[3][0] = m[12];
	mat[3][1] = m[13];
	mat[3][2] = m[14];
	mat[3][3] = m[15];
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
// Trace
//********************************************
void
CArcball::Trace() const
{
  TRACE("\n");
  TRACE("** Arcball **\n");
  TRACE("Address : %x\n",this);
  TRACE("Center  : (%g %g %g)\n",m_Center.x(),m_Center.y(),m_Center.z());
  TRACE("Radius  : %g\n",m_Radius);
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
    double zDisp    = m_VecDown[1] - currentVec[1];
    double zCurrent = m_Center[2];
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
  double mat[16];
  CQuaternion q = m_QuatCurrent * m_QuatTotal;
  q.Conjugate().GetMatrix(mat);
  glTranslated(m_Center[0],m_Center[1],m_Center[2]);
  glMultMatrixd(mat);
  if(m_Show)
  {
    glDisable(GL_LIGHTING);
    glBegin(GL_POLYGON);
      glVertex3d(-m_Radius, 0     ,0);
      glVertex3d( 0     ,-m_Radius,0);
      glVertex3d( m_Radius, 0     ,0);
      glVertex3d( 0     , m_Radius,0);
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
  double a = rayDir.Dot(rayDir);
  double b = 2.f*rayDir.Dot(q);
  double c = q.Dot(q) - m_Radius*m_Radius;

  double d = b*b - 4.f*a*c;
  if(d < 0.f)
  {
    IntersectPlane(rayStart,rayDir,result);
    result -= m_Center;
    result.Normalize();
    TRACE("mVec   : (%g,%g,%g)\n",result.x(),result.y(),result.z());
    return 0;
  }
  double t1  = (-b+sqrt(d))/(2.f*a);
  double t2  = (-b-sqrt(d))/(2.f*a);
  if (t2 < 0.f) t2 = t1;
  result = rayStart + t2 * rayDir - m_Center;
  result.Normalize();
  TRACE("mVec   : (%g,%g,%g)\n",result.x(),result.y(),result.z());
  return 1;
}

//********************************************
// IntersectPlane
//********************************************
int 
CArcball::IntersectPlane(CVector3d& rayStart, CVector3d& rayDir, 
			 CVector3d& result, int whichPlane /* = 2 */)
{
  double denom = rayDir[whichPlane];
  if(denom == 0.f)
  {
    result = m_VecDown;
    return 0;
  }

  CVector3d numerVec = m_Center - rayStart;
  double t = (numerVec[whichPlane])/denom;
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
  double r  = xRes - xc;
  r = (yRes - yc) < r ? (yRes - yc) : r;

  CVector3d theVec, p;
  double screenWidth;
  if(m_Mode == ARCBALL_ROTATE)
  {
    theVec.Set((x - xc)/r, (y - yc)/r, 0.f);
    double len = theVec.Length();
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
    theVec.Set(((double)x)/xRes,((double)y)/yRes,0.f);
    p = m_Center - theCam.GetPosition();
	screenWidth = vp.GetAspectRatio() * 2.f
		*tan((theCam.GetHeightAngle()/2.f)*(PI/180.));
    theVec *= p.Dot(theCam.GetToward())*screenWidth;
  }
  else if(m_Mode == ARCBALL_TRANSLATE_Z)
  {
    theVec.Set(0.f,((double)y)/yRes,0.f);
	theVec *= .3f*(theCam.GetFarDistance() - theCam.GetNearDistance());
  }
  return theVec;
}
// ** EOF **

