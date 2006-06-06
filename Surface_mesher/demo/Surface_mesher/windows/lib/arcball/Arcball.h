///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CArcball                                                      //
//                                                                       //
//  An arcball class used to translate and rotate an object or scene.    //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _ARCBALL_
#define _ARCBALL_

#include "Camera.h"
#include "Quaternion.h"
#include "Vector3d.h"
#include "Viewport.h"

#define ARCBALL_NO_MOTION    0
#define ARCBALL_ROTATE       1
#define ARCBALL_TRANSLATE_XY 2
#define ARCBALL_TRANSLATE_Z  3
#define PI 3.1415926535897932

class CArcball
{
private :
  // Data
  CQuaternion m_QuatTotal, m_QuatCurrent;
  CVector3d   m_VecDown;
  CVector3d   m_Center;
  double      m_Radius;
  int         m_Dragging;
  int         m_Mode;
  int         m_Show;

public :

  // Constructors
  CArcball() { Clear(); }
  CArcball(const double x, const double y, const double z, const double r = 1.0f)
    { Clear(); Set(x,y,z,r); }
  CArcball(CVector3d center, const double r = 1.0f)
    { Clear(); Set(center,r); }
  CArcball(CArcball& a)  { Set(a); }
  CArcball(CArcball *pA) { Set(pA); }

  virtual ~CArcball() { }

  // Data setting
  void Clear();
  void Set(const double x, const double y, const double z, 
	   const double r = 1.0f);
  void Set(CVector3d center, const double r = 1.0f);
  void Set(CArcball &arcball);
  void Set(CArcball *pArcball);

  void SetMode(int mode) { m_Mode = mode; }
  int  GetMode()         { return m_Mode; }
  void Show()		 { m_Show = 1; }
  void Hide()		 { m_Show = 0; }

  void      SetCenter(const double x, const double y, const double z)
    { m_Center.Set(x,y,z); }
  void      SetCenter(CVector3d center)
    { m_Center = center; }
  void      SetRadius(const double r) { m_Radius = r; }

  CVector3d GetCenter() const { return m_Center; }
  double GetRadius() const { return m_Radius; }

  void GetMatrix(double *mat);
	void GetMatrix(double mat[4][4]);
  CMatrix44 GetMatrix();

  // Debug
  void Trace() const;

  void BeginDrag(CVector3d& currentVec);
  void EndDrag(CVector3d& currentVec);
  void Motion(CVector3d& currentVec);
  int  glDraw();

  int IntersectSphere(CVector3d& rayStart, CVector3d& rayDir, 
		      CVector3d& result);
  int IntersectPlane(CVector3d& rayStart, CVector3d& rayDir, 
		     CVector3d& result, int whichPlane = 2);
  int Intersect(CVector3d& rayStart, CVector3d& rayDir, CVector3d& result);
  CVector3d Intersect(int x, int y, CCamera& theCam, CViewport& vp);
};


#endif // 

