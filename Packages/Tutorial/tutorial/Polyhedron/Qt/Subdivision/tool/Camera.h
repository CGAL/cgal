///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CCamera                                                       //
//                                                                       //
//  Camera class to represent virtual cameras.                           //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _CAMERA_
#define _CAMERA_

#include "Quaternion.h"
#include "Vector3d.h"
class CViewport;

class CCamera
{

private :

protected:
  CVector3d position;
  float     radAngle;
  CVector3d rotAxis;
  float heightAngle;
  float nearDistance;
  float farDistance;

public :

  // Constructor
  CCamera()                   { }
  CCamera(const CVector3d& p) { SetPosition(p); }
  CCamera(CCamera& c)         { Set(c); }
  CCamera(CCamera *pC)        { Set(pC); }

  ~CCamera() { }

  // Data setting
  void Set(CCamera& c);
  void Set(CCamera *pC);
  void SetPosition(const CVector3d& p)
    { position = p; }
  void SetPosition(const float x, const float y, const float z)
    { position.Set(x,y,z); }
  void SetOrientation(float ax, float ay, float az, float ra)
    { rotAxis.Set(ax,ay,az); radAngle = ra; }
  void SetHeightAngle(const float ha)
    { heightAngle = ha; }
  void SetNearDistance(const float nd)
    { nearDistance = nd; }
  void SetFarDistance(const float fd)
    { farDistance = fd; }

  // Data access
  CVector3d GetPosition()        const { return position; }
  CVector3d GetOrientationAxis() const { return rotAxis; }
  float GetOrientationAngle()    const { return radAngle; }
  float GetHeightAngle()         const { return heightAngle; }
  float GetNearDistance()        const { return nearDistance; }
  float GetFarDistance()         const { return farDistance; }

  // Coordinate Vectors
  CVector3d GetRight();
  CVector3d GetUp();
  CVector3d GetToward();

  // Ray Casting
  CVector3d GetRayDirection(int x, int y, CViewport *vp);

  // Viewing
  void ViewAll(float xMin, float xMax,
	       float yMin, float yMax,
	       float zMin, float zMax,
	       CViewport *vp);

  // OpenGL
  void glDraw(CViewport *vp) const;
};

#endif // _CAMERA_


