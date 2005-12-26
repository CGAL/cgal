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
#include "Viewport.h"

#define PI 3.1415926535897932

class CCamera
{
protected:
  CVector3d position;
  double     radAngle;
  CVector3d rotAxis;
  double heightAngle;
  double nearDistance;
  double farDistance;

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
  void SetPosition(const double x, const double y, const double z)
    { position.Set(x,y,z); }
  void SetOrientation(double ax, double ay, double az, double ra)
    { rotAxis.Set(ax,ay,az); radAngle = ra; }
  void SetHeightAngle(const double ha)
    { heightAngle = ha; }
  void SetNearDistance(const double nd)
    { nearDistance = nd; }
  void SetFarDistance(const double fd)
    { farDistance = fd; }

  // Data access
  CVector3d GetPosition()        const { return position; }
  CVector3d GetOrientationAxis() const { return rotAxis; }
  double GetOrientationAngle()    const { return radAngle; }
  double GetHeightAngle()         const { return heightAngle; }
  double GetNearDistance()        const { return nearDistance; }
  double GetFarDistance()         const { return farDistance; }

  // Coordinate Vectors
  CVector3d GetRight();
  CVector3d GetUp();
  CVector3d GetToward();

  // Ray Casting
  CVector3d GetRayDirection(int x, int y, CViewport& vp);

  // Viewing
  void ViewAll(double xMin, double xMax,
	       double yMin, double yMax,
	       double zMin, double zMax,
	       CViewport& vp);

  // OpenGL
  void glDraw(CViewport& vp) const;
};

#endif // _CAMERA_


