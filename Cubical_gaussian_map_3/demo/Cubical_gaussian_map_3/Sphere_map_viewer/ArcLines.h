// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 *
 * part of the viewer of spherical arrangements
 */


#ifndef _ARC_LINES_H_
#define _ARC_LINES_H_

#include <math.h>

/*
 normalize a vector to have a unit length at the same direction

 v - the vector to be normalized
*/
inline
void normalize(float v[3]) {
  float d = sqrtf(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (d == 0.0) {
//      error("zero length vector");
    return;
   }
   v[0] /= d; v[1] /= d; v[2] /= d;
}

/*
 this class handles the segmentation of a spherical arc to line segments
*/
class ArcLines
{
public:

  /*
   constructor
  */
  ArcLines(void)
  {
  }

  /*
   constructor,
   stDir, enDir - the arc endpoints
  */
  ArcLines(float stDir[], float enDir[]) {
    putArc(stDir, enDir);
  }

  ~ArcLines(void)
  {
  }

  /*
   update internal arc to be between points stPnt to enPnt,
   assuming spPnt and enPnt have unit length

   stPnt, enPnt - the arc endpoints
  */
  void putArc(float stPnt[], float enPnt[]) {
    m_u[0]=stPnt[0]; m_u[1]=stPnt[1]; m_u[2]=stPnt[2];
    float arcNorm[3];
    normcrossprod(stPnt, enPnt, arcNorm);
    float normLen = sqrtf(arcNorm[0]*arcNorm[0]+arcNorm[1]*arcNorm[1]+
      arcNorm[2]*arcNorm[2]);
    arcNorm[0]/=normLen; arcNorm[1]/=normLen; arcNorm[2]/=normLen;
    float cosLen = enPnt[0]*stPnt[0]+enPnt[1]*stPnt[1]+
      enPnt[2]*stPnt[2];
    m_angle = atan2f(normLen, cosLen);
    normcrossprod(arcNorm, m_u, m_v);
  }

  /*
   get a point on the arc with ang angle compare to start direction
   at ang = 0 the start endpoint is received,
   at ang = angle between endpoints directions, the ending endpoint is reeived

   ang - angle between point to start point on the arc surface
   pnt - will hold the found point
  */
  inline void getPoint(float ang, float pnt[3]) {
    float cosAng = cosf(ang),
      sinAng = sinf(ang);
    pnt[0]=cosAng*m_u[0]+sinAng*m_v[0];
    pnt[1]=cosAng*m_u[1]+sinAng*m_v[1];
    pnt[2]=cosAng*m_u[2]+sinAng*m_v[2];
  }

  /*
   get the angle between arc endpoints directions

   return value - the angle between arc endpoints directions
  */
  inline float getAngle() {
    return m_angle;
  }
private:
  /*
   calculate the cross product of v1 with v2 to out,
   out := v1xv2

   v1,v2 - vectors to cross product
   out - will hold the cross product value, v1xv2
  */
  inline
  void normcrossprod(float v1[3], float v2[3], float out[3]) {
    out[0] = v1[1]*v2[2] - v1[2]*v2[1];
    out[1] = v1[2]*v2[0] - v1[0]*v2[2];
    out[2] = v1[0]*v2[1] - v1[1]*v2[0];
  }

  // the u and v vectors of the surface containing the two directions
  float m_u[3], m_v[3];
  float m_angle; // the angle between the end points of the arc
};

#endif
