// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 */

#ifndef CGAL_SPHERE_ARC_H
#define CGAL_SPHERE_ARC_H

#include <CGAL/basic.h>

#include <ostream>

CGAL_BEGIN_NAMESPACE

template <class Kernel_> class Sphere_traits;

#include <CGAL/Direction_3.h>
#include <CGAL/Cartesian/Direction_3.h>

#ifdef WIN32
/*
 to make this project compile on windows platform it was required for the
 output and input operators of Direction_3 to be redefined
*/
#include <CGAL/Cartesian.h>

template <class _NumberType>
std::ostream &operator<<(std::ostream &os, 
             const CGAL::Direction_3<CGAL::Cartesian<_NumberType> > &dir) {
  return os << dir.dx() << "," << dir.dy() << "," << dir.dz();
}

template <class _NumberType>
std::istream &operator>>(std::istream &is, 
             CGAL::Direction_3<CGAL::Cartesian<_NumberType> > &dir) {
  _NumberType dx, dy, dz;
  is >> dx >> dy >> dz;
  dir = CGAL::Direction_3<CGAL::Cartesian<_NumberType> >(dx, dy, dz);
  return is;
}

#endif

/*
  represents an arc of a great circle on the unit sphere,
  an arc is represented by the directions of it's end points and the direction
  normal to the plane containing the arc, the arc is the part on the circle
  plane when moving from the start point to the end point counter clock-wise
  in relation to the direction normal to the plane
 */
template <class Kernel_>
class Sphere_arc {
public:
  // public definitions
  typedef typename Kernel_::Direction_3   Direction_3;
  typedef typename Kernel_::Vector_3     Vector_3;
  typedef typename Kernel_::Plane_3     Plane_3;
  typedef typename Kernel_::Point_3     Point_3;
  typedef typename Kernel_::Line_3     Line_3;

  /*
     constructor,
   gets the two arc's endpoints directions and a boolean stating if the arc
   is the shorter or longer arc of great circle between these directions

   st, en - arc's endpoints directions
   shortArc - true if the arc is the short arc between these two endpoints
  */
  Sphere_arc(Direction_3 st, Direction_3 en, bool shortArc=true);

  /*
     constructor,
   gets the two arc's endpoints directions and a direction for the 
   normal to the arc's plane, the arc will be between the endpoints
   and it's orientation is starting from the start endpoint to the second
   endpoint going counter clock wise with regards to the plane's normal direction

   st, en - arc's endpoints direction
   planeDir - the direction of the normal to the arc's plane
  */ 
  Sphere_arc(Direction_3 st, Direction_3 en, Direction_3 planeDir);
  
  /*
   empty constructor
  */
  Sphere_arc() {}


  /*
   print the arc to an os output stream with format 
   {normal_to_plane_direction}: (start_direction) -> (end_direction)

   os - output stream to write arc values to
   return value - the output stream with the written arc's information
  */  
  std::ostream & print(std::ostream & os) const;
  
  
  /*
   get the direction of the start point  

   return value - the direction of the arc's starting endpoint
  */
  inline Direction_3 getStDir() const;

  
  /*
   get the direction of the end point

   return value - the direction of the arc's target endpoint
  */
  inline Direction_3 getEnDir() const;

  /*
   get the direction of the normal to the arc's circle plane

   return value - the direction of the normal to the arc's plane
  */
  inline Direction_3 getPlaneDir() const;

private:  
  Direction_3 m_dirSt; // direction of the curve start
  Direction_3 m_dirEn; // direction of the curve end
  Direction_3 m_planeDir; // direction of the normal to the plane of the arc
};

/*
 get the direction of the start point  

 return value - the direction of the arc's starting endpoint
*/
template <class Kernel_>
inline
typename Sphere_arc<Kernel_>::Direction_3 Sphere_arc<Kernel_>::
getStDir() const
{
  return m_dirSt;
}

/*
 get the direction of the end point

 return value - the direction of the arc's target endpoint
*/
template <class Kernel_>
inline
typename Sphere_arc<Kernel_>::Direction_3 Sphere_arc<Kernel_>::
getEnDir() const
{
  return m_dirEn;
}

/*
 get the direction of the normal to the arc's circle plane

 return value - the direction of the normal to the arc's plane
*/
template <class Kernel_>
inline
typename Sphere_arc<Kernel_>::Direction_3 Sphere_arc<Kernel_>::
getPlaneDir() const
{
  return m_planeDir;
}


/*
    constructor,
  gets the two arc's endpoints directions and a boolean stating if the arc
  is the shorter or longer arc of great circle between these directions (if
  endpoints are in opposite directions, use c'tor with plane normal)

  st, en - arc's endpoints direction
  shortArc - true if the arc is the short arc between these two endpoints
*/
template <class Kernel_>
Sphere_arc<Kernel_>::Sphere_arc
(Direction_3 st, Direction_3 en, bool shortArc) :
  m_dirSt(st), m_dirEn(en) {
  
  if (st == en || st == -en) {
    // opposite endpoint directions, can't determine the plane of the arc
    std::cerr << "ERROR: arc's plane cannot be decided" << std::endl;
    CGAL_assertion(false);    
  }

  Vector_3 stVec = m_dirSt.vector();
  Vector_3 enVec = m_dirEn.vector();
  // a direction perpendicular to arc plane, if short arc plane in this direction else
  // plane in the opposite direction
  Vector_3 plNorm = CGAL::cross_product(stVec, enVec);
  if (shortArc) {
    m_planeDir = Direction_3(plNorm);
  } else {
    m_planeDir = Direction_3(-plNorm);
  }
};

/*
 constructor,
 gets the two arc's endpoints directions and a direction for the 
 normal to the arc's plane, the arc will be between the endpoints
 and it's orientation is starting from the start endpoint to the second
 endpoint going counter clock wise with regards to the plane's normal direction

 st, en - arc's endpoints direction
 planeDir - the direction of the normal to the arc's plane
*/ 
template <class Kernel_>
Sphere_arc<Kernel_>::Sphere_arc
(Direction_3 st, Direction_3 en, Direction_3 planeDir) :
  m_dirSt(st), m_dirEn(en), m_planeDir(planeDir) {
    
  // check feasability of arc parameters
  if (st != -en) {
    // not opposite directions, check if normal plane is arc's plane
    Direction_3 dir = CGAL::cross_product(st.vector(), en.vector());
    if (dir == planeDir || dir == -planeDir) {
    } else {
      std::cerr << "ERROR: the arc points are not on the plane" << std::endl;
      CGAL_assertion(false);
    }
  } else { // the two points are in opposite directions
    Plane_3 plane(Point_3(0,0,0), planeDir);
    // check if the points are on the plane
    if (!plane.has_on(Line_3(Point_3(0,0,0), st)) ||
      !plane.has_on(Line_3(Point_3(0,0,0), en))) {
      
      std::cerr << "ERROR: the arc points are not on the plane" << std::endl;
      CGAL_assertion(false);
    }
  }
};

/*
  print the arc to os stream with format 
  {normal_to_plane_direction}: (start_direction) -> (end_direction)

  os - output stream to write arc to
  return value - the output stream with the appended arc information
*/
template <class Kernel_>
std::ostream & Sphere_arc<Kernel_>::print(std::ostream & os) const {  
  return os << "{" << m_planeDir << "} :(" << m_dirSt << ") -> (" << 
    m_dirEn << ")";
}

/*
 overloading of the output operator to display the arc
 displaying an arc will show the direction normal to arc's circle plane
 and the directions of the start and end points of the arc

 os - output stream to write the arc into
 arc - the arc to be written
 return value - the output stream with the appended arc information
*/
template <class Kernel_>
std::ostream & operator<<(std::ostream &os,
        const Sphere_arc<Kernel_> &arc
) {
  return arc.print(os);
}

CGAL_END_NAMESPACE

#endif
