/****************************************************************************
 * Core Library Version 1.7, August 2004                                     
 * Copyright (c) 1995-2004 Exact Computation Project                         
 * All rights reserved.                                                      
 *                                                                           
 * This file is part of CGAL (www.cgal.org).                
 * You can redistribute it and/or modify it under the terms of the GNU       
 * Lesser General Public License as published by the Free Software Foundation,      
 * either version 3 of the License, or (at your option) any later version.   
 *                                                                           
 * Licensees holding a valid commercial license may use this file in         
 * accordance with the commercial license agreement provided with the        
 * software.                                                                 
 *                                                                           
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE   
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. 
 *                                                                           
 *                                                                           
 * $URL$                                                                     
 * $Id$                                                                      
 ***************************************************************************/
/*****************************************************************
 * File: triangle3d.h
 * Synopsis:
 *      Basic 3-dimensional geometry
 * Author: Shubin Zhao (shubinz@cs.nyu.edu), 2001.
 *
 *****************************************************************
 * CORE Library Version 1.4 (July 2001)
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * Copyright (c) 1995, 1996, 1998, 1999, 2000, 2001 Exact Computation Project
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/


#ifndef _TRIANGLE3D_H_
#define _TRIANGLE3D_H_

#include <CGAL/CORE/geom3d/point3d.h>
#include <CGAL/CORE/geom3d/line3d.h>
#include <CGAL/CORE/geom3d/segment3d.h>
#include <CGAL/CORE/geom3d/plane3d.h>
#include <CGAL/CORE/geom3d/polygon3d.h>

class Triangle3d : public GeomObj{

private:
  
  // three vertices
  Point3d p0; 
  Point3d p1; 
  Point3d p2;

public:

  /************************************************************
   *   constructors
   ************************************************************/

  Triangle3d(const Point3d& v1, const Point3d& v2, const Point3d& v3);
  // given three vertices

  Triangle3d(const Triangle3d& T);
  // given a triangle

  Triangle3d(): p0(ORIGIN_3D), p1(ORIGIN_3D), p2(ORIGIN_3D) {}
  //trivial triangle

  virtual ~Triangle3d() {}
  /************************************************************
   *   member functions
   ************************************************************/

   // view a triangle as a surface
  virtual int dim() const { return 2; }
   
  Point3d V1() const { return p0; }
  Point3d V2() const { return p1; }
  Point3d V3() const { return p2; }
  
  Vector normal() const { return (p1 - p0).cross( p2 - p0); }
  // return normal of the plane containing this triangle

  Plane3d toPlane() const { return Plane3d(p0,p1,p2); }  
 
  Polygon3d* toPolygon() const;
  /************************************************************
   *   predicates
   ************************************************************/
   
  inline bool isCoplanar( const Point3d& p ) const   { 
                return orientation3d(p0, p1, p2, p) == 0; }
  inline bool isCoplanar( const Segment3d& s ) const { 
                return isCoplanar(s.startPt()) && isCoplanar(s.stopPt()); }
  inline bool isCoplanar( const Line3d& l )   const  { 
                return isCoplanar(l.startPt()) && isCoplanar(l.stopPt()); }
  inline bool isCoplanar( const Triangle3d& T ) const {
                return isCoplanar(T.V1()) && isCoplanar(T.V2()) && isCoplanar(T.V3()); }
  inline bool isCoplanar( const Plane3d& pl )  const { 
                return pl.contains(p0) && pl.contains(p1) && pl.contains(p2); }

   // test if p is on triangle 
  bool contains( const Point3d& p ) const;
   // test if s is on triangle 
  inline bool contains( const Segment3d& s ) const  {
               return contains( s.startPt() ) && contains( s.stopPt() ); }

  // test if T is on triangle 
  inline bool contains( const Triangle3d& T ) const  {
               return contains( T.V1() ) && contains( T.V2() ) && contains( T.V3() ); }


  bool isOnEdge( const Point3d& p ) const;
  // test if p is on the edge

  bool inside( const Point3d& p ) const;
  // test if p is inside the triangle
  
  /************************************************************
   *  Intersection
   ************************************************************/

  /** all intersect predicates return the dimension of the intersection 
    * -1 if disjoint (i.e., parallel but distinct lines)
    * 0  if coincident
    * 0  if intersect in a point.  In this case, the
    *
    // 	-1 if disjoint (i.e., parallel but distinct lines)
    //  1  if coincident
    //  0  if intersect in a point.  In this case, the
    //		intersection point is assigned to p if this is available.
  **/
  
   // intersect predicates
  bool do_intersect( const Segment3d& s ) const;
  bool do_intersect( const Line3d& l ) const;
  bool do_intersect( const Plane3d& pl ) const;
  bool do_intersect( const Triangle3d& t ) const;
  
   // these are consistent with other classes
   // they return the dimension of the intersection
   // return -1 if no intersection
  int  intersects( const Segment3d& s ) const;
  int  intersects( const Line3d& l ) const;
  int  intersects( const Plane3d& pl ) const;
  int  intersects( const Triangle3d& T ) const;

   // general intersections
  GeomObj* intersection( const Segment3d& s ) const;
  GeomObj* intersection( const Line3d& l ) const;
  GeomObj* intersection( const Plane3d& pl ) const;
  GeomObj* intersection( const Triangle3d& t ) const;

   // coplanar intersections
  GeomObj* coplanar_intersection( const Segment3d& s ) const;
  GeomObj* coplanar_intersection( const Line3d& l ) const;
  GeomObj* coplanar_intersection( const Triangle3d& T ) const;

  Polygon3d* in_half_plane( const Point3d& pa, 
                           const Point3d& pb, 
                           const Point3d& pSide,
                           Polygon3d& plg ) const;

  int coplanar_orientation( const Point3d& pa, const Point3d& pb, 
                            const Point3d& ps, const Point3d& p ) const;
  /************************************************************
   *   I/O 
   ************************************************************/

  friend std::istream& operator>>(std::istream& in, Triangle3d& T);
  friend std::ostream& operator<<(std::ostream & out, const Triangle3d & T);
  
}; //class Triangle3d

#endif
