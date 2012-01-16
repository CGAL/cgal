/****************************************************************************
 * Core Library Version 1.7, August 2004                                     
 * Copyright (c) 1995-2004 Exact Computation Project                         
 * All rights reserved.                                                      
 *                                                                           
 * This file is part of CORE (http://cs.nyu.edu/exact/core/).                
 * You can redistribute it and/or modify it under the terms of the GNU       
 * General Public License as published by the Free Software Foundation,      
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
 * File: plane3d.h
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

#ifndef _PLANE3D_H_
#define _PLANE3D_H_

#include <CGAL/CORE/geom3d/point3d.h>
#include <CGAL/CORE/geom3d/line3d.h>

class Segment3d;

class Plane3d : public GeomObj{

private:

  // ax + by + cz + d = 0
  double a; 
  double b; 
  double c; 
  double d; 
  Vector n;

public:
  /************************************************************
   *   constructors
   ************************************************************/

  Plane3d(): a(0.0), b(0.0), c(0.0), d(0.0), n(0.0, 0.0, 0.0) {}
  //trivial plane

  Plane3d(const Plane3d & plane);
  //copy constructor

  Plane3d(const Point3d & p, const Vector &v);
  // plane with direction v passes through point p 

  Plane3d(const Point3d &p1, const Point3d &p2, const Point3d &p3);
  //plane passes through points p1, p2, p3

  Plane3d(const Point3d &p, const Line3d &l);
  //plane passes through point p and line l

  Plane3d(const Point3d &p, const Segment3d &s);
  //plane passes through point p and segment s

  Plane3d(const Vector &v1, double d1);
  // plane determined by vector and displacement

  // plane determined by equation
  Plane3d(double a1, double b1, double c1, double d1);  
  
  virtual ~Plane3d() {}
/************************************************************
  *   member functions
 ************************************************************/

  virtual int dim() const { return 2; }
  
  double* coeffients() const;
  double A() const { return a; }
  double B() const { return b; }
  double C() const { return c; }
  double displacement() const { return d; }
  const Vector& normal() const { return n; }

   // test if plane is trivial
  bool isTrivial() const { return a==double(0) && b==double(0) && c==double(0); }

   // apply equation of plane to a point
  double apply( const Point3d& p ) const { return a*p.X()+b*p.Y()+c*p.Z()+d; }
  
   // plane against plane predicates
  bool isCoincident(const Plane3d& pl) const;
  bool isParallel(const Plane3d& pl) const; 
   // test parallel
  bool isParallel(const Line3d& l) const; 

  bool contains( const Point3d& p ) const;
  bool contains( const Line3d& l ) const;
  bool contains( const Segment3d& s ) const;

   // returns the projection of p on this line
  Point3d projection(const Point3d& p) const;

   /** be careful of line(segment) projection
    *  It could be a line(segment) or point
    *  The function returns degenerated line(segment) in point case
    **/
  Line3d projection(const Line3d& l) const;
  Segment3d projection(const Segment3d& s) const;

   //distance
  double distance( const Point3d& p ) const;
  double distance( const Line3d& l ) const;
  double distance( const Segment3d& s ) const;

   /** intersect predicates
    * later implementation may return like this:
    * return dimension of intersection 
    * return -1 if not intersect 
    * return 0 if intersect on a point ... and so on.
    **/
  int intersects( const Line3d& l ) const;
  int intersects( const Point3d& p ) const;
  int intersects( const Segment3d& s ) const;
  int intersects( const Plane3d& pl ) const;

   // return intersection 
  GeomObj* intersection( const Line3d& l ) const;
  GeomObj* intersection( const Segment3d& s ) const;
  GeomObj* intersection( const Plane3d& pl ) const;

  bool operator==(const Plane3d& pl) const { return isCoincident(pl); }
  bool operator!=(const Plane3d& pl) { return !operator==(pl); }

  friend std::istream& operator>>(std::istream& in, Plane3d& pl);
  friend std::ostream& operator<<(std::ostream& out, const Plane3d& pl);

}; //class Plane3d

#endif
