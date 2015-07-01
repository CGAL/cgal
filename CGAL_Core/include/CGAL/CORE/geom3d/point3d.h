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
 * File: point3d.h
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

#ifndef _POINT3D_H
#define _POINT3D_H

#include <CGAL/CORE/CORE.h>
#include <CGAL/CORE/linearAlgebra.h>
#include <CGAL/CORE/geombase.h>

 // class defination for 3d points
class Point3d : public GeomObj{

private:
  double x, y, z;
 
public:
 
  /************************************************************
   * constructors and destructors
   ************************************************************/
  Point3d();  //initialized to origin(0,0,0)
  Point3d(double x, double y, double z);
  Point3d(const Point3d & p);
  Point3d(const Vector& v);
    //create a point initialized to the point $(v[0], v[1], v[2])$ 
    //precondition: v.dim() >= 3 (only the first 3 components are used)

  //destructor
  virtual ~Point3d() {}

  /************************************************************
   * Methods
   ************************************************************/
  Point3d& operator=(const Point3d&);
  
  double X() const { return x; }
  double Y() const { return y; } 
  double Z() const { return z; } 

  Vector toVector() const { return Vector(x, y, z); } 

  virtual int dim() const { return 0; }

  double distance(const Point3d& p) const;
    // returns the Euclidean distance between p and this

  double distance() const { return distance(Point3d(0, 0, 0)); }
    // returns distance between this and origin

  Point3d negate() const { return Point3d(-x, -y, -z); }
  Vector operator-(const Point3d &p) const;
  Point3d operator+(const Vector &v) const;
  Point3d operator-(const Vector &v) const;
  Point3d operator*(const double& d) const;
 
 /************************************************************
   * predicates
   ************************************************************/

  bool operator==(const Point3d&) const;
  bool operator!=(const Point3d& p) const {return !operator==(p); }
  
  /************************************************************
   * I/O, debugging
   ************************************************************/
  friend std::ostream& operator<< (std::ostream&, const Point3d&);
    // write point p to output stream

  friend std::istream& operator>>(std::istream&, Point3d&);
    // reads the x and y coordinates of point p from the input stream

  // routines to display point:
  void dump() const {
    std::cout << "(" << x <<", " << y << z << ")" ; // simply outputs "(x, y)"
  }

  void dump(const char* s) const {
    std::cout << s << "(" << x <<", " << y << z << ")" ; // s is the prefix message
  }

  void dump(const char* s, const char* ss) const {
    std::cout << s << "(" << x <<", " << y << z << ss ; // ss is the suffix message
  }

   // compute signed volume of a tetrahedron
  friend double signed_volume(const Point3d& a, const Point3d& b, 
                     const Point3d& c, const Point3d& d);

};//class Point3d


/************************************************************
 *	Inline implementation and some 3d predicates
 ************************************************************/
// removed inline implementation for compile under visual c++
// Zilin Du

// midPt(p, q) returns (p+q)/2:
Point3d midPt3d ( Point3d& a, Point3d& b);

/* orientation3d(a, b, c, d) 
 *   computes the orientation of points a, b, c, d as the sign
 *   of the determinant
 *              | ax  ay  az 1 |
 *              | bx  by  bz 1 |
 *              | cx  cy  cz 1 |
 *              | dx  dy  dz 1 |
 *   i.e., it returns +1 if d lies in the opposite side w.r.t. the 
 *   counter-clockwise side of plane formed by a, b, c
 */
int orientation3d(const Point3d& a, const Point3d& b, 
                         const Point3d& c, const Point3d& d); 

/* area(a, b, c) returns 1/2 times the determinant of orientation(a,b,c)
 * above.  This is the signed area of the triangle determined by a, b, c,
 * positive if orientation(a,b,c) > 0, and negative otherwise.  */

double volume(const Point3d& a, const Point3d& b, 
                     const Point3d& c, const Point3d& d);


/* returns true if points a, b, c and d are coplanar, i.e.,
 * orientation(a, b, c, d) = 0, and false otherwise. 
 */
bool coplanar(const Point3d& a, const Point3d& b, 
                     const Point3d& c, const Point3d& d);

/************************************************************
 *  CONSTANTS 
 ************************************************************/

static Point3d ORIGIN_3D(0.0, 0.0, 0.0);
static Point3d X_UNIT_3D(1.0, 0.0, 0.0);
static Point3d Y_UNIT_3D(0.0, 1.0, 0.0);
static Point3d Z_UNIT_3D(0.0, 0.0, 1.0);

#endif
