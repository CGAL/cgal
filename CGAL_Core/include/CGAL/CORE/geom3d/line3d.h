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
 * File: line3d.h
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
 
#ifndef _LINE3D_H_
#define _LINE3D_H_

#include <CGAL/CORE/geom3d/point3d.h>

class Line3d : public GeomObj{

 /************************************************************
  *  An instance l of the class Line3d is a directed line
  *  in the three dimensional plane. The angle between a right oriented 
  *  horizontal line and $l$ is called the direction of $l$.
  *  We assume that l is defined by two points, l.startPt()
  *  and l.stopPt().  We do not assure that these points are distinct.
  *  So the line could be "improper".  But most operators assume
  *  that lines are proper (it is the user's responsibility to check).
  *
  *  In the future, we may generalize this to allow the dual 
  *  representation in terms of a linear equation.
  *
  *  Member Vector is not used in this class. It is intended for use in
  *  the operator+/-, etc.
  ************************************************************/
private:

  Point3d p0; // = startPt
  Point3d p1; // = stopPt
  Vector V;   // = stopPt - startPt

public:

  /************************************************************
   *   constructors
   ************************************************************/

  Line3d(const Point3d & p, const Vector &v);
  // line initialized to pass through points p and p+v

  Line3d(const Point3d &p, const Point3d &q);
  //line is initialized to pass through points p and q directed from p to q

  Line3d(const Line3d &l);
  //copy constructor

  //Line3d(const Segment3d & s): p0(s.startPt()), p1(s.stopPt()), V(s.direction()) {}
  //construct from a segment

  Line3d();
  //horizontal line passes through the origin with direction 0.

  virtual ~Line3d() {}
  /************************************************************
   *   MEMBERS
   ************************************************************/

  virtual int dim() const { return 1; }

  Point3d startPt() const { return p0; }  
  Point3d stopPt() const { return p1; }
  const Vector direction() const { return V; }
  
  double distance(const Point3d& q) const;
  // returns the Euclidean distance between this line and point q
  
  Point3d projection(const Point3d& p) const;
  // returns the projection of p on this line

  /************************************************************
   *   PREDICATES
   ************************************************************/

  bool isProper() {return p0 == p1; }   
  bool contains(const Point3d& p) const;

  bool isCoincident(const Line3d& g) const;
  bool isParallel(const Line3d& g) const;  // same slope
  bool isSkew(const Line3d& l2) const;

  bool operator==(const Line3d& g) const { return isCoincident(g); }
  bool operator!=(const Line3d& g) { return !operator==(g); }

  int intersects(const Line3d& g ) const;
    // return the dimension of the intersection of g with this line:
    // 	-1 if disjoint (i.e., parallel but distinct lines)
    //  1  if coincident
    //  0  if intersect in a point.  In this case, the
    //		intersection point is assigned to p if this is available.

  GeomObj* intersection(const Line3d &l) const;

  /************************************************************
   *   I/O 
   ************************************************************/

  friend std::istream& operator>>(std::istream& in, Line3d& l);
  friend std::ostream& operator<<(std::ostream & out, const Line3d & l);

}; //class Line3d

#endif
