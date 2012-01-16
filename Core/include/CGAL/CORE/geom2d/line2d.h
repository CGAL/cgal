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
 * File: line2d.h
 * Synopsis:
 *      Basic 2-dimensional geometry
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

#ifndef _LINE2D_H_
#define _LINE2D_H_

#include <CGAL/CORE/geom2d/point2d.h>

class Line2d : public GeomObj {

  /* An instance l of the data type $line$ is a directed straight line
     in the two dimensional plane. The angle between a right oriented 
     horizontal line and $l$ is called the direction of $l$.
   */
  /* member Vector is not used in this class, it's intended for use in
       the operator+,- etc 
     need to do: assure p0 != p1
   */
private:

  Point2d p0;
  Point2d p1;
  Vector V;

public:

  /*************************************************************
   *  constructors
   *************************************************************/

  Line2d(const Point2d & p, const Vector &v);
  // line initialized to pass through points p and p+v

  Line2d(const Point2d &p, const Point2d &q);
  //line is initialized to pass through points p and q directed from p to q

//  Line2d(const point& p, double alpha);
  //line passes through point p with direction alpha

  Line2d(const double& a, const double& b, const double& c);

  Line2d(const Line2d &);

  Line2d();
  //line passes through the origin with direction 0.

  virtual ~Line2d() {}
  /*************************************************************
   *   member functions
   *************************************************************/

  Vector direction() const { return p1-p0; }
   // returns the direction as a vector 

  Point2d startPt() const { return p0; }  
  Point2d stopPt() const  { return p1; }

  double distance(Point2d q) const;
   // returns the Euclidean distance between this line and point q

  Point2d projection(const Point2d& p) const;
   // returns the projection of p on this line

  int orientation( const Point2d& p ) const;
   // orientation of p0, p1 and p

   // the sine/cosine of the angle made with positive x-direction
  double sine() const { return (p1.Y() - p0.Y()) / p0.distance(p1); }
  double cosine() const { return (p1.X() - p0.X()) / p0.distance(p1); }

  Line2d rotate90( const Point2d& q)
  { return Line2d(startPt().rotate90(q), stopPt().rotate90(q)); }

  double y_abs() const;
  // returns the y-abscissa of the line 

  double slope() const ;    
  //precond: is not vertical

  /*************************************************************
   *   predicates
   *************************************************************/

  bool isVertical() const { return p0.X() == p1.X(); }
  bool isHorizontal() const { return p0.Y() == p1.Y(); }
  bool is_trivial() const {return p0 == p1; }   //meaning for a line?

  bool contains( const Point2d& p) const { 
          return orientation2d(p0, p1, p) == 0; }
  bool isCoincident( const Line2d& g) const { 
          return contains(g.p0) && contains(g.p1); }  
  bool isParallel(const Line2d& l) const {
    return det(V, l.direction()) == 0; }

  bool operator==( const Line2d& g ) const { return isCoincident(g); }
  bool operator!=( const Line2d& g ) const { return !operator==(g); }

  /*************************************************************
   *   intersection
   *************************************************************/

  int intersects(const Line2d& t) const;
   // decides whether *this and t intersects
   // return dim of intersection. 
   // return -1 if no intersection

  GeomObj* intersection(const Line2d& g) const;
   //if this line and g intersect in a single point, this point is 
   // assigned to p and the result is true, otherwise the result is false

  /*************************************************************
   *   angles and others
   *************************************************************/
  friend int orientation2d( const Line2d& l, const Point2d& p);
  // computes the orientation (a, b, p), where a!=b and a and b appear 
  // in this order on line l

  friend int cmp_slopes(const Line2d& l1, const Line2d& l2) 
  //l1.slope > l2.slope: +1; equal: 0; otherwise: -1
  {
     if (l1.slope() == l2.slope())
         return 0;
     else
         return (l1.slope() > l2.slope()) ? +1 : -1;
  } 
  
  /*************************************************************
   *   I/O
   *************************************************************/

  friend std::istream& operator>>(std::istream& in, Line2d& l);
  friend std::ostream &operator<<(std::ostream & out, const Line2d & l);
}; // class Line2d

extern Line2d p_bisector(const Point2d& p, const Point2d& q);

#endif
