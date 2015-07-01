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
 * File: circle2d.h
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

#ifndef _CIRCLE2D_H
#define _CIRCLE2D_H

#include <CGAL/CORE/geom2d/point2d.h>
#include <CGAL/CORE/geom2d/line2d.h>

class Circle2d : public GeomObj {

   /* An instance C of the data type circle is an oriented circle 
      in the plane passing through three points p1, p2, p3. The 
      orientation of C is equal to the orientation of the three defining
      points. i.e. orientation(p1, p2, p3).
      If \Labs{\{p1, p2, p3\}} = 1, C is the empty circle with center p1.
      If p1, p2 and p3 are collinear, C is a straight line passing through
      p1, p2 and p3 in this order and the center of C is undefined.    
    */

private:
 
   Point2d p1;  // the 3 points defining the circle
   Point2d p2;
   Point2d p3;

   int orient;  //orientation(p1, p2, p3)

   Point2d* cp;  //pointer to center
   double * rp;  //pointer to radius

public:

   Circle2d( const Point2d& p1, const Point2d& p2,  const Point2d& p3);
   //initialized to the oriented circle through points p1, p2, p3 

   Circle2d(const Point2d& a, const Point2d& b0);
   //initialized to the counter-clockwise oriented circle with center a 
   //passing through b0 

   Circle2d(const Point2d& p);
   //initialized to the trivial circle with center p

   Circle2d();
   //initialized to the trivial circle with center (0,0)

   Circle2d(const Point2d& c, double r);
   //initialized to the circle with center c and radius r with positive
   //(i.e. counter-clockwise) orientation

   Circle2d(const Circle2d& c);
   //copy constructor 

   virtual ~Circle2d(); 

   Circle2d& operator=(const Circle2d& C);

   //operations

   Point2d center();
   //return the center of the circle

   double radius();
   //returns the radius. 
   //precond: the orientation of the circle is not 0

   Point2d point1() const { return p1; }
   Point2d point2() const { return p2; }
   Point2d point3() const { return p3; }

//   Point2d point_on_circle(float alpha);
   //returns a point p on the circle with angle of alpha

   bool is_degerate() const { return orient == 0; }
   //returns true if the defining points are collinear

   bool is_trivial() const {return p1 == p2; }
   //returns true if radius is zero

   int orientation() const { return orient; }

   int side_of(const Point2d& p) const; 
   // returns -1, +1 or 0 if p lies right of, left of or on the circle
   // respectively

   bool inside(const Point2d& p);
   //returns true if p lies inside of the circle

   bool outside(const Point2d& p); 
   
   bool contains(const Point2d& p) const ;
   //returns true if p lies on the circle, false otherwise

   double distance(const Point2d& p);
   //returns the distance between p and the circle: distance to center - radius

   double distance(const Line2d& l);
   //returns the distance between l and the circle
   //distance from center to l minus radius

   double distance(Circle2d& D);
   //returns the distance between this circle and circle D
   //distance between two centers minus two radius

   bool operator==(const Circle2d& D) const ;
   bool operator!=(const Circle2d& D) const 
	{ return !operator==(D); }

   friend std::ostream& operator<<(std::ostream& out, Circle2d& c);
   friend std::istream& operator>>(std::istream& in, Circle2d c); //?? Circle2d &

}; // class Circle2d


#endif
