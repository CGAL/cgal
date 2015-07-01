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
 * File: point2d.h
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


#ifndef _POINT2D_H
#define _POINT2D_H

#ifndef CORE_LEVEL
#  define CORE_LEVEL 3
#endif

#include <CGAL/CORE/CORE.h>
#include <CGAL/CORE/linearAlgebra.h>
#include <CGAL/CORE/geombase.h>

class Point2d : public GeomObj {

private:
  double x, y;

public:
 
  //CONSTRUCTORS
  //
  Point2d();  //initialized to origin(0,0)
  Point2d(double, double);
  Point2d(const Point2d &);
  Point2d(Vector v);
  //create a point initialized to the point $(v[0], v[1])$ 
  //precondition: v.dim() = 2

  //DESTRUCTOR
  virtual ~Point2d() {}

  //ASSIGNMENT AND QUERY
  //
  Point2d& operator=(const Point2d&);
  
  double X() const { return x; }
  double Y() const { return y; } 
  void setX( const double a){ x = a; }
  void setY( const double a){ y = a; } 
  void set( const double a, const double b){ x = a; y = b;} 
  
  int dim() const { return 2; }

  //CONVERSION
  //
  Vector toVector() const { return Vector(X(), Y()); } 

  //DISTANCES
  //
  double distance(const Point2d) const;
  // returns the Euclidean distance between p and this

  double distance() const { return distance(Point2d(0, 0)); }
  // returns distance between this and origin

  //VECTOR OPERATIONS
  //
  Vector operator-(const Point2d &) const;
  Point2d operator+(const Vector &) const;

  //TRANSFORMATIONS
  //
  Point2d rotate90( const Point2d& q);
  // returns the point rotated about q by angle of 90 degrees

  //COMPARISONS
  //
  bool operator==(const Point2d&) const;
  bool operator!=(const Point2d& p) const {return !operator==(p); }
  
  //INPUT-OUTPUT
  //
  friend std::ostream& operator<< (std::ostream&, const Point2d);
  // write point p to output stream
  // The format is, e.g.,  Point2d(1.0, 23)

  friend std::istream& operator>> (std::istream&, Point2d&);
  // reads the x and y coordinates of point p from the input stream
  // The format is " ( x , y ) " where the white spaces are optional.
  // Even the comma and the "(" and ")" are optional.
  // The comment char '#' is allowed, and the rest of
  // the line is then treated as white space.
  // However, you must not use other kinds of parenthesis
  // E.g., the following are all equivalent:
  // 	1.0 -0.2 # comment
  // 	( +1.0, -0.2)
  // 	1.0, -.2 

  friend int readPoints(std::istream &iS,
		  Point2d *pA, int MaxN = 1000, int N = 0);
  // reads a sequence of points from input stream iS into Point2d array pA.
  // The input stream constains a sequence of 2K+1 numbers of the form
  //     [NN]   ( x1 ,  y1 )  ( x2 , y2 )  ... ( xK , yK )
  // The i-th point is (xi, yi). 
  //    (0) NN is optional if N is given as argument (then N is set to NN)
  //    (1) Any of the "(", "," and ")" are optional
  //    (2) Newlines, extra white spaces, '#' are all ignored.
  //    (3) Also, everything after '#' is discarded.
  // If N > MaxN, nothing is read and 0 is returned.  
  // Returns the number of points actually read, i.e, min(K, N).
  
}; //Point2d Class

// //////////////////////////////////////////////////
// AUXILLIARY FUNCTIONS:
// //////////////////////////////////////////////////

Point2d midPoint(const Point2d& a, const Point2d& b); 
// returns midpoint between a and b

Point2d aCenter(const Point2d& a, const Point2d& b, machine_double alpha =0.5); 
// returns the "asymmetric Center" point 
// that is alpha-fraction of the distance from a to b

double area(const Point2d& a, const Point2d& b, const Point2d& c);
  // returns twice (!) the signed area of triangle (a,b,c)
  
int orientation2d(const Point2d& a, const Point2d& b, const Point2d& c);
  // returns sign of area(a,b,c)
  
bool leftTurn(const Point2d& a, const Point2d& b, const Point2d& c);
  // returns true iff orientation2d(a,b,c) = +1
  
bool rightTurn(const Point2d& a, const Point2d& b, const Point2d& c);
  // returns true iff orientation2d(a,b,c) = -1

bool collinear(const Point2d& a, const Point2d& b, const Point2d& c);
  // returns true iff orientation2d(a,b,c) = 0
  
bool between(const Point2d& a, const Point2d& b, const Point2d& c);
  // returns true iff orientation2d(a,b,c) = 0 and b is strictly
  //    between a and c.

//variant of between:
bool betweenVar(const Point2d& a, const Point2d& b, const Point2d& c);
  // returns true iff the scalar product (a-b, c-b) is positive.
  //    In case orientation2d(a,b,c)=0, then this is equivalent to
  //    b being strictly between a and c.

// THE FOLLOWING ARE CALLED by
// 	operator>>(..) and readPoints(..)
// bool getToNum( std::istream& in, char mark, bool strict=false) ;
// bool getToChar( std::istream& in, char mark) ;
// bool startNum(char c) ;

#endif
