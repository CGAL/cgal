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

#include "CORE.h"
#include "../linearAlgebra.h"
#include "../geombase.h"

class Point2d : public GeomObj {

private:
  double x, y;

public:
 
  //constructors
  Point2d();  //initialized to origin(0,0)
  Point2d(double, double);
  Point2d(const Point2d &);
  Point2d(Vector v);
  //create a point initialized to the point $(v[0], v[1])$ 
  //precondition: v.dim() = 2

  //destrutor
  virtual ~Point2d() {}

  Point2d& operator=(const Point2d&);
  
  double X() const { return x; }
  double Y() const { return y; } 

  Vector toVector() const { return Vector(X(), Y()); } 

  int dim() const { return 2; }

  double distance(const Point2d) const;
  // returns the Euclidean distance between p and this

  double distance() const { return distance(Point2d(0, 0)); }
  // returns distance between this and origin

  Vector operator-(const Point2d &) const;
  Point2d operator+(const Vector &) const;

  Point2d rotate90( const Point2d& q);
  // returns the point rotated about q by angle of 90 degrees

  bool operator==(const Point2d&) const;
  bool operator!=(const Point2d& p) const {return !operator==(p); }
  
  friend std::ostream& operator<< (std::ostream&, const Point2d);
  // write point p to output stream

  friend std::istream& operator>>(std::istream&, Point2d&);
  // reads the x and y coordinates of point p from the input stream

};


Point2d center(Point2d& a, Point2d& b); 
int orientation2d(const Point2d& a, const Point2d& b, const Point2d& c);
double area(const Point2d& a, const Point2d& b, const Point2d& c);
bool collinear(const Point2d& a, const Point2d& b, const Point2d& c);

#endif
