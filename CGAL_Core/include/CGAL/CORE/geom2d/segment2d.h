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
 * File: segment2d.h
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


#ifndef _SEGMENT2D_H_
#define _SEGMENT2D_H_

#include <CGAL/CORE/geom2d/point2d.h>
#include <CGAL/CORE/geom2d/line2d.h>

/************************************************************
 *  Class Segment2d:
 *
 *  An instance s of Segment2d is a finite or infinite line segment
 *  in the two dimensional plane, defined by a start point
 *  s.startPt() and a stop point s.stopPt().  It can be regarded
 *  as an open or a closed segment (default: open), and directed
 *  or not (default: directed).
 *
 *  We do not necessarily assume that startPt() != stopPt().
 ************************************************************/

class Segment2d : public GeomObj {

private:

  Point2d p0;
  Point2d p1;
  bool directed;  // segments can be directed or not (default is directed)
  bool open;     // segments can be open or closed (default is open)

public:

  /************************************************************
   *   constructors
   ************************************************************/

  Segment2d(const Segment2d &);

  Segment2d(const Point2d &p, const Point2d &q);
  	//finite segment with endpoints p and q

  Segment2d(const Point2d & p, const Vector & v);
	//ray segment

  Segment2d();
  	//unit segment from (0,0) to (1,0) 

  virtual ~Segment2d() {}
  /*************************************************************
   *   member functions
   *************************************************************/

  Point2d startPt() const { return p0; }  
  Point2d stopPt() const { return p1; }

  void reverse() { Point2d pTmp = p0; p0=p1; p1=pTmp; }
   //  reverses the direction of the segment

  Line2d toLine() const { return Line2d(p0,p1); }
 
  double length() const { return p0.distance(p1); }
    //length of segment

  double distance( const Point2d& p ) const;
   // returns the Euclidean distance between this segment and point q

  Point2d nearPt( const Point2d& p ) const;
   // returns the point on segment closest to q;

  void setDirected( bool _directed ) { directed = _directed; }
  void setOpen( bool _open ) { directed = open; }

   // orientation of p0, p1 and p
  int orientation( const Point2d& p ) const {
    return toLine().orientation(p); }

  /*************************************************************
   *   predicates
   *************************************************************/
  bool isOpen() const {return open; }
  bool isDirected() const {return directed; }
  bool isTrivial() const {return p0 == p1; }  
  bool isVertical() const { return p0.X() == p1.X(); }
  bool isHorizontal() const { return p0.Y() == p1.Y(); }
  bool isCollinear( Point2d& p ) const { return toLine().contains(p); }
  bool isCoincident( const Segment2d& s) const;
  bool isParallel( const Segment2d& s ) {
    return toLine().isParallel( s.toLine() );  }

  bool contains( const Point2d& p ) const;
  bool contains( const Segment2d& s ) const { 
     return contains(s.startPt()) && contains(s.stopPt()); }  

  bool operator==(const Segment2d& s) const { return isCoincident(s); }

  bool operator!=(const Segment2d& s) const { return !isCoincident(s); }

  /*************************************************************
   *   intersection
   *************************************************************/

  int intersects( const Line2d& l ) const;
   //decides whether *this and t intersect in one point
   // return dim of intersetion
  
  int intersects( const Segment2d& s ) const;
   //decides whether *this and t intersect in one point
   // return dim of intersetion

  GeomObj* intersection( const Line2d& l ) const;
   // return intersection point if this segment and l intersect at a single point
   // the intersection point is returned 
  
  GeomObj* intersection( const Segment2d& s ) const;
   // return intersection point if this segment and s intersect at a single point
   // the intersection point is returned 
 
  /*************************************************************
   *   angles
   *************************************************************/

   // the sine/cosine of the angle made with positive x-direction
  double sine() const { return (p1.Y() - p0.Y()) / length() ; }
  double cosine() const { return (p1.X() - p0.X()) / length() ; }

  Line2d rotate90(const Point2d& q)
  { return Line2d(startPt().rotate90(q), stopPt().rotate90(q)); }

   // computes the orientation (a, b, p), where a!=b and a and b appear 
   // in this order on segment l
  friend int orientation2d( const Segment2d& s, const Point2d& p) {
    return orientation2d( s.toLine(), p ); 
  }

  friend int cmp_slopes( const Segment2d& s1, const Segment2d& s2) 
   //l1.slope > l2.slope: +1; equal: 0; otherwise: -1
  {
     Line2d l1 = s1.toLine();
     Line2d l2 = s2.toLine();
     if (l1.slope() == l2.slope())
         return 0;
     else
         return (l1.slope() > l2.slope()) ? +1 : -1;
  } 
  
  /*************************************************************
   *   I/O
   *************************************************************/


  friend std::istream& operator>>(std::istream& in, Segment2d& l);
  friend std::ostream &operator<<(std::ostream & out, const Segment2d & l);
   // syntax: {[} p {===} q {]}

}; //class Segment2d

#endif
