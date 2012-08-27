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
 * File: segment3d.h
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

#ifndef _SEGMENT3D_H_
#define _SEGMENT3D_H_

#include <CGAL/CORE/geom3d/point3d.h>
#include <CGAL/CORE/geom3d/line3d.h>

/************************************************************
 *  Class Segment3d:
 *
 *  An instance s of Segment3d is a finite or infinite line segment
 *  in the three dimensional plane, defined by a start point
 *  s.startPt() and a stop point s.stopPt().  It can be regarded
 *  as an open or a closed segment (default: open), and directed
 *  or not (default: directed).
 *
 *  We do not necessarily assume that startPt() != stopPt().
 ************************************************************/

#define S_TYPE_FINITE 0
#define S_TYPE_RAY    1
#define S_TYPE_LINE   2

class Plane3d;

class Segment3d : public GeomObj{

private:

  Point3d p0;
  Point3d p1;
  bool directed;     // segments can be directed or not (default is directed)
  bool open;        // segments can be open or closed (default is open)
  //int finite;     // 0=finite, 1=ray, 2=line (default is 0)

public:

  /************************************************************
   *   constructors
   ************************************************************/

  Segment3d(const Segment3d &s);

  Segment3d(const Point3d &p, const Point3d &q);
  	//finite segment with endpoints p and q

  Segment3d(const Point3d & p, const Vector & v);
	//ray segment

  Segment3d();
  	//trivial segment from (0,0) to (0,0) 

  virtual ~Segment3d() {}
  /*************************************************************
   *   member functions
   *************************************************************/
  
  virtual int dim() const { return 1; }

  Point3d startPt() const { return p0; }  
  Point3d stopPt()  const { return p1; }
  Vector direction() const { return p1 - p0; }

  void reverse();
  //  reverses the direction of the segment

  void setDirected( bool beDirected ) { directed = beDirected; }
  void setOpen(  bool beOpen ) { open = beOpen; }
  
  void setStartPt( Point3d& p ) { p0 = p; }
  void setStopPt ( Point3d& p ) { p1 = p; }

  double length() const { return p0.distance(p1); }
    //length of segment
    
  Line3d toLine() const { return Line3d(p0,p1); }
 
  double distance( const Point3d& p ) const;
  // returns the Euclidean distance between this segment and point q

  Point3d nearPt( const Point3d& q ) const;
   // returns the point on segment closest to q;

  /*************************************************************
   *   predicates
   *************************************************************/

  bool isDirected() const { return directed; }
  bool isOpen() const {return open; }
  bool isTrivial() const {return p0 == p1; }  
  bool isCollinear( const Point3d& p ) const {return toLine().contains(p); }
  bool contains( const Segment3d& s ) const { return contains(s.startPt()) && contains(s.stopPt()); }
  bool isCoincident( const Segment3d& s) const;
  
  bool isCoplanar( const Line3d& s) const;
  bool isCoplanar( const Segment3d& s) const;
  
  bool contains( const Point3d& p ) const;
  
  bool operator==( const Segment3d& s ) { return isCoincident( s ); }

  bool operator!=( const Segment3d& s ) { return !operator==(s); }

  /*************************************************************
   *   intersection
   *************************************************************/

  int intersects( const Line3d& l ) const;
  //decides whether *this and t intersect in one point
  // return dim of intersetion
  
  int intersects( const Segment3d& s ) const;
  //decides whether *this and t intersect in one point
  // return dim of intersetion

  GeomObj* intersection( const Line3d& l ) const;
  // return intersection point if this segment and l intersect at a single point
  // the intersection point is returned 
  
  GeomObj* intersection( const Segment3d& s ) const;
  // return intersection point if this segment and s intersect at a single point
  // the intersection point is returned 
 
  Plane3d bisect_plane() const;
   // return bisector plane
   
  /*************************************************************
   *   I/O
   *************************************************************/

  friend std::istream& operator>>(std::istream& in, Segment3d& l);
  friend std::ostream& operator<<(std::ostream& out, const Segment3d& l);

}; //class Segment3d


#endif
