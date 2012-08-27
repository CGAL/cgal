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
 * File: polygon3d.h
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

#ifndef _POLYGON3D_H_
#define _POLYGON3D_H_

#include <CGAL/CORE/geom3d/point3d.h>
#include <CGAL/CORE/geom3d/segment3d.h>
#include <CGAL/CORE/geom3d/plane3d.h>
#include <CGAL/CORE/geom3d/line3d.h>

class Polygon3d : public GeomObj {

private:
  class PointNode {
    public:
      Point3d  *p;
      PointNode *prev;
      PointNode *next;

      PointNode( const Point3d& _p ) { p = new Point3d(_p); prev=NULL; next=NULL; }
      ~PointNode()  { delete p; }
  };

  PointNode* headN;       //head of linked list
  int size;               //length of the list

public:
  class Iterator;
  friend class Iterator;
  class Iterator {
    private:
      PointNode* pointer;

    public:
      Iterator( PointNode* node ) { pointer = node; }

      PointNode* getPointer() { return pointer; }
      //bool hasNext() { return pointer->next != Polygon3d::headN; }
      //Point3d* nextPoint() { pointer=pointer->next; return pointer->prev->p; }

      Iterator& operator =(Iterator & it) { pointer=it.getPointer(); return *this; }

       // postfix only
      Iterator& operator ++(int) { pointer = pointer->next; return *this; }
      Iterator& operator --(int) { pointer = pointer->prev; return *this; } 
      Point3d* getPoint() { return pointer->p; }

       // remove current node
      void remove( )      { 
          PointNode* temp = pointer;
          pointer->next->prev = pointer->prev;
          pointer->prev->next = pointer->next;
          pointer = pointer->prev;
          delete temp;
      }

  };

   // default
  Polygon3d();

  // initialize given a triangle
  // not included, use Polygon3d::toPolygon() instead
  //Polygon3d(const Triangle3d& T);
  


 //copy constructor
  Polygon3d(const Polygon3d& plg);

  virtual ~Polygon3d();

  /************************************************************
   *   member functions
   ************************************************************/

   // view a polygon as a surface
  virtual int dim() const { return 2; }
   
  int getSize() const { return size; }

  Iterator getIterator() { return Iterator(headN); }

  bool searchPoint( const Point3d& p ) const;
  Point3d* getPoint( int index ) const;
   // append to the end of list
  void addPoint( const Point3d& p );  

   // remove point at given position
  void removePoint( int index );
   // return false if point doesn't exist
  bool removePoint( const Point3d& p );
  void removeAllPoints();
  
   // get next point of p
   // return NULL if p doesn't exist
  Point3d* nextPoint( const Point3d& p ) const;

   // get previous point of p
   // return NULL if p doesn't exist
  Point3d* prevPoint( const Point3d& p ) const;

   //operators
  Point3d* operator []( int index ) const { return getPoint(index); }
  Polygon3d& operator =(const Polygon3d& plg);

   // test identity
  bool operator ==(Polygon3d& plg) const;
  bool operator !=(Polygon3d& plg) const  { return !(*this == plg); }

   // verify if polygon is valid
  bool verify();
   
   // test coplanarity 
  bool isCoplanar( const Point3d& p ) const;
  bool isCoplanar( const Segment3d& s )  const { 
         return isCoplanar(s.startPt()) && isCoplanar(s.stopPt()); }
  bool isCoplanar( const Line3d& l )     const { 
         return isCoplanar(l.startPt()) && isCoplanar(l.stopPt()); }

  bool isCoplanar( const Plane3d& pl ) const;

   // test if p is on any edge
  bool isOnEdge( const Point3d& p ) const;
 
protected:  
  void freeMemory();      // delete the point list

   // copy point list from other polygon
   // make sure old list has been deleted before calling this
  void copy( const Polygon3d& plg );

  /************************************************************
   *   I/O 
   ************************************************************/

  friend std::ostream& operator<<(std::ostream& in, const Polygon3d& plg);
  
}; //class Polygon3d

#endif
