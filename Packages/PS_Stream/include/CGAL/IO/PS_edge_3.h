// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/IO/PS_edge_3.h
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_PS_EDGE_3_H
#define CGAL_PS_EDGE_3_H

#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Point_3.h>
#include <CGAL/leda_real.h>
#include <CGAL/Aff_transformation_3.h>

CGAL_BEGIN_NAMESPACE

enum VISIBILITY{INVISIBLE,VISIBLE};

typedef double coord_type;
//typedef leda_real coord_type;

typedef CGAL::Cartesian< coord_type > CT;
typedef CGAL::Point_3< CT > Point3;
typedef CGAL::Aff_transformation_3< CT > Transformation;

class PS_edge_3 {

public:

  //Constructor
  PS_edge_3(Point3 &PA,Point3 &PB,Color color=BLACK,VISIBILITY v=VISIBLE);

  // Copy constructor
  PS_edge_3(const PS_edge_3& ar);

  // Inline functions

  //Accessors
  //Return the visibility of the edge 
  inline VISIBILITY visibility() const {return _visibility;};
  //Return the color of the edge
  inline Color color() const {return _edge_color;};
  //Return the first point of the edge
  inline Point3 first_point() const { return _first_point;};
  //Return the second point of the edge
  inline Point3 second_point() const { return _second_point;};
  

  //Settings
  //Modify the visibility of the edge
  inline void set_visibility(VISIBILITY v) {_visibility=v;};
  //Modify the color of the edge
  inline void set_color(Color c) {_edge_color=c;}
  //Modify the first point coordinate of the edge
  inline void set_first_point(Point3 pa) {_first_point=pa;};
  //Modify the second point coordinate of the edge
  inline Point3 set_second_point(Point3 pb) {_second_point=pb;};
  
  
  //Operators
  //Setting operator
  PS_edge_3 operator= (const PS_edge_3& ar);

  //Friend Function to display the parameters of the PS_Edge_3
  friend ostream& operator<< (ostream& os, const PS_edge_3& ar);
 
  //Others functions
  //Return the minimum or the maximun of the two points of the edge
  //according the x,y,z 
  coord_type xmin(); coord_type ymin(); coord_type zmin();
  coord_type xmax(); coord_type ymax(); coord_type zmax();
  //Function to transform the two points of the edge according the
  //Transformation t
  void transformation(Transformation &t);

private:
 
  //First point of the edge
  Point3 _first_point;
  //Second point of the edge
  Point3 _second_point;
  //Edge color
  Color _edge_color;
  //VISIBILITY of the edge :  0 for INVISIBLE, 1 for VISIBLE
  VISIBILITY _visibility;
};

CGAL_END_NAMESPACE

#endif // CGAL_PS_EDGE_3_H
