// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Drawable_object_3.h
// revision      : $Revision$
//
// author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef DRAWABLE_OBJECT_H
#define DRAWABLE_OBJECT_H

#include "PS_Stream_3.C"

typedef CGAL::Cartesian<double> D;
typedef CGAL::Bbox_3 PS_BBox3;
typedef CGAL::Direction_3< D > Direction;
typedef CGAL::Point_3< D > Point3;
typedef CGAL::Point_2< D > Point2;
typedef CGAL::Plane_3< D > Plane3;typedef CGAL::Direction_3< D > Direction;
typedef CGAL::Line_3< D > Line3;


#include <CGAL/basic.h>
#include <CGAL/IO/Color.h>
#include <GL/gl.h>
#define DRAWABLE
//enum Style {FILL=1, WIRE, RAW, FOS1, FOS2, FOS3, FOS4, FOS5};
CGAL_BEGIN_NAMESPACE

enum Style {FILL=1, WIRE, RAW, FOS1, FOS2, FOS3, FOS4, FOS5};
typedef int           Size;
typedef unsigned char Precision;



//############################################################################
//#### GENERIC CLASS FOR DRAWABLE OBJECT #########
class Drawable_object {

protected:
  Color color;
  Color col2;
  Style style;
  Precision precision;
  Size     size;
  double o_center[3];
  int lind;

public:
  char* type;
  // Drawable_object(){type="Undefined";}

  Drawable_object()
     {}

  virtual ~Drawable_object(){}
  virtual void draw()  {std::cerr << "virtual draw object()" << std::endl;}
  void set_center() {std::cerr << "virtual set_center" <<
				     std::endl;}
  double get_center(int i) 
    {
      if ((i<1) || (i>3))
        std::cerr << "bad indice value : " << i <<std::endl;
      else if (i==1)
        return o_center[0];
      else if (i==2)
          return o_center[1];
      else
        return o_center[2];
      return 0;
    }

  void set_style(Style s) 
    {
      style=s; 
      glDeleteLists(lind,1);
      lind=0;
    }

  void set_colors(Color c1, Color c2) 
    {
      color=c1; 
      col2=c2; 
      glDeleteLists(lind,1);
      lind=0;
    }

  void set_color1(Color c) {set_colors(c, col2);}
  void set_color2(Color c) {set_colors(color, c);}

  virtual void add_point(double x, double y, double z) {}

  virtual void to_ps(PS_Stream_3 &ps){
  }


};

CGAL_END_NAMESPACE

#endif
