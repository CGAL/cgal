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



  void set_color1(Color c) {color=c;}
  void set_color2(Color c) {col2=c;}
  virtual void add_point(double x, double y, double z) {}

  //POUR le postscript : mettre la bonne signature
  //  virtual void to_ps(Ps_stream &ps){}


};

CGAL_END_NAMESPACE
