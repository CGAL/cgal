// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
//
// file          : include/CGAL/LEDA/color.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================

#ifndef CGAL_WINDOW_COLOR_H
#define CGAL_WINDOW_COLOR_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <CGAL/LEDA/system.h>
#include <CGAL/LEDA/basic.h>
#include <CGAL/LEDA/impl/x_window.h>

/*{\Manpage {color} {} {Colors}}*/


namespace CGAL {

class __exportC color {

/*{\Mdefinition
The data type $color$ is the type of all colors available for
drawing operations in windows (cf. \ref{Windows}).
Each color is defined by a triple of integers $(r,g,b)$ with
$0\le r,g,b \le 255$, the so-called {\em rgb-value} of the color.
The number of available colors is restricted and depends on the 
underlying hardware. Colors can be created from rgb-values,
from names in a color data base (X11), or from the 16 integer 
constants (enumeration in \<CGAL/LEDA/incl/impl/x\_window.h\>)
$black$, $white$, $red$, $green$, $blue$, $yellow$, $violet$, $orange$;
$cyan$, $brown$, $pink$, $green2$, $blue2$, $grey1$, $grey2$, $grey3$. 
}*/

  int  col_index;
  bool ok;

public:

/*{\Mcreation col}*/

  color();
/*{\Mcreate  creates a color with rgb-value $(0,0,0)$ (i.e. black). }*/

  color(int r, int g, int b);
/*{\Mcreate  creates a color with rgb-value $(r,g,b)$. }*/

  color(const char* name);
/*{\Mcreate  creates a color and initializes it with the rgb-value
             of color $name$ from the X11 color data base ($(0,0,0)$ 
             if $name$ does not exist). }*/

  color(int i);
/*{\Mcreate  creates a color and initializes it with one of the
             16 predefined colors. Here $i$ is one of the 16 integer
             constants $black$, $white$, $red$, $green$, $blue$, 
             $yellow$, $violet$, $orange$, $cyan$, $brown$, $pink$,                          $green2$, $blue2$, $grey1$, $grey2$, or $grey3$.  }*/


/*{\Moperations 1.5 6}*/

  void set_rgb(int r, int g, int b);
/*{\Mop    sets the rgb-value of |\Mvar| to $(r,g,b)$. }*/

  void get_rgb(int& r, int& g, int& b) const;
/*{\Mop    assigns the rgb-value of |\Mvar| to $(r,g,b)$. }*/

  void set_red(int r);
/*{\Mop    changes the  r-value of |\Mvar| to $r$. }*/

  void set_green(int g);
/*{\Mop    changes the  g-value of |\Mvar| to $g$. }*/

  void set_blue(int b);
/*{\Mop    changes the  b-value of |\Mvar| to $b$. }*/

  bool is_good() const { return ok; }

  int  get_index() const { return col_index; }

  operator int() const  { return col_index; }


#if defined(__SUNPRO_CC) && (__SUNPRO_CC >= 0x500)

friend bool operator==(const color& c1 , const color& c2)
{ return c1.col_index == c2.col_index; }

friend bool operator!=(const color& c1 , const color& c2)
{ return c1.col_index != c2.col_index; }

#endif

};

extern __exportF std::istream& operator>>(std::istream& in, color& c);
extern __exportF std::ostream& operator<<(std::ostream& out, const color& c);

}

#endif
