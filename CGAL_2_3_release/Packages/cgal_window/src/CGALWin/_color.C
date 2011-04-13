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
// file          : src/CGALWin/_color.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================



#include <CGAL/LEDA/basic.h>
#include <CGAL/LEDA/color.h>
#include <CGAL/LEDA/impl/x_basic.h>


namespace CGAL {


color::color() : ok(true)
{ col_index = 0; }

color::color(const char* name) 
{ col_index = x_new_color(name);  
  if (col_index == -1) col_index = white;
}

color::color(int i ) : ok(true) { col_index = i;  }

color::color(int r, int g, int b) 
{ if (r <   0) r = 0;  
  if (r > 255) r = 255; 
  if (g <   0) g = 0;  
  if (g > 255) g = 255; 
  if (b <   0) b = 0;  
  if (b > 255) b = 255; 
  col_index = x_new_color(r,g,b); 
  ok = true;
  if (col_index == -1) 
  { col_index = white;
    ok = false;
   }
}


void color::get_rgb(int& r,int& g,int& b) const
{ x_get_rgb(col_index,&r,&g,&b); }


void color::set_rgb(int r, int g, int b)
{ ok = x_set_rgb(col_index,r,g,b) != 0; }


void color::set_red(int x)
{ int r,g,b;
  x_get_rgb(col_index,&r,&g,&b); 
  ok = x_set_rgb(col_index,x,g,b) != 0;
}

void color::set_green(int x)
{ int r,g,b;
  x_get_rgb(col_index,&r,&g,&b); 
  ok = x_set_rgb(col_index,r,x,b) != 0;
}

void color::set_blue(int x)
{ int r,g,b;
  x_get_rgb(col_index,&r,&g,&b); 
  ok = x_set_rgb(col_index,r,g,x) != 0;
}



std::ostream& operator<<(std::ostream& out, const color& c)
{ int r,g,b; 
  if (c == invisible)
    r = g = b = -1;
  else
    c.get_rgb(r,g,b);
  return out << r << " " << g << " " << b;
}

std::istream& operator>>(std::istream& in, color& c)
{ int r,g,b; 
  in >> r >> g >> b;
  if (r == -1)
     c = invisible;
  else
     c = color(r,g,b);
  return in;
}

}
