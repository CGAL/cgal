// Copyright (c) 1999  Martin-Luther-University Halle-Wittenberg (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Matthias Baesken, Algorithmic Solutions


#ifndef CGAL_WINDOW_PANEL_H 
#define CGAL_WINDOW_PANEL_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <CGAL/LEDA/window.h>

//------------------------------------------------------------------------------
//   PANELS
//------------------------------------------------------------------------------

/*{\Manpage {panel} {} {Panels}}*/


namespace CGAL {


class __exportC panel : public window {

/*{\Mdefinition
Panels are windows consisting of a panel section only (cf. section
\ref{Windows}). They are used for displaying text messages and updating the
values of variables.
}*/

public:

/*{\Mcreation P }*/

 panel() : window(-1,-1) {}
/*{\Mcreate creates an empty panel $P$.}*/

 panel(string s) : window(-1,-1,s.c_str()) {}
/*{\Mcreate creates an empty panel $P$ with header $s$.}*/

 panel(int w, int h) : window(w,h) {}
/*{\Mcreate creates an empty panel $P$ of width $w$ and height $h$.}*/


 panel(int w, int h,string s) : window(w,h,s.c_str()) {}
/*{\Mcreate creates an empty panel $P$ of width $w$ and 
            height $h$ with header $s$.}*/

// for backward compatibility
panel(string s, int w, int h) : window(w,h,s.c_str()) {}


~panel() {}


/*{\Moperations 1.2 4.8 }*/

/*{\Mtext       
All window operations for displaying, reading, closing and adding panel 
items are available (see section \ref{panel-operations}). There are two 
additional operations for opening and reading panels.
}*/

// open = display + read + close

int  open(int x = window::center, int y = window::center)
{ return panel_open(x,y,0); } 
/*{\Mopl   |P.display(x,y)| $+$ |P.read_mouse()| $+$ |P.close()|.}*/


int open(window& W, int x=window::center, int y=window::center)
{ return panel_open(x,y,&W); } 
/*{\Mopl   |P.display(W,x,y)| $+$ |P.read_mouse()| $+$ |P.close()|.}*/


};

}


#endif
