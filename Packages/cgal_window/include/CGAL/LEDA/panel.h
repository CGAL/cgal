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
// file          : include/CGAL/LEDA/panel.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


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
