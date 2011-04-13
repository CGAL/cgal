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
// file          : include/CGAL/LEDA/menu.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


#ifndef CGAL_WINDOW_MENU_H 
#define CGAL_WINDOW_MENU_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <CGAL/LEDA/window.h>

//------------------------------------------------------------------------------
//   MENUES
//------------------------------------------------------------------------------

namespace CGAL {


/*{\Manpage {menu} {} {Menues}}*/

class __exportC menu : public window {

/*{\Mdefinition
Menues are special panels consisting only of a vertical list of buttons. }*/

public:

/*{\Mcreation M }*/

 menu() : window(-1,-1) { buttons_per_line(1); }
/*{\Mcreate creates an empty menu $M$.}*/

 menu(string s) : window(-1,-1,s.c_str()) { buttons_per_line(1); }
/*{\Xcreate creates an empty menu $M$ with header $s$.}*/

~menu() {}


/*{\Moperations 1.2 4.8 }*/

int button(string s, int n) { return window::button(s,n); }
/*{\Mopl     adds a button with label $s$ and number $n$ to $M$.}*/

int  button(string s) { return window::button(s); }
/*{\Mopl     adds a new button to $M$ with label $s$ and number equal to its
             position in the list of all buttons (starting with $0$).}*/


int  button(string s, int n, void (*F)(int)) { return window::button(s,n,F); }
/*{\Mopl     adds a button with label $s$, number $n$ and action 
             function $F$ to $M$. Function $F$ is called with actual
             parameter $n$  whenever the button is pressed. }*/

int  button(string s, void (*F)(int)) { return window::button(s,F); }
/*{\Mopl     adds a button with label $s$, number equal to its rank and action 
             function $F$ to $M$. Function $F$ is called with the number of the
             button as argument whenever the button is pressed. }*/

int  button(string s, int n, window& W) { return window::button(s,n,W); }
/*{\Mopl     adds a button with label $s$, number $n$, and attached window 
             $W$ to $M$.Whenever the button is pressed $W$ is opened. }*/


int  button(string s, window& W) { return window::button(s,W); }
/*{\Mopl     adds a button with label $s$ and attached window $W$ to $M$.
             Whenever the button is pressed $W$ is opened and 
             $W$.read\_mouse() is returned. }*/

void separator() { window::separator(); }
/*{\Mopl     inserts a separator (horizontal line) at the current position. }*/ 


int open(window& W, int x, int y) { return menu_open(x,y,&W); } 
/*{\Mopl   open and read menu $M$ at position $(x,y)$ in window $W$. }*/

};


}

#endif
