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
