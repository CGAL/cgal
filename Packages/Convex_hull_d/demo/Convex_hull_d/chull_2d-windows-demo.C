// Copyright (c) 2001, 2004  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifdef CGAL_USE_LEDA
#include <CGAL/basic.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/IO/Convex_hull_d_window_stream.h>
#include <iostream>

typedef leda_integer RT;
typedef CGAL::Homogeneous_d<RT> Kernel;
typedef CGAL::Convex_hull_d<Kernel> Convex_hull_d;
typedef Convex_hull_d::Point_d Point;
typedef Convex_hull_d::Simplex_handle Simplex_handle;

int main() {
  CGAL::set_pretty_mode ( std::cerr );
  SETDTHREAD(191);
  leda_string startmess = "input points with left mouse button and ";
  startmess += "exit program with right mouse button!";
  CGAL::Window_stream W; 
  W.set_grid_mode(5);
  W.set_show_coordinates(true);
  W.display(); 
  W.message(startmess);
  double a,b;  // coordinates of a point in the window
  int mouse = W.read_mouse(a,b); 
  // variable to indicate which mouse button was pressed 
  W.del_messages();

  Convex_hull_d T(2);  
  // we are working in the plane

  std::ofstream To("ch2-demo.log");
  CGAL_assertion(To);
  while (mouse != MOUSE_BUTTON(3)) {
    // while mouse click is not the right button
    RT ia(a), ib(b);

    To << a << "," << b << std::endl; To.flush(); 
    Point x(ia,ib); 

    T.insert(x); 
    T.is_valid(true);
    W.clear(); 
    CGAL::d2_show(T,W); 
    
    mouse = W.read_mouse(a,b);  
    // read the window coordinates into a and b
  }
  return 0;
}

#else
#include <iostream>

int main()
{ 
  std::cout << "this program requires LEDA" << std::endl;
  return 0;
}

#endif // CGAL_USE_LEDA


