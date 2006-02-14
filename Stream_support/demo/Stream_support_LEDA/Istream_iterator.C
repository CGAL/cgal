// Copyright (c) 2001, 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>


#ifndef CGAL_USE_LEDA
#include <iostream>
int main(){ std::cout << "This demo needs LEDA" << std::endl; return 0;}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Istream_iterator.h>
#include <CGAL/IO/Window_stream.h>
#include <iostream>
#include <algorithm>

typedef CGAL::Cartesian<double>::Point_2                   Point;
typedef CGAL::Istream_iterator<Point, CGAL::Window_stream> Iterator;

#ifdef CGAL_USE_CGAL_WINDOW
#define leda_window CGAL::window
#define leda_green  CGAL::green
#endif

void init_window( leda_window& W) {
    CGAL::cgalize( W);
    W.set_fg_color( leda_green);
    W.display();
    W.init(-1.0, 1.0, -1.0);
}

int main () {
    CGAL::Window_stream window( 512, 512);
    init_window(window);
    std::copy( Iterator(window), Iterator(),
               std::ostream_iterator<Point>(std::cout,"\n"));
    return 0;
}

#endif
