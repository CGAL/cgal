// Copyright (c) 2018 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_SURFACE_MESH_H
#define CGAL_DRAW_SURFACE_MESH_H

#ifdef DOXYGEN_RUNNING

/*!
\ingroup PkgDrawNefS2

Open a new window and draw `asm`, an instance of the `CGAL::Nef_polyhedron_S2` class. The function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam SM an instance of the `CGAL::Nef_polyhedron_S2` class.
\param asm the surface mesh to draw.

*/
template<class Nef_polyhedron>
void draw(const Nef_polyhedron& asm);

/*!
\ingroup PkgDrawNefS2

Open a new window and draw a set of instances of the `CGAL::Nef_polyhedron_S2` class. The function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam SM an instance of the `CGAL::Nef_polyhedron_S2` class.
\param asm the surface mesh to draw.

*/
template<class iterator>
void draw(iterator start, iterator end);
#else // DOXYGEN_RUNNING
  
#include <CGAL/license/Nef_S2.h>

#include <CGAL/Nef_S2/Qt_widget_Nef_S2.h>

namespace CGAL
{


    template<class Nef_polyhedron,class iterator>
        void draw(iterator start, iterator end)
        { 
#if defined(CGAL_TEST_SUITE)
            bool cgal_test_suite=true;
#else
            bool cgal_test_suite=false;
#endif

            if (!cgal_test_suite)
            {
                int argc=1;
                const char* argv[2]={"Nef_S2 viewer","\0"};
                QApplication app(argc,const_cast<char**>(argv));
                CGAL::Qt_widget_Nef_S2 mainwindow;
                unsigned int i=0;
                for (iterator it=start;it!=end;it++,i++) {
                    char lbl[128];
                    sprintf(lbl,"NefS2 %d",i);
                    mainwindow.addObject<Nef_polyhedron>(*it,lbl); 
                }
                mainwindow.show();
                app.exec();
            }

        }

    template<class Nef_polyhedron>
        void draw(const Nef_polyhedron& S)
        { 
#if defined(CGAL_TEST_SUITE)
            bool cgal_test_suite=true;
#else
            bool cgal_test_suite=false;
#endif

            if (!cgal_test_suite)
            {
                int argc=1;
                const char* argv[2]={"Nef_S2 viewer","\0"};
                QApplication app(argc,const_cast<char**>(argv));
                CGAL::Qt_widget_Nef_S2 mainwindow;
                mainwindow.addObject<Nef_polyhedron>(S); 
                mainwindow.show();
                app.exec();
            }

        }
} // End namespace CGAL


#endif // DOXYGEN_RUNNING

#endif // CGAL_DRAW_SURFACE_MESH_H
