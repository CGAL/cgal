// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// 
//
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_QT_WIDGET_NEF_3_H
#define CGAL_QT_WIDGET_NEF_3_H

#include <CGAL/IO/Qt_widget_OpenGL.h>
#include <CGAL/Nef_3/OGL_helper.h>

namespace CGAL {

template <typename Nef_polyhedron>
class Qt_widget_Nef_3 : public Qt_widget_OpenGL {

public:
  Qt_widget_Nef_3(const Nef_polyhedron& N) : 
    Qt_widget_OpenGL(600,600,0.5) {

    object_ = new CGAL::OGL::Polyhedron();
    CGAL::OGL::Nef3_Converter<Nef_polyhedron>::convert_to_OGLPolyhedron(N,
                                   static_cast<CGAL::OGL::Polyhedron*>(object_));
    resize(window_width,window_height);
    
    main = new QPopupMenu;
    sub1 = new QPopupMenu;
    sub2 = new QPopupMenu;
    sub3 = new QPopupMenu;

    sub1->insertItem("Reset", RESET_CONTROL);
    sub1->insertItem("Rotate", ROTATE);
    sub1->insertItem("Scale", SCALE);
    sub1->insertItem("Translate in XY", TRANSLATE);
    //    sub1->insertItem("Translate in Z", TRANS_Z);
    QObject::connect(sub1, SIGNAL(activated(int)), this, SLOT(slotControlMenu(int)));

    sub2->insertItem("Boundary", CGAL::OGL::SNC_BOUNDARY);
    sub2->insertItem("Skeleton", CGAL::OGL::SNC_SKELETON);
    QObject::connect(sub2, SIGNAL(activated(int)), this, SLOT(slotRenderMenu(int)));
    
    sub3->insertItem("Toggle Axes", CGAL::OGL::SNC_AXES);
    QObject::connect(sub3, SIGNAL(activated(int)), this, SLOT(slotOptionsMenu(int)));

    main->insertItem("&Control", sub1);
    main->insertItem("&Render", sub2);
    main->insertItem("&Options", sub3);
    //    main->insertItem("&Persp/Ortho", this, SLOT(slotPerspective()));
    main->insertItem("&Toggle Fullscreen", this, SLOT(slotFullscreen()));
    main->insertItem("&Quit", qApp, SLOT(quit()));
  }

  ~Qt_widget_Nef_3() {
    delete object_;
  }
};

} //namespace CGAL
#endif // CGAL_QT_WIDGET_NEF_3_H
