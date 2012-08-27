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

#ifndef CGAL_QT_WIDGET_NEF_S2_H
#define CGAL_QT_WIDGET_NEF_S2_H

#include <CGAL/IO/Qt_widget_OpenGL.h>
#include <CGAL/Nef_S2/Sphere_geometry_OGL.h>

namespace CGAL {

template <typename Nef_polyhedron>
class Qt_widget_Nef_S2 : public Qt_widget_OpenGL {

 public:
  Qt_widget_Nef_S2(const typename Nef_polyhedron::Const_decorator& N) : 
    Qt_widget_OpenGL(300,300,1.5) {
    
    object_ = new CGAL::OGL::Unit_sphere(CGAL::OGL::NefS2_to_UnitSphere<Nef_polyhedron>::convert(N));
    resize(window_width, window_height);
  
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
    
    sub2->insertItem("Faces", CGAL::OGL::SM_FACES);
    sub2->insertItem("Skeleton", CGAL::OGL::SM_SKELETON);
    sub2->insertItem("Triangulation",CGAL::OGL::SM_TRIANGULATION);
    QObject::connect(sub2, SIGNAL(activated(int)), this, SLOT(slotRenderMenu(int)));
    
    sub3->insertItem("Toggle Axes", CGAL::OGL::SM_AXES);
    sub3->insertItem("Toggle Unity Cube", CGAL::OGL::SM_CUBE);
    QObject::connect(sub3, SIGNAL(activated(int)), this, SLOT(slotOptionsMenu(int)));
    
    main->insertItem("&Control", sub1);
    main->insertItem("&Render", sub2);
    main->insertItem("&Options", sub3);
    //    main->insertItem("&Persp/Ortho", this, SLOT(slotPerspective()));
    main->insertItem("&Toggle Fullscreen", this, SLOT(slotFullscreen()));
    main->insertItem("&Quit", qApp, SLOT(quit()));
  }

  ~Qt_widget_Nef_S2() {
    delete object_;
  }
};

} //namespace CGAL
#endif // CGAL_QT_WIDGET_NEF_S2_H
