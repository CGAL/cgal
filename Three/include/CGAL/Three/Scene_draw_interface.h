// Copyright (c) 2009,2010,2012,2015  GeometryFactory Sarl (France)
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
// Author(s)     : Laurent RINEAU

#ifndef SCENE_DRAW_INTERFACE_H
#define SCENE_DRAW_INTERFACE_H

class QKeyEvent;
namespace CGAL
{
namespace Three {
  class Viewer_interface;


//! Base class to interact with the scene from the viewer.
class Scene_draw_interface {
public:
  virtual ~Scene_draw_interface(){}
  virtual void initializeGL() = 0;
  virtual void draw() = 0;
  virtual void draw(CGAL::Three::Viewer_interface*) { draw(); };
  virtual void drawWithNames() = 0;
  virtual void drawWithNames(CGAL::Three::Viewer_interface*) { drawWithNames(); }
  virtual bool keyPressEvent(QKeyEvent* e) = 0;
  virtual float get_bbox_length() const = 0;
};
}
}
#endif // SCENE_DRAW_INTERFACE_H;
