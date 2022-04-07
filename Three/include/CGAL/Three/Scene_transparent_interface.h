// Copyright (c) 2009,2010,2012,2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime GIMENO
#ifndef SCENE_TRANSPARENT_INTERFACE_H
#define SCENE_TRANSPARENT_INTERFACE_H
#include <CGAL/license/Three.h>
#include <QtPlugin>
namespace CGAL
{
namespace Three {
  class Viewer_interface;


//! Base class to allow an item to draw transparent faces.
class Scene_transparent_interface {
public:
  virtual ~Scene_transparent_interface(){}
 //! Draw transparent faces. It is the last drawing call in the Scene, so all
 //! the items are already in place and the blending
 //! takes them all, their points and their edges into account, whatever the
 //! position of the transparent item is in the scene's entries.
 //!
 virtual void drawTransparent(CGAL::Three::Viewer_interface*)const = 0;
};
}
}
Q_DECLARE_INTERFACE(CGAL::Three::Scene_transparent_interface, "com.geometryfactory.PolyhedronDemo.TransparentInterface/1.0")
#endif // SCENE_TRANSPARENT_INTERFACE_H
