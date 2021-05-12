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

#ifndef SCENE_ZOOMABLE_ITEM_INTERFACE_H
#define SCENE_ZOOMABLE_ITEM_INTERFACE_H

#include <CGAL/license/Three.h>
#include <QtPlugin>
#include <QPoint>
namespace CGAL
{
namespace Three {
  class Viewer_interface;


//! This class provides a function to move the camera orthogonaly to the wanted region
class Scene_zoomable_item_interface {
public:
  virtual ~Scene_zoomable_item_interface(){}
 //! Move the camera orthogonaly to the region defined by `point`
 virtual void zoomToPosition(const QPoint& point, CGAL::Three::Viewer_interface*)const = 0;
};
}
}

Q_DECLARE_INTERFACE(CGAL::Three::Scene_zoomable_item_interface, "com.geometryfactory.PolyhedronDemo.ZoomInterface/1.0")
#endif // SCENE_ZOOMABLE_ITEM_INTERFACE_H
