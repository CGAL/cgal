// Copyright (c) 2009,2014,2015  GeometryFactory Sarl (France)
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
//
// Author(s)     : Laurent RINEAU

//! \file Polyhedron_demo_plugin_helper.h 
#ifndef POLYHEDRON_DEMO_OPERATION_HELPER_H
#define POLYHEDRON_DEMO_OPERATION_HELPER_H

#include <CGAL/license/Three.h>


#include "Scene_item_config.h" //defines SCENE_ITEM_EXPORT

#include <QString>
#include <QStringList>
#include <QMap>

class QAction;
struct QMetaObject;
class QMainWindow;
class QDockWidget;

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
namespace CGAL {
namespace Three {
  /*! \brief Provides convenient functions for a plugin.
   * This class provides convenient functions to manage dock_widgets and to access a certain type of items in the scene.
   * It also provides member variables for a Scene_interface and a QMainWindow.
   */
class SCENE_ITEM_EXPORT Polyhedron_demo_plugin_helper
  : public Polyhedron_demo_plugin_interface
{
public:

  /*! \brief Gets an item of the templated type.
   * \returns The currently selected `SceneType` item
   * \returns the first `SceneType` item found in the scene's list of items if the selected item is not a `SceneType`
   * \returns NULL if there is no `SceneType` in the list.
   */
  template<class SceneType>
  SceneType* getSelectedItem() const{
   int item_id = scene->mainSelectionIndex();
   SceneType* scene_item = qobject_cast<SceneType*>(scene->item(item_id));
   if(!scene_item) {
     // no selected SceneType - if there is only one in list return it, otherwise NULL
     int counter = 0;
     int last_selected = 0;
     for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end && counter < 2; ++i) {
       if(SceneType* tmp = qobject_cast<SceneType*>(scene->item(i))) {
         scene_item = tmp;
         counter++;
         last_selected=i;
       }
     }
     if(counter != 1) { return NULL; }
     scene->setSelectedItem(last_selected);
   }
   return scene_item;
 }

  /*! \brief Adds a dock widget to the interface
   *
   * Adds a dock widget in the left section of the MainWindow. If the slot is already taken, the dock widgets will be tabified.
   */
  void addDockWidget(QDockWidget* dock);
  /*! \brief Automatically connects each action of the plugin to the corresponding slot.
   *
   * \attention Each action named `ActionName` in the plugin's `actions()` list must have a corresponding slot named `on_ActionsName_triggered()`
   * in the plugin.
   */
  void autoConnectActions();
protected:
  //!The reference to the scene
  CGAL::Three::Scene_interface* scene;
  //!The reference to the main window
  QMainWindow* mw;
};
}
}
#endif // POLYHEDRON_DEMO_OPERATION_HELPER_H
