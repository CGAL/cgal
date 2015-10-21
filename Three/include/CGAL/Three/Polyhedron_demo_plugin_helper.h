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
//
//
// Author(s)     : Laurent RINEAU

//! \file Polyhedron_demo_plugin_helper.h 
#ifndef POLYHEDRON_DEMO_OPERATION_HELPER_H
#define POLYHEDRON_DEMO_OPERATION_HELPER_H

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
  /*!
   * This class provides a base for creating a new plugin.
   */
class SCENE_ITEM_EXPORT Polyhedron_demo_plugin_helper
  : public Polyhedron_demo_plugin_interface
{
public:
  //! get action object from its name
  static QAction* getActionFromMainWindow(QMainWindow*, QString action_name);
  
  //! Init plugin
  virtual void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface);
  
  //! Get list of actions supported by this plugin
  virtual QStringList actionsNames() const;
  //!List of the actions of the plugin
  virtual QList<QAction*> actions() const;

  //! To get a selected item with the type of SceneType
  template<class SceneType>
  SceneType* get_selected_item() const {
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

  //!To add a dock widget to the interface
  void add_dock_widget(QDockWidget* dock);

  //! Auto-connects actions to slots. Called by init().
  void autoConnectActions();
  
protected:
  //!The list of actions
  QMap<QString, QAction*> actions_map;
  //!The reference to the scene
  CGAL::Three::Scene_interface* scene;
  //!The reference to the main window
  QMainWindow* mw;
};
}
}
#endif // POLYHEDRON_DEMO_OPERATION_HELPER_H
